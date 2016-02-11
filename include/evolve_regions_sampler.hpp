#ifndef FWDPY_EVOLVE_REGIONS_SAMPLER_HPP
#define FWDPY_EVOLVE_REGIONS_SAMPLER_HPP

#include <vector>
#include <type_traits>
#include <types.hpp>
#include <reserve.hpp>
#include <internal/internal.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/sugar/sampling.hpp>

//Namespace pollution!!
struct detailed_deme_sample
{
  KTfwd::sep_sample_t genotypes;
  std::vector<std::pair<double,double> > sh;
  template<typename T1,typename T2>
  detailed_deme_sample(T1 && t1, T2 && t2) : genotypes(std::forward<T1>(t1)),
					     sh(std::forward<T2>(t2))
  {
  }
};

struct selected_mut_data
{
  unsigned generation;
  double pos,freq,esize;
  explicit selected_mut_data(unsigned g,double p,double f,double e) :
    generation(g),pos(p),freq(f),esize(e)
  {
  }
};

namespace fwdpy
{
  /*
    This is the generic function.

    The type "sampler" must define an operator() taking a const singlepop_t * and a gsl_rng *, and unsigned as arguments.

    The unsigned represents the generation, and the intent is that it is in the return value of sampler(...) somewhere.
  
    Other arguments may be bound by the caller, etc.

    The object s will be applied every "interval" generations.

    If result_t is the return type of s(const singlepop_t *, gsl_rng *), then the return type
    of this function is vector<pair<unsigned,result_t> >, where the unsigned values refer
    to the generation in which a specific result_t was obtained.

    Design issues to work out:

    1. The result_t is an object is considered below.  What happens when it is 
    a vector of objects?

    For example, a sample is a single object.  The selected mutation frequencies are a vector
    of data on mutations.

    Perhaps that is just ok? It isn't super-friendly, and would need conversion to pandas.DataFrame,
    and I guess we'd have to supply those functions.

    The problem is one of space.  Doing frequency trajectories this way will be much more RAM-intensive
    Than the current implementation.  The "qtrait stats" is probably fine, with some overhead from strings.
  */
  template<typename sampler>
  inline auto
  evolve_regions_sampler_details(singlepop_t * pop,
				 const unsigned long seed,
				 const unsigned * Nvector,
				 const size_t Nvector_len,
				 const double neutral,
				 const double selected,
				 const double recrate,
				 const double f,
				 const char * fitness,
				 const sampler & s,
				 const int interval,
				 KTfwd::extensions::discrete_mut_model && __m,
				 KTfwd::extensions::discrete_rec_model && __recmap) -> std::vector< typename std::result_of<sampler(const singlepop_t *, gsl_rng *, const unsigned)>::type>
  {
    std::vector< typename std::result_of<sampler(const singlepop_t *, gsl_rng *, const unsigned)>::type> rv;

    const size_t simlen = Nvector_len;
    auto x = std::max_element(Nvector,Nvector+Nvector_len);
    assert(x!=Nvector+Nvector_len);
    reserve_space(pop->gametes,pop->mutations,*x,neutral+selected);
    const double mu_tot = neutral + selected;
    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
    KTfwd::extensions::discrete_mut_model m(std::move(__m));
    KTfwd::extensions::discrete_rec_model recmap(std::move(__recmap));
    //Recombination policy: more complex than the standard case...
    const auto recpos = KTfwd::extensions::bind_drm(recmap,pop->gametes,pop->mutations,
						    rng,recrate);

    /*
      The fitness model.

      Normally, we'd declare dipfit as "auto", but that won't work here b/c there is the chance
      that we have to re-assign it using an additive model based on input from calling environment.

      The std::bind signature has a different type for the two models, and thus we must coerce it to
      the proper function signature, which is a member typedef provided by the fwdpp sugar type
      from which fwdpy::singlepop_t inherits
    */
    fwdpy::singlepop_t::fitness_t dipfit = std::bind(KTfwd::multiplicative_diploid(),
						     std::placeholders::_1,
						     std::placeholders::_2,
						     std::placeholders::_3,
						     2.);

    if( std::string(fitness) == "additive" )
      {
	dipfit = std::bind(KTfwd::additive_diploid(),
			   std::placeholders::_1,
			   std::placeholders::_2,
			   std::placeholders::_3,
			   2.);
      }

    for( size_t g = 0 ; g < simlen ; ++g, ++pop->generation )
      {
	const unsigned nextN = 	*(Nvector+g);
	if (interval && pop->generation &&pop->generation%interval==0.)
	  {
	    rv.emplace_back(s(pop,rng,pop->generation));
	  }
	KTfwd::sample_diploid(rng,
			      pop->gametes,
			      pop->diploids,
			      pop->mutations,
			      pop->mcounts,
			      pop->N,
			      nextN,
			      mu_tot,
			      KTfwd::extensions::bind_dmm(m,pop->mutations,pop->mut_lookup,rng,neutral,selected,pop->generation),
			      recpos,
			      dipfit,
			      pop->neutral,pop->selected,
			      f);
	pop->N=nextN;
	KTfwd::update_mutations(pop->mutations,pop->fixations,pop->fixation_times,pop->mut_lookup,pop->mcounts,pop->generation,2*nextN);
	assert(KTfwd::check_sum(pop->gametes,2*nextN));
      }
    if (interval && pop->generation &&pop->generation%interval==0.)
      {
	rv.emplace_back(s(pop,rng,pop->generation));
      }
    //Update population's size variable to be the current pop size
    pop->N = unsigned(pop->diploids.size());
    //cleanup
    gsl_rng_free(rng);
    //restore data
    return rv;
  }

  //These are specific samplers
  
  struct sample_n //take a sample of size n from a population
  {
    using result_type = std::pair<unsigned,detailed_deme_sample>;
    using bound_t = std::function<result_type(const singlepop_t * ,gsl_rng * r,const unsigned)>;
    inline result_type operator()(const singlepop_t * pop,gsl_rng * r,
				  const unsigned generation, const unsigned nsam) const
    {
      auto s = KTfwd::sample_separate(r,*pop,nsam,true);
      std::vector< std::pair<double,double> > sh;
      for( const auto & i : s.second )
	{
	  auto itr = std::find_if(pop->mutations.begin(),pop->mutations.end(),[&i](const singlepop_t::mutation_t & m) noexcept
				  {
				    return m.pos == i.first;
				  });
	  sh.emplace_back(itr->s,itr->h);
	}
      return std::make_pair(generation,detailed_deme_sample(std::move(s),std::move(sh)));
    }
    static inline bound_t make_bound_t(unsigned nsam) noexcept
    {
      return std::bind(sample_n(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,nsam);
    }
  };

  struct get_selected_mut_data //record info on selected mutations in population, including fixations
  {
    using result_type = std::vector<selected_mut_data>;
    using bound_t = std::function<result_type(const singlepop_t * ,gsl_rng * r,const unsigned)>;
    static inline bound_t make_bound_t() noexcept
    {
      return std::bind(get_selected_mut_data(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
    }
    inline result_type operator()(const singlepop_t * pop,gsl_rng * r,
				  const unsigned generation) const
    {
      result_type rv;
      double twoN=2.0*double(pop->diploids.size());
      for(std::size_t i=0;i<pop->mcounts.size();++i)
	{
	  if(pop->mcounts[i])
	    {
	      rv.emplace_back(generation,pop->mutations[i].pos,
			      double(pop->mcounts[i])/twoN,pop->mutations[i].s);
	    }
	}
      return rv;
    }
  };
  //Prototypes for functions using samplers

  std::vector<std::vector<sample_n::result_type> >
  evolve_regions_sample_async( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
			       const unsigned * Nvector,
			       const size_t Nvector_length,
			       const double mu_neutral,
			       const double mu_selected,
			       const double littler,
			       const double f,
			       const int sample,
			       const unsigned nsam,
			       const internal::region_manager * rm,
			       const char * fitness);

  //This function will use moves to collapse the ugly return type to something nicer
  std::vector<get_selected_mut_data::result_type>
  evolve_regions_track_async( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
			       const unsigned * Nvector,
			       const size_t Nvector_length,
			       const double mu_neutral,
			       const double mu_selected,
			       const double littler,
			       const double f,
			       const int sample,
			       const internal::region_manager * rm,
			       const char * fitness);

} //ns fwdpy
#endif
