#ifndef __FWDPY_QTRAIT_DETAILS_HPP__
#define __FWDPY_QTRAIT_DETAILS_HPP__

#include "types.hpp"
#include <fwdpp/experimental/sample_diploid.hpp>
#include <fwdpp/extensions/regions.hpp>

//Namespace pollution due to Cython/struct/namespace bug:
struct qtrait_sample_info_t
{
  KTfwd::sep_sample_t genotypes;
  std::vector<std::pair<double, double> > sh;
  template<typename T1,typename T2>
  qtrait_sample_info_t( T1 && a, T2 && b ) noexcept :
					    genotypes(std::forward<T1>(a)),
						sh(std::forward<T2>(b))
  {
  }
  qtrait_sample_info_t() noexcept : genotypes(KTfwd::sep_sample_t()),
				    sh(std::vector<std::pair<double,double>>())
  {
  }
};

namespace fwdpy {
  namespace qtrait {
    template<typename rules>
    void
    qtrait_sim_details_t( unsigned long seed,
			  fwdpy::singlepop_t * pop,
			  const unsigned * Nvector,
			  const size_t Nvector_len,
			  const double neutral,
			  const double selected,
			  const double recrate,
			  const double f,
			  const double sigmaE,
			  const double optimum ,
			  KTfwd::extensions::discrete_mut_model && __m,
			  KTfwd::extensions::discrete_rec_model && __recmap,
			  rules && __model_rules)
    {
      gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
      gsl_rng_set(rng,seed);
      const unsigned simlen = Nvector_len;
      const double mu_tot = neutral + selected;

      KTfwd::extensions::discrete_mut_model m(std::move(__m));
      KTfwd::extensions::discrete_rec_model recmap(std::move(__recmap));
      rules model_rules(std::move(__model_rules));
      const auto recpos = KTfwd::extensions::bind_drm(recmap,pop->gametes,pop->mutations,
						      rng,recrate);

      //We use an empty fitness fxn here b/c the rules policies keep track of it separately.
      const auto ff = []( const fwdpy::singlepop_t::diploid_t &,
			  const fwdpy::singlepop_t::gcont_t &,
			  const fwdpy::singlepop_t::mcont_t ) noexcept { return 0.; };

      for( unsigned g = 0 ; g < simlen ; ++g, ++pop->generation )
	{
	  const unsigned nextN = *(Nvector+g);
	  KTfwd::experimental::sample_diploid(rng,
					      pop->gametes,
					      pop->diploids,
					      pop->mutations,
					      pop->mcounts,
					      pop->N,
					      nextN,
					      mu_tot,
					      KTfwd::extensions::bind_dmm(m,pop->mutations,pop->mut_lookup,rng,neutral,selected,pop->generation),
					      recpos,
					      ff,
					      pop->neutral,pop->selected,
					      f,
					      model_rules,
					      KTfwd::remove_nothing());
	  KTfwd::update_mutations_n(pop->mutations,pop->fixations,pop->fixation_times,pop->mut_lookup,pop->mcounts,pop->generation,2*nextN);
	  assert(KTfwd::check_sum(pop->gametes,2*nextN));
	}
      gsl_rng_free(rng);
      //Update population's size variable to be the current pop size
      pop->N = pop->diploids.size();
    }


    inline std::vector<std::pair<double, double> >
    get_sh(const fwdpy::singlepop_t::mcont_t & mutations,
	   const KTfwd::sep_sample_t & genotypes)
    {
      std::vector<std::pair<double, double> > rv;
      for( const auto & site : genotypes.second )
	{
	  auto itr = std::find_if(mutations.begin(),mutations.end(),
				  [&site](const fwdpy::singlepop_t::mutation_t & mut) {
				    return site.first == mut.pos;
				  });
	  rv.emplace_back(itr->s,itr->h);
	}
      return rv;
    }
    
    //Return type means that should be run using std::async
    template<typename rules>
    std::vector< std::pair<unsigned,qtrait_sample_info_t > > //rv is pair<generation,sample>
    qtrait_sim_details_samples_t( unsigned long seed,
				  fwdpy::singlepop_t * pop,
				  const unsigned * Nvector,
				  const size_t Nvector_len,
				  const double neutral,
				  const double selected,
				  const double recrate,
				  const double f,
				  const double sigmaE,
				  const double optimum ,
				  const int trackSamples,  //do we want to track the trajectories of all mutations and how often?
				  const unsigned nsam, //Sample size
				  KTfwd::extensions::discrete_mut_model && __m,
				  KTfwd::extensions::discrete_rec_model && __recmap,
				  rules && __model_rules)
    {
      gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
      gsl_rng_set(rng,seed);
      const unsigned simlen = Nvector_len;
      const double mu_tot = neutral + selected;

      KTfwd::extensions::discrete_mut_model m(std::move(__m));
      KTfwd::extensions::discrete_rec_model recmap(std::move(__recmap));
      rules model_rules(std::move(__model_rules));
      const auto recpos = KTfwd::extensions::bind_drm(recmap,pop->gametes,pop->mutations,
						      rng,recrate);

      //We use an empty fitness fxn here b/c the rules policies keep track of it separately.
      const auto ff = []( const fwdpy::singlepop_t::diploid_t &,
			  const fwdpy::singlepop_t::gcont_t &,
			  const fwdpy::singlepop_t::mcont_t ) noexcept { return 0.; };
      std::vector< std::pair<unsigned,qtrait_sample_info_t > > rv;
      for( unsigned g = 0 ; g < simlen ; ++g, ++pop->generation )
	{
	  if(trackSamples&&pop->generation&&pop->generation%trackSamples==0)
	    {
	      //Get a sample
	      auto x = KTfwd::ms_sample_separate(rng,pop->mutations,pop->gametes,pop->diploids,nsam);
	      //get effect size and dominance for each mutation in sample
	      rv.emplace_back(pop->generation,qtrait_sample_info_t(std::move(x),get_sh(pop->mutations,x)));
	    }
	  const unsigned nextN = *(Nvector+g);
	  KTfwd::experimental::sample_diploid(rng,
					      pop->gametes,
					      pop->diploids,
					      pop->mutations,
					      pop->mcounts,
					      pop->N,
					      nextN,
					      mu_tot,
					      KTfwd::extensions::bind_dmm(m,pop->mutations,pop->mut_lookup,rng,neutral,selected,pop->generation),
					      recpos,
					      ff,
					      pop->neutral,pop->selected,
					      f,
					      model_rules,
					      KTfwd::remove_nothing());
	  KTfwd::update_mutations_n(pop->mutations,pop->fixations,pop->fixation_times,pop->mut_lookup,pop->mcounts,pop->generation,2*nextN);
	  assert(KTfwd::check_sum(pop->gametes,2*nextN));
	}

      //make sure we update in the last generation if needed
      if(trackSamples&&pop->generation&&pop->generation%trackSamples==0)
	{
	  //Get a sample
	  auto x = KTfwd::ms_sample_separate(rng,pop->mutations,pop->gametes,pop->diploids,nsam);
	  //get effect size and dominance for each mutation in sample
	  rv.emplace_back(pop->generation,qtrait_sample_info_t(std::move(x),get_sh(pop->mutations,x)));
	}
      gsl_rng_free(rng);
      //Update population's size variable to be the current pop size
      pop->N = pop->diploids.size();
      return rv;
    }
  }
} //namespace

#endif
