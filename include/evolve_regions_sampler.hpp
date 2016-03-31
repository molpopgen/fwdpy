#ifndef FWDPY_EVOLVE_REGIONS_SAMPLER_HPP
#define FWDPY_EVOLVE_REGIONS_SAMPLER_HPP

#include <future>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <fwdpp/diploid.hh>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/sugar/sampling.hpp>

#include "types.hpp"
#include "reserve.hpp"
#include "internal_region_manager.hpp"
#include "sampler_pop_properties.hpp"
#include "sampler_sample_n.hpp"
#include "sampler_selected_mut_tracker.hpp"

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
  template<typename sampler,class... Args>
  inline typename sampler::final_t
  evolve_regions_sampler_details(singlepop_t * pop,
				 const unsigned long seed,
				 const unsigned * Nvector,
				 const size_t Nvector_len,
				 const double neutral,
				 const double selected,
				 const double recrate,
				 const double f,
				 const char * fitness,
				 const int interval,
				 KTfwd::extensions::discrete_mut_model && __m,
				 KTfwd::extensions::discrete_rec_model && __recmap,
				 Args&&... args) 
  {
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

    //create the sampler
    sampler s(std::forward<Args>(args)...);
    
    for( size_t g = 0 ; g < simlen ; ++g, ++pop->generation )
      {
	const unsigned nextN = 	*(Nvector+g);
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
	if (interval && pop->generation &&(pop->generation+1)%interval==0.)
	  {
	    s(pop,rng,pop->generation+1);
	  }
	KTfwd::update_mutations(pop->mutations,pop->fixations,pop->fixation_times,pop->mut_lookup,pop->mcounts,pop->generation,2*nextN);
	assert(KTfwd::check_sum(pop->gametes,2*nextN));
      }
    //Update population's size variable to be the current pop size
    pop->N = unsigned(pop->diploids.size());
    //cleanup
    gsl_rng_free(rng);
    return s.final();
  }

  template<typename sampler,class... Args>
  inline
  typename std::enable_if<!std::is_void<typename sampler::final_t>::value,
			  std::vector<typename sampler::final_t> >::type
  evolve_regions_async_wrapper( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
				const unsigned * Nvector,
				const size_t Nvector_len,
				const double mu_neutral,
				const double mu_selected,
				const double littler,
				const double f,
				const int sample,
				const internal::region_manager * rm,
				const char * fitness,
				Args&&... args)
  {
    using future_t = std::future<typename sampler::final_t>;
    std::vector<future_t> futures;
    for(std::size_t i=0;i<pops->size();++i)
      {
	futures.emplace_back( std::async(std::launch::async,
					 evolve_regions_sampler_details<sampler,Args&&...>,
					 pops->operator[](i).get(),gsl_rng_get(rng->get()),Nvector,Nvector_len,
					 mu_neutral,mu_selected,littler,f,fitness,sample,
					 std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
					 std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)),
					 std::forward<Args...>(args)...
					 )
			      );	
      }
    std::vector<typename sampler::final_t> rv(futures.size());
    for(std::size_t i=0;i<futures.size();++i ) rv[i]=futures[i].get();
    return rv;
  }

  template<typename sampler,class... Args>
  inline
  typename std::enable_if<std::is_void<typename sampler::final_t>::value,
			  void>::type
  evolve_regions_async_wrapper( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
				const unsigned * Nvector,
				const size_t Nvector_len,
				const double mu_neutral,
				const double mu_selected,
				const double littler,
				const double f,
				const internal::region_manager * rm,
				const char * fitness,
				Args&&... args)
  {
    using future_t = std::future<typename sampler::final_t>;
    std::vector<future_t> futures;
    for(std::size_t i=0;i<pops->size();++i)
      {
	futures.emplace_back( std::async(std::launch::async,
					 evolve_regions_sampler_details<sampler,Args&&...>,
					 pops->operator[](i).get(),gsl_rng_get(rng->get()),Nvector,Nvector_len,
					 mu_neutral,mu_selected,littler,f,fitness,0,//0 will mean not to sample
					 std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
					 std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)),
					 std::forward<Args...>(args)...
					 )
			      );	
      }
    for(std::size_t i=0;i<futures.size();++i) futures[i].get();
  }
  

  //Prototype for 'default' function that doesn't sample
  void evolve_regions_no_sampling_async(GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
					const unsigned * Nvector,
					const size_t Nvector_length,
					const double mu_neutral,
					const double mu_selected,
					const double littler,
					const double f,
					const internal::region_manager * rm,
					const char * fitness);
  
  //Prototypes for functions using samplers
  std::vector<sample_n::final_t>
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
  std::vector<selected_mut_tracker::final_t>
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
