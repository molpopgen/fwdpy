#ifndef FWDP_QTRAIT_EVOLVE_QTRAIT_SAMPLER_HPP
#define FWDP_QTRAIT_EVOLVE_QTRAIT_SAMPLER_HPP

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
  namespace qtrait
  {
    template<typename sampler,typename rules_t,class... Args>
    inline typename sampler::final_t
    evolve_qtrait_sampler_details(singlepop_t * pop,
				  const unsigned long seed,
				  const unsigned * Nvector,
				  const size_t Nvector_len,
				  const double neutral,
				  const double selected,
				  const double recrate,
				  const double f,
				  const double sigmaE,
				  const double optimum,
				  const int interval,
				  KTfwd::extensions::discrete_mut_model && __m,
				  KTfwd::extensions::discrete_rec_model && __recmap,
				  rules_t && rules,
				  Args&&... args)
    {
      gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
      gsl_rng_set(rng,seed);
      const unsigned simlen = unsigned(Nvector_len);
      const double mu_tot = neutral + selected;

      KTfwd::extensions::discrete_mut_model m(std::move(__m));
      KTfwd::extensions::discrete_rec_model recmap(std::move(__recmap));
      rules_t model_rules(std::forward<rules_t>(rules));
      const auto recpos = KTfwd::extensions::bind_drm(recmap,pop->gametes,pop->mutations,
						      rng,recrate);
      //create the sampler
      sampler s(std::forward<Args>(args)...);
      //We use an empty fitness fxn here b/c the rules policies keep track of it separately.
      const auto ff = []( const fwdpy::singlepop_t::diploid_t &,
			  const fwdpy::singlepop_t::gcont_t &,
			  const fwdpy::singlepop_t::mcont_t ) noexcept { return 0.; };

      for( unsigned g = 0 ; g < simlen ; ++g, ++pop->generation )
	{
	  const unsigned nextN = *(Nvector+g);
	  if (interval && pop->generation &&pop->generation%interval==0.)
	    {
	      s(pop,rng,pop->generation);
	    }
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
					      rules,
					      KTfwd::remove_neutral());

	  KTfwd::update_mutations_n(pop->mutations,pop->fixations,pop->fixation_times,pop->mut_lookup,pop->mcounts,pop->generation,2*nextN);
	  assert(KTfwd::check_sum(pop->gametes,2*nextN));
	}
      if (interval && pop->generation &&pop->generation%interval==0.)
	{
	  s(pop,rng,pop->generation);
	}
      gsl_rng_free(rng);
      //Update population's size variable to be the current pop size
      pop->N = KTfwd::uint_t(pop->diploids.size());
      return s.final();
    }
    
    template<typename sampler,typename rules_t,class... Args>
    inline
    typename std::enable_if<!std::is_void<typename sampler::final_t>::value,
			    std::vector<typename sampler::final_t> >::type
    evolve_qtrait_async_wrapper( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
				 const unsigned * Nvector,
				 const size_t Nvector_len,
				 const double mu_neutral,
				 const double mu_selected,
				 const double littler,
				 const double f,
				 const double sigmaE,
				 const double optimum,
				 const int sample,
				 const internal::region_manager * rm,
				 const rules_t & rules,
				 Args&&... args)
    {
      using future_t = std::future<typename sampler::final_t>;
      std::vector<future_t> futures;
      for(std::size_t i=0;i<pops->size();++i)
	{
	  rules_t rules_thread(rules);
	  futures.emplace_back( async(std::launch::async,
				      evolve_qtrait_sampler_details<sampler,rules_t,Args&&...>,
				      pops->operator[](i).get(),gsl_rng_get(rng->get()),Nvector,Nvector_len,
				      mu_neutral,mu_selected,littler,f,sigmaE,optimum,sample,
				      std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
				      std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)),
				      std::move(rules_thread),
				      std::forward<Args...>(args)...
				      )
				);	
	}
      std::vector<typename sampler::final_t> rv(futures.size());
      for(std::size_t i=0;i<futures.size();++i ) rv[i]=futures[i].get();
      return rv;
    }

    template<typename sampler,typename rules_t,class... Args>
    inline
    typename std::enable_if<std::is_void<typename sampler::final_t>::value,
			    void>::type
    evolve_qtrait_async_wrapper( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
				 const unsigned * Nvector,
				 const size_t Nvector_len,
				 const double mu_neutral,
				 const double mu_selected,
				 const double littler,
				 const double f,
				 const double sigmaE,
				 const double optimum,
				 const internal::region_manager * rm,
				 const rules_t & rules,
				 Args&&... args)
    {
      using future_t = std::future<typename sampler::final_t>;
      std::vector<future_t> futures;
      for(std::size_t i=0;i<pops->size();++i)
	{
	  rules_t rules_thread(rules);
	  futures.emplace_back( async(std::launch::async,
				      evolve_qtrait_sampler_details<sampler,rules_t,Args&&...>,
				      pops->operator[](i).get(),gsl_rng_get(rng->get()),Nvector,Nvector_len,
				      mu_neutral,mu_selected,littler,f,sigmaE,optimum,0,
				      std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
				      std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)),
				      std::move(rules_thread),
				      std::forward<Args...>(args)...
				      )
				);	
	}
      for(std::size_t i=0;i<futures.size();++i ) futures[i].get();
    }

    //No sampling
    void evolve_qtrait_no_sampling_async( GSLrng_t * rng,
					  std::vector<std::shared_ptr<singlepop_t> > * pops,
					  const unsigned * Nvector,
					  const size_t Nvector_length,
					  const double mu_neutral,
					  const double mu_selected,
					  const double littler,
					  const double f,
					  const double sigmaE,
					  const double optimum,
					  const double VS,
					  const internal::region_manager * rm);
    
    //Take samples over time
    std::vector<sample_n::final_t>
    evolve_qtrait_sample_async( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
				const unsigned * Nvector,
				const size_t Nvector_length,
				const double mu_neutral,
				const double mu_selected,
				const double littler,
				const double f,
				const double sigmaE,
				const double optimum,
				const double VS,
				const int sample,
				const unsigned nsam,
				const internal::region_manager * rm);

    //Track VG, etc.
    std::vector<pop_properties::final_t>
    evolve_qtrait_popstats_async( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
				  const unsigned * Nvector,
				  const size_t Nvector_length,
				  const double mu_neutral,
				  const double mu_selected,
				  const double littler,
				  const double f,
				  const double sigmaE,
				  const double optimum,
				  const double VS,
				  const int sample,
				  const internal::region_manager * rm);
    
    //track mutation frequencies
    std::vector<selected_mut_tracker::final_t>
    evolve_qtrait_track_async( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
			       const unsigned * Nvector,
			       const size_t Nvector_length,
			       const double mu_neutral,
			       const double mu_selected,
			       const double littler,
			       const double f,
			       const double sigmaE,
			       const double optimum,
			       const double VS,
			       const int track,
			       const internal::region_manager * rm);

    void evolve_gbr_no_sampling_async( GSLrng_t * rng,
				       std::vector<std::shared_ptr<singlepop_t> > * pops,
				       const unsigned * Nvector,
				       const size_t Nvector_length,
				       const double mu_neutral,
				       const double mu_selected,
				       const double littler,
				       const double f,
				       const double sigmaE,
				       const double optimum,
				       const double VS,
				       const internal::region_manager * rm);
    
    //Take samples over time
    std::vector<sample_n::final_t>
    evolve_gbr_sample_async( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
			     const unsigned * Nvector,
			     const size_t Nvector_length,
			     const double mu_neutral,
			     const double mu_selected,
			     const double littler,
			     const double f,
			     const double sigmaE,
			     const double optimum,
			     const double VS,
			     const int sample,
			     const unsigned nsam,
			     const internal::region_manager * rm);

    //Track VG, etc.
    std::vector<pop_properties::final_t>
    evolve_gbr_popstats_async( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
			       const unsigned * Nvector,
			       const size_t Nvector_length,
			       const double mu_neutral,
			       const double mu_selected,
			       const double littler,
			       const double f,
			       const double sigmaE,
			       const double optimum,
			       const double VS,
			       const int sample,
			       const internal::region_manager * rm);

    //track mutation frequencies
    std::vector<selected_mut_tracker::final_t>
    evolve_gbr_track_async( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
			       const unsigned * Nvector,
			       const size_t Nvector_length,
			       const double mu_neutral,
			       const double mu_selected,
			       const double littler,
			       const double f,
			       const double sigmaE,
			       const double optimum,
			       const double VS,
			       const int track,
			       const internal::region_manager * rm);
  }
}


#endif
