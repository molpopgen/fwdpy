/*
  Models of quantitative traits
  Trait values are additive over 1, 1+hs, 1+2s, where s is a Gaussian deviate

  The infinitely-many sites stuff is an Cython/fwdpp-based re-implementation of the
  code used to generate preliminary data for R01GM115564.
*/

#include <future>
#include <thread>
#include <algorithm>
#include <memory>
#include <limits>
#include <vector>
#include <utility>
#include <fwdpp/diploid.hh>
#include <fwdpp/experimental/sample_diploid.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/popgenmut.hpp>

#include <types.hpp>
#include <internal_region_manager.hpp>
#include <qtrait_evolve_rules.hpp>
#include <qtrait_details.hpp>
#include <qtrait_evolve.hpp>
#include <no_sampling.hpp>
#include <gsl/gsl_statistics_double.h>

using namespace std;

using poptype = fwdpy::singlepop_t;

namespace fwdpy
{
  namespace qtrait
  {
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
					  const internal::region_manager * rm)
    {
      qtrait_model_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      evolve_qtrait_async_wrapper<no_sampling,qtrait_model_rules>(rng,pops,Nvector,Nvector_length,
								  mu_neutral,mu_selected,littler,f,
								  sigmaE,optimum,rm,rules);
    }
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
				const internal::region_manager * rm)
    {
      qtrait_model_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      return evolve_qtrait_async_wrapper<sample_n,qtrait_model_rules,decltype(nsam)>(rng,pops,Nvector,Nvector_length,
										     mu_neutral,mu_selected,littler,
										     f,sigmaE,optimum,sample,rm,rules,
										     std::forward<decltype(nsam)>(nsam));
    }

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
				  const internal::region_manager * rm)
    {
      qtrait_model_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      return evolve_qtrait_async_wrapper<pop_properties,qtrait_model_rules>(rng,pops,Nvector,Nvector_length,
									    mu_neutral,mu_selected,littler,
									    f,sigmaE,optimum,sample,rm,rules);
									
    }

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
			       const internal::region_manager * rm)
    {
      qtrait_model_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      return evolve_qtrait_async_wrapper<selected_mut_tracker,qtrait_model_rules>(rng,pops,Nvector,Nvector_length,
										   mu_neutral,mu_selected,littler,
										   f,sigmaE,optimum,track,rm,rules);
    }
    
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
				       const internal::region_manager * rm)
    {
      gbr_model_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      evolve_qtrait_async_wrapper<no_sampling,gbr_model_rules>(rng,pops,Nvector,Nvector_length,
							       mu_neutral,mu_selected,littler,f,
							       sigmaE,optimum,rm,rules);
    }
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
			     const internal::region_manager * rm)
    {
      gbr_model_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      return evolve_qtrait_async_wrapper<sample_n,gbr_model_rules,decltype(nsam)>(rng,pops,Nvector,Nvector_length,
										  mu_neutral,mu_selected,littler,
										  f,sigmaE,optimum,sample,rm,rules,
										  std::forward<decltype(nsam)>(nsam));
    }

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
			       const internal::region_manager * rm)
    {
      gbr_model_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      return evolve_qtrait_async_wrapper<pop_properties,gbr_model_rules>(rng,pops,Nvector,Nvector_length,
									 mu_neutral,mu_selected,littler,
									 f,sigmaE,optimum,sample,rm,rules);
									
    }

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
			    const internal::region_manager * rm)
    {
      gbr_model_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      return evolve_qtrait_async_wrapper<selected_mut_tracker,gbr_model_rules>(rng,pops,Nvector,Nvector_length,
										mu_neutral,mu_selected,littler,
										f,sigmaE,optimum,track,rm,rules);
									
    }
  } //ns qtrait
} //ns fwdpy


