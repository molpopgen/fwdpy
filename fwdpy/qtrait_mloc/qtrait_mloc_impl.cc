
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
#include <gsl/gsl_statistics_double.h>

#include "types.hpp"
#include "qtrait_evolve_mlocus.hpp"
#include "qtrait_mloc_rules.hpp"

using namespace std;

using poptype = fwdpy::singlepop_t;

namespace fwdpy
{
  namespace qtrait
  {
    std::vector<KTfwd::extensions::shmodel> make_GaussianDFE(const std::vector<double> & sigma_mus)
    //! Convenience function for what is below: makes additive Gaussian DFE for all elements in sigma_mus
    {
      std::vector<KTfwd::extensions::shmodel> rv;
      for(  auto & s : sigma_mus )
	{
	  rv.emplace_back(KTfwd::extensions::gaussian(s),KTfwd::extensions::constant(1.0));
	}
      return rv;
    }
    
    std::vector<KTfwd::extensions::shmodel> make_ExponentialDFE(const std::vector<double> & sigma_mus)
    //! Convenience function for what is below: makes additive Gaussian DFE for all elements in sigma_mus
    {
      std::vector<KTfwd::extensions::shmodel> rv;
      for(  auto & s : sigma_mus )
	{
	  rv.emplace_back(KTfwd::extensions::exponential(s),KTfwd::extensions::constant(1.0));
	}
      return rv;
    }

    std::vector<sample_n<multilocus_t>::final_t>
    evolve_qtrait_mloc_sample_async( GSLrng_t * rng,
				     GSLrng_t * rng_sampling,
				     std::vector<std::shared_ptr<multilocus_t> > * pops,
				     const unsigned * Nvector,
				     const size_t Nvector_length,
				     const std::vector<double> & neutral_mutation_rates,
				     const std::vector<double> & selected_mutation_rates,
				     const std::vector<double> & sigma_mus,
				     const std::vector<double> & within_region_rec_rates,
				     const std::vector<double> & between_region_rec_rates,
				     const double f,
				     const double sigmaE,
				     const double optimum,
				     const double VS,
				     const int sample,
				     const unsigned nsam)
    {
      qtrait_mloc_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      return evolve_qtrait_mloc_async_wrapper<fwdpy::sample_n<multilocus_t>,
					      qtrait_mloc_rules,
					      decltype(nsam),
					      const gsl_rng *>(rng,
							       pops,
							       Nvector,
							       Nvector_length,
							       neutral_mutation_rates,
							       selected_mutation_rates,
							       make_GaussianDFE(sigma_mus),
							       within_region_rec_rates,
							       between_region_rec_rates,
							       f,
							       sample,
							       rules,
							       std::forward<decltype(nsam)>(nsam),rng_sampling->get());
    }

    //Sample nothing
    void evolve_qtrait_mloc_no_sampling_async( GSLrng_t * rng,
					       std::vector<std::shared_ptr<multilocus_t> > * pops,
					       const unsigned * Nvector,
					       const size_t Nvector_length,
					       const std::vector<double> & neutral_mutation_rates,
					       const std::vector<double> & selected_mutation_rates,
					       const std::vector<double> & sigma_mus,
					       const std::vector<double> & within_region_rec_rates,
					       const std::vector<double> & between_region_rec_rates,
					       const double f,
					       const double sigmaE,
					       const double optimum,
					       const double VS)
    {
      qtrait_mloc_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      evolve_qtrait_mloc_async_wrapper<fwdpy::no_sampling,
				       qtrait_mloc_rules>(rng,
							  pops,
							  Nvector,
							  Nvector_length,
							  neutral_mutation_rates,
							  selected_mutation_rates,
							  make_GaussianDFE(sigma_mus),
							  within_region_rec_rates,
							  between_region_rec_rates,
							  f,
							  0,
							  rules);
    }
    
    //Sample quant. genetics params from pop
    std::vector<pop_properties::final_t>
    evolve_qtrait_mloc_popstats_async( GSLrng_t * rng,
				       std::vector<std::shared_ptr<multilocus_t> > * pops,
				       const unsigned * Nvector,
				       const size_t Nvector_length,
				       const std::vector<double> & neutral_mutation_rates,
				       const std::vector<double> & selected_mutation_rates,
				       const std::vector<double> & sigma_mus,
				       const std::vector<double> & within_region_rec_rates,
				       const std::vector<double> & between_region_rec_rates,
				       const double f,
				       const double sigmaE,
				       const double optimum,
				       const double VS,
				       const int sample)
    {
      qtrait_mloc_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      return  evolve_qtrait_mloc_async_wrapper<fwdpy::pop_properties,
					       qtrait_mloc_rules,
					       decltype(optimum)>(rng,
								  pops,
								  Nvector,
								  Nvector_length,
								  neutral_mutation_rates,
								  selected_mutation_rates,
								  make_GaussianDFE(sigma_mus),
								  within_region_rec_rates,
								  between_region_rec_rates,
								  f,
								  sample,
								  rules,
								  std::forward<decltype(optimum)>(optimum));

    }

    //Causative mutation frequency trajectories
    std::vector<selected_mut_tracker::final_t>
    evolve_qtrait_mloc_track_async( GSLrng_t * rng,
				    std::vector<std::shared_ptr<multilocus_t> > * pops,
				    const unsigned * Nvector,
				    const size_t Nvector_length,
				    const std::vector<double> & neutral_mutation_rates,
				    const std::vector<double> & selected_mutation_rates,
				    const std::vector<double> & sigma_mus,
				    const std::vector<double> & within_region_rec_rates,
				    const std::vector<double> & between_region_rec_rates,
				    const double f,
				    const double sigmaE,
				    const double optimum,
				    const double VS,
				    const int sample)
    {
      qtrait_mloc_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      return  evolve_qtrait_mloc_async_wrapper<fwdpy::selected_mut_tracker,
					       qtrait_mloc_rules>(rng,
								  pops,
								  Nvector,
								  Nvector_length,
								  neutral_mutation_rates,
								  selected_mutation_rates,
								  make_GaussianDFE(sigma_mus),
								  within_region_rec_rates,
								  between_region_rec_rates,
								  f,
								  sample,
								  rules);
    }

    std::vector<additive_variance::final_t>
    evolve_qtrait_mloc_VA_async( GSLrng_t * rng,
				 std::vector<std::shared_ptr<multilocus_t> > * pops,
				 const unsigned * Nvector,
				 const size_t Nvector_length,
				 const std::vector<double> & neutral_mutation_rates,
				 const std::vector<double> & selected_mutation_rates,
				 const std::vector<double> & sigma_mus,
				 const std::vector<double> & within_region_rec_rates,
				 const std::vector<double> & between_region_rec_rates,
				 const double f,
				 const double sigmaE,
				 const double optimum,
				 const double VS,
				 const int sample)
    {
      qtrait_mloc_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      return  evolve_qtrait_mloc_async_wrapper<fwdpy::additive_variance,
					       qtrait_mloc_rules>(rng,
								  pops,
								  Nvector,
								  Nvector_length,
								  neutral_mutation_rates,
								  selected_mutation_rates,
								  make_GaussianDFE(sigma_mus),
								  within_region_rec_rates,
								  between_region_rec_rates,
								  f,
								  sample,
								  rules);
    }
    //POWER MEAN FUNCTIONS
        std::vector<sample_n<multilocus_t>::final_t>
    evolve_qtrait_mloc_pm_sample_async( GSLrng_t * rng,
				     GSLrng_t * rng_sampling,
				     std::vector<std::shared_ptr<multilocus_t> > * pops,
				     const unsigned * Nvector,
				     const size_t Nvector_length,
				     const std::vector<double> & neutral_mutation_rates,
				     const std::vector<double> & selected_mutation_rates,
				     const std::vector<double> & sigma_mus,
				     const std::vector<double> & within_region_rec_rates,
				     const std::vector<double> & between_region_rec_rates,
				     const double f,
				     const double sigmaE,
				     const double optimum,
				     const double VS,
				     const std::vector<double> & SLd,
					 const double & SLp,
			  		 const std::vector<double> & MLd,
			  		 const double & MLp,
				     const int sample,
				     const unsigned nsam)
    {
      qtrait_mloc_pm_rules rules(sigmaE,optimum,VS,SLd,SLp,MLd,MLp,*std::max_element(Nvector,Nvector+Nvector_length));
      return evolve_qtrait_mloc_async_wrapper<fwdpy::sample_n<multilocus_t>,
					      qtrait_mloc_pm_rules,
					      decltype(nsam),
					      const gsl_rng *>(rng,
							       pops,
							       Nvector,
							       Nvector_length,
							       neutral_mutation_rates,
							       selected_mutation_rates,
							       make_ExponentialDFE(sigma_mus),
							       within_region_rec_rates,
							       between_region_rec_rates,
							       f,
							       sample,
							       rules,
							       std::forward<decltype(nsam)>(nsam),rng_sampling->get());
    }

    //Sample nothing
    void evolve_qtrait_mloc_pm_no_sampling_async( GSLrng_t * rng,
					       std::vector<std::shared_ptr<multilocus_t> > * pops,
					       const unsigned * Nvector,
					       const size_t Nvector_length,
					       const std::vector<double> & neutral_mutation_rates,
					       const std::vector<double> & selected_mutation_rates,
					       const std::vector<double> & sigma_mus,
					       const std::vector<double> & within_region_rec_rates,
					       const std::vector<double> & between_region_rec_rates,
					       const double f,
					       const double sigmaE,
					       const double optimum,
					       const double VS,
					       const std::vector<double> & SLd,
					 	   const double & SLp,
			  		 	   const std::vector<double> & MLd,
			  		 	   const double & MLp)
    {
      qtrait_mloc_pm_rules rules(sigmaE,optimum,VS,SLd,SLp,MLd,MLp,*std::max_element(Nvector,Nvector+Nvector_length));
      evolve_qtrait_mloc_async_wrapper<fwdpy::no_sampling,
				       qtrait_mloc_pm_rules>(rng,
							  pops,
							  Nvector,
							  Nvector_length,
							  neutral_mutation_rates,
							  selected_mutation_rates,
							  make_ExponentialDFE(sigma_mus),
							  within_region_rec_rates,
							  between_region_rec_rates,
							  f,
							  0,
							  rules);
    }
    
    //Sample quant. genetics params from pop
    std::vector<pop_properties::final_t>
    evolve_qtrait_mloc_pm_popstats_async( GSLrng_t * rng,
				       std::vector<std::shared_ptr<multilocus_t> > * pops,
				       const unsigned * Nvector,
				       const size_t Nvector_length,
				       const std::vector<double> & neutral_mutation_rates,
				       const std::vector<double> & selected_mutation_rates,
				       const std::vector<double> & sigma_mus,
				       const std::vector<double> & within_region_rec_rates,
				       const std::vector<double> & between_region_rec_rates,
				       const double f,
				       const double sigmaE,
				       const double optimum,
				       const double VS,
				       const std::vector<double> & SLd,
					   const double & SLp,
			  		   const std::vector<double> & MLd,
			  		   const double & MLp,
				       const int sample)
    {
      qtrait_mloc_pm_rules rules(sigmaE,optimum,VS,SLd,SLp,MLd,MLp,*std::max_element(Nvector,Nvector+Nvector_length));
      return  evolve_qtrait_mloc_async_wrapper<fwdpy::pop_properties,
					       qtrait_mloc_pm_rules,
					       decltype(optimum)>(rng,
								  pops,
								  Nvector,
								  Nvector_length,
								  neutral_mutation_rates,
								  selected_mutation_rates,
								  make_ExponentialDFE(sigma_mus),
								  within_region_rec_rates,
								  between_region_rec_rates,
								  f,
								  sample,
								  rules,
								  std::forward<decltype(optimum)>(optimum));

    }

    //Causative mutation frequency trajectories
    std::vector<selected_mut_tracker::final_t>
    evolve_qtrait_mloc_pm_track_async( GSLrng_t * rng,
				    std::vector<std::shared_ptr<multilocus_t> > * pops,
				    const unsigned * Nvector,
				    const size_t Nvector_length,
				    const std::vector<double> & neutral_mutation_rates,
				    const std::vector<double> & selected_mutation_rates,
				    const std::vector<double> & sigma_mus,
				    const std::vector<double> & within_region_rec_rates,
				    const std::vector<double> & between_region_rec_rates,
				    const double f,
				    const double sigmaE,
				    const double optimum,
				    const double VS,
				    const std::vector<double> & SLd,
					const double & SLp,
			  		const std::vector<double> & MLd,
			  		const double & MLp,
				    const int sample)
    {
      qtrait_mloc_pm_rules rules(sigmaE,optimum,VS,SLd,SLp,MLd,MLp,*std::max_element(Nvector,Nvector+Nvector_length));
      return  evolve_qtrait_mloc_async_wrapper<fwdpy::selected_mut_tracker,
					       qtrait_mloc_pm_rules>(rng,
								  pops,
								  Nvector,
								  Nvector_length,
								  neutral_mutation_rates,
								  selected_mutation_rates,
								  make_ExponentialDFE(sigma_mus),
								  within_region_rec_rates,
								  between_region_rec_rates,
								  f,
								  sample,
								  rules);
    }

    std::vector<additive_variance::final_t>
    evolve_qtrait_mloc_pm_VA_async( GSLrng_t * rng,
				 std::vector<std::shared_ptr<multilocus_t> > * pops,
				 const unsigned * Nvector,
				 const size_t Nvector_length,
				 const std::vector<double> & neutral_mutation_rates,
				 const std::vector<double> & selected_mutation_rates,
				 const std::vector<double> & sigma_mus,
				 const std::vector<double> & within_region_rec_rates,
				 const std::vector<double> & between_region_rec_rates,
				 const double f,
				 const double sigmaE,
				 const double optimum,
				 const double VS,
				 const std::vector<double> & SLd,
				 const double & SLp,
			  	 const std::vector<double> & MLd,
			  	 const double & MLp,
				 const int sample)
    {
      qtrait_mloc_pm_rules rules(sigmaE,optimum,VS,SLd,SLp,MLd,MLp,*std::max_element(Nvector,Nvector+Nvector_length));
      return  evolve_qtrait_mloc_async_wrapper<fwdpy::additive_variance,
					       qtrait_mloc_pm_rules>(rng,
								  pops,
								  Nvector,
								  Nvector_length,
								  neutral_mutation_rates,
								  selected_mutation_rates,
								  make_ExponentialDFE(sigma_mus),
								  within_region_rec_rates,
								  between_region_rec_rates,
								  f,
								  sample,
								  rules);
    }

  }
}
