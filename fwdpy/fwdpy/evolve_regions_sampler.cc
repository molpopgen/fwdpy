#include <future>
#include <iterator>
#include <functional>
#include "evolve_regions_sampler.hpp"
#include "sampler_no_sampling.hpp"

using namespace std;

namespace fwdpy
{
  void evolve_regions_no_sampling_async(GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
					const unsigned * Nvector,
					const size_t Nvector_length,
					const double mu_neutral,
					const double mu_selected,
					const double littler,
					const double f,
					const internal::region_manager * rm,
					const char * fitness)
  {
    evolve_regions_async_wrapper<no_sampling>(rng,pops,Nvector,Nvector_length,mu_neutral,mu_selected,littler,f,
					      0,
					      rm,fitness);
  }
  
  std::vector<sample_n<singlepop_t>::final_t>
  evolve_regions_sample_async( GSLrng_t * rng, GSLrng_t * rng_sampler,
			       std::vector<std::shared_ptr<singlepop_t> > * pops,
			       const unsigned * Nvector,
			       const size_t Nvector_len,
			       const double mu_neutral,
			       const double mu_selected,
			       const double littler,
			       const double f,
			       const int sample,
			       const unsigned nsam,
			       const internal::region_manager * rm,
			       const char * fitness)
  {
    return evolve_regions_async_wrapper<sample_n<singlepop_t>,
					decltype(nsam),
					const gsl_rng *>(rng,pops,Nvector,Nvector_len,
							 mu_neutral,mu_selected,littler,f,sample,
							 rm,fitness,
							 std::forward<decltype(nsam)>(nsam),
							 std::forward<const gsl_rng *>(rng_sampler->get()));
  }

  std::vector<selected_mut_tracker::final_t>
  evolve_regions_track_async( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
			      const unsigned * Nvector,
			      const size_t Nvector_len,
			      const double mu_neutral,
			      const double mu_selected,
			      const double littler,
			      const double f,
			      const int sample,
			      const internal::region_manager * rm,
			      const char * fitness)
  {
    return evolve_regions_async_wrapper<selected_mut_tracker>(rng,pops,Nvector,Nvector_len,mu_neutral,mu_selected,littler,f,sample,
							      rm,fitness);
  }
}
