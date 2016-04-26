
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
#include "sampler_sample_n.hpp"
#include "qtrait_mloc_rules.hpp"

using namespace std;

using poptype = fwdpy::singlepop_t;

namespace fwdpy
{
  namespace qtrait
  {
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
							       sigma_mus,
							       within_region_rec_rates,
							       between_region_rec_rates,
							       f,
							       sample,
							       rules,
							       std::forward<decltype(nsam)>(nsam),rng_sampling->get());
    }
  }
}
