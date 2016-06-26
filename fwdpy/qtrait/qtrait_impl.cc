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
#include <gsl/gsl_statistics_double.h>

#include "types.hpp"
#include "internal_region_manager.hpp"
#include "qtrait_evolve_rules.hpp"
#include "qtrait_details.hpp"
#include "qtrait_evolve.hpp"
#include "sampler_no_sampling.hpp"
#include "sampler_additive_variance.hpp"

using namespace std;

using poptype = fwdpy::singlepop_t;

namespace fwdpy
{
  namespace qtrait
  {
    void evolve_regions_qtrait_cpp( GSLrng_t * rng,
				    std::vector<std::shared_ptr<singlepop_t> > * pops,
				    std::vector<std::unique_ptr<sampler_base> > & samplers,
				    const unsigned * Nvector,
				    const size_t Nvector_length,
				    const double neutral,
				    const double selected,
				    const double recrate,
				    const double f,
				    const double sigmaE,
				    const double optimum,
				    const double VS,
				    const int interval,
				    const internal::region_manager * rm,
				    const singlepop_fitness & fitness)
    {
      if(neutral < 0. || selected < 0. || recrate < 0.)
	{
	  throw std::runtime_error("mutation and recombination rates must all be non-negative.");
	}
      if(samplers.size()!=pops->size())
	{
	  throw std::runtime_error("length of samplers != length of population container");
	}
      if(f<0.||f>1.) throw std::runtime_error("selfing probabilty must be 0<=f<=1.");
      if(interval<0) throw std::runtime_error("sampling interval must be non-negative");
      std::vector<std::thread> threads;
      qtrait_model_rules rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length));
      for(std::size_t i=0;i<pops->size();++i)
	{
	  threads.emplace_back( std::thread(evolve_regions_qtrait_sampler_cpp_details<qtrait_model_rules>,
					    pops->operator[](i).get(),
					    gsl_rng_get(rng->get()),
					    Nvector,
					    Nvector_length,
					    neutral,
					    selected,
					    recrate,
					    f,
					    sigmaE,
					    optimum,
					    VS,
					    fitness,
					    interval,
					    KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks),
					    KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw),
					    std::ref(*samplers[i]),
					    rules
					    )
				);	
	}
      for(auto & t : threads) t.join();
    }
  } //ns qtrait
} //ns fwdpy


