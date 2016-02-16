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
#include <internal/internal.hpp>
#include <qtrait/rules.hpp>
#include <qtrait/details.hpp>
#include <qtrait/evolve_qtrait_sampler.hpp>

#include <gsl/gsl_statistics_double.h>

using namespace std;

using poptype = fwdpy::singlepop_t;

namespace fwdpy
{
  namespace qtrait
  {
    void evolve_qtraits_t( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
			   const unsigned * Nvector,
			   const size_t Nvector_length,
			   const double mu_neutral,
			   const double mu_selected,
			   const double littler,
			   const double f,
			   const double sigmaE,
			   const double optimum,
			   const double VS,
			   const fwdpy::internal::region_manager * rm)
    {
      std::vector<std::thread> threads(pops->size());
      for(unsigned i=0;i<pops->size();++i)
	{
	  threads[i]=std::thread(fwdpy::qtrait::qtrait_sim_details_t<qtrait_model_rules>,
				 gsl_rng_get(rng->get()),
				 pops->operator[](i).get(),
				 Nvector,Nvector_length,
				 mu_neutral,mu_selected,littler,f,sigmaE,optimum,
				 std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
				 std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)),
				 std::move(qtrait_model_rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length))));
	}
      for(unsigned i=0;i<threads.size();++i) threads[i].join();
    }

    // std::vector< std::vector< std::pair<unsigned,qtrait_sample_info_t > > >
    // evolve_qtraits_sample_t( GSLrng_t * rng, std::vector<std::shared_ptr<fwdpy::singlepop_t> > * pops,
    // 			     const unsigned * Nvector,
    // 			     const size_t Nvector_length,
    // 			     const double mu_neutral,
    // 			     const double mu_selected,
    // 			     const double littler,
    // 			     const double f,
    // 			     const double sigmaE,
    // 			     const double optimum,
    // 			     const double VS,
    // 			     const int trackSamples,
    // 			     const unsigned nsam,
    // 			     const fwdpy::internal::region_manager * rm)
    // {
    //   using future_t = std::future<std::vector< std::pair<unsigned,qtrait_sample_info_t > > >;
    //   std::vector<future_t> futures;
    //   for(unsigned i=0;i<pops->size();++i)
    // 	{
    // 	  auto as = std::async(std::launch::async,
    // 			       fwdpy::qtrait::qtrait_sim_details_samples_t<qtrait_model_rules>,
    // 			       gsl_rng_get(rng->get()),
    // 			       pops->operator[](i).get(),
    // 			       Nvector,Nvector_length,
    // 			       mu_neutral,mu_selected,littler,f,sigmaE,optimum,trackSamples,nsam,
    // 			       std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
    // 			       std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)),
    // 			       std::move(qtrait_model_rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length)))
    // 			       );
    // 	  futures.emplace_back( std::move(as) );
    // 	}
    //   std::vector< std::vector< std::pair<unsigned,qtrait_sample_info_t > > > rv(pops->size());
    //   for(unsigned i=0;i<pops->size();++i)
    // 	{
    // 	  rv[i]=std::move(futures[i].get());
    // 	}
    //   return rv;
    // }

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
  } //ns qtrait
} //ns fwdpy


