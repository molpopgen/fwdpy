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
#include <types.hpp>
#include <internal/internal.hpp>
#include <qtrait/rules.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/experimental/sample_diploid.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/popgenmut.hpp>

#include <qtrait/details.hpp>


#include <gsl/gsl_statistics_double.h>

using namespace std;

using poptype = fwdpy::singlepop_t;

namespace fwdpy
{
  //Get properties out from the population
  std::vector<qtrait_stats_cython> convert_qtrait_stats( const fwdpy::singlepop_t * pop )
  {
    using namespace fwdpy;
    std::vector<qtrait_stats_cython> rv;
    for( const auto & i : pop->qstats )
      {
	rv.emplace_back( "VG",i[std::size_t(qtrait_stat_names::VG)],i[std::size_t(qtrait_stat_names::GEN)] );
	rv.emplace_back( "VE",i[std::size_t(qtrait_stat_names::VE)],i[std::size_t(qtrait_stat_names::GEN)] );
	rv.emplace_back( "leading_q",i[std::size_t(qtrait_stat_names::PLF)],i[std::size_t(qtrait_stat_names::GEN)] );
	rv.emplace_back( "leading_e",i[std::size_t(qtrait_stat_names::LE)],i[std::size_t(qtrait_stat_names::GEN)] );
	rv.emplace_back( "max_expl",i[std::size_t(qtrait_stat_names::MAXEXP)],i[std::size_t(qtrait_stat_names::GEN)] );
	rv.emplace_back( "ebar",i[std::size_t(qtrait_stat_names::EBAR)],i[std::size_t(qtrait_stat_names::GEN)] );
	rv.emplace_back( "wbar",i[std::size_t(qtrait_stat_names::WBAR)],i[std::size_t(qtrait_stat_names::GEN)] );
      }
    return rv;
  };

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
			   const int track,
			   const int trackStats,
			   const fwdpy::internal::region_manager * rm)
    {
      std::vector<std::thread> threads(pops->size());
      for(unsigned i=0;i<pops->size();++i)
	{
	  threads[i]=std::thread(fwdpy::qtrait::qtrait_sim_details_t<qtrait_model_rules>,
				 gsl_rng_get(rng->get()),
				 pops->operator[](i).get(),
				 Nvector,Nvector_length,
				 mu_neutral,mu_selected,littler,f,sigmaE,optimum,track,trackStats,
				 std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
				 std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)),
				 std::move(qtrait_model_rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length))));
	}
      for(unsigned i=0;i<threads.size();++i) threads[i].join();
    }

    std::vector< std::vector< std::pair<unsigned,qtrait_sample_info_t > > >
    evolve_qtraits_sample_t( GSLrng_t * rng, std::vector<std::shared_ptr<fwdpy::singlepop_t> > * pops,
			     const unsigned * Nvector,
			     const size_t Nvector_length,
			     const double mu_neutral,
			     const double mu_selected,
			     const double littler,
			     const double f,
			     const double sigmaE,
			     const double optimum,
			     const double VS,
			     const int trackSamples,
			     const unsigned nsam,
			     const fwdpy::internal::region_manager * rm)
    {
      using future_t = std::future<std::vector< std::pair<unsigned,qtrait_sample_info_t > > >;
      std::vector<future_t> futures;
      for(unsigned i=0;i<pops->size();++i)
	{
	  auto as = std::async(std::launch::async,
			       fwdpy::qtrait::qtrait_sim_details_samples_t<qtrait_model_rules>,
			       gsl_rng_get(rng->get()),
			       pops->operator[](i).get(),
			       Nvector,Nvector_length,
			       mu_neutral,mu_selected,littler,f,sigmaE,optimum,trackSamples,nsam,
			       std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
			       std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)),
			       std::move(qtrait_model_rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length)))
			       );
	  futures.emplace_back( std::move(as) );
	}
      std::vector< std::vector< std::pair<unsigned,qtrait_sample_info_t > > > rv(pops->size());
      for(unsigned i=0;i<pops->size();++i)
	{
	  rv[i]=std::move(futures[i].get());
	}
      return rv;
    }
  } //ns qtrait
} //ns fwdpy


