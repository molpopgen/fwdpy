#include <future>
#include <functional>
#include <evolve_regions_sampler.hpp>

using namespace std;

namespace fwdpy
{
  std::vector< std::vector<std::pair<unsigned,KTfwd::sep_sample_t> > >
  evolve_regions_sample_async( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
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
      using rv_t = std::vector<std::pair<unsigned,KTfwd::sep_sample_t> >;
      using future_t = std::future<rv_t>;
      vector<future_t> futures;
      auto sampler_t = sample_n::make_bound_t(nsam);
      for( unsigned i=0;i<pops->size() ; ++i )
	{
	  futures.emplace_back( async(launch::async,
	  			      evolve_regions_sampler_details<sample_n::bound_t>,
	  			      pops->operator[](i).get(),gsl_rng_get(rng->get()),Nvector,Nvector_len,
	  			      mu_neutral,mu_selected,littler,f,fitness,sampler_t,sample,
	  			      std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
	  			      std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)))
	  			);
	}
      vector<rv_t> rv(futures.size());
      for(unsigned i=0;i<futures.size();++i) rv[i]=futures[i].get();
      return rv;
    }
}
