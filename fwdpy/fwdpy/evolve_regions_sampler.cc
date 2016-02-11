#include <future>
#include <iterator>
#include <functional>
#include <evolve_regions_sampler.hpp>

using namespace std;

namespace fwdpy
{
  std::vector<sample_n::final_t>
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
      using future_t = std::future<sample_n::final_t>;
      vector<future_t> futures;
      for( unsigned i=0;i<pops->size() ; ++i )
	{
	  futures.emplace_back( async(launch::async,
	  			      evolve_regions_sampler_details<sample_n,decltype(nsam)>,
	  			      pops->operator[](i).get(),gsl_rng_get(rng->get()),Nvector,Nvector_len,
	  			      mu_neutral,mu_selected,littler,f,fitness,sample,
	  			      std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
	  			      std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)),
				      std::forward<decltype(nsam)>(nsam))
	  			);
	}
      vector<sample_n::final_t> rv(futures.size());
      for(unsigned i=0;i<futures.size();++i) rv[i]=futures[i].get();
      return rv;
    }

    std::vector<get_selected_mut_data::final_t>
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
      using future_t = std::future<get_selected_mut_data::final_t>;
      vector<future_t> futures;
      for( unsigned i=0;i<pops->size() ; ++i )
    	{
    	  futures.emplace_back( async(launch::async,
    	  			      evolve_regions_sampler_details<get_selected_mut_data>,
    	  			      pops->operator[](i).get(),gsl_rng_get(rng->get()),Nvector,Nvector_len,
    	  			      mu_neutral,mu_selected,littler,f,fitness,sample,
    	  			      std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
    	  			      std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw))
    				      )
    	  			);
    	}
      vector<get_selected_mut_data::final_t> rv;
      for(std::size_t i=0;i<futures.size();++i ) rv.emplace_back(futures[i].get());
      return rv;
    }
}
