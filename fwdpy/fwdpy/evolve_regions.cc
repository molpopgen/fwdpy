#include <thread>
#include <future>
#include <functional>
#include <fwdpp/diploid.hh>
#include <fwdpp/extensions/regions.hpp>
#include "types.hpp"
#include "reserve.hpp"
#include "metapop.hpp"
#include "evolve_regions.hpp"

namespace fwdpy {
  void split_and_evolve_details(metapop_t * mpop,
				unsigned long seed,
				const unsigned * Nvector_A,
				const size_t Nvector_A_len,
				const unsigned * Nvector_B,
				const size_t Nvector_B_len,
				const double neutral,
				const double selected,
				const double recrate,
				const std::vector<double> fs,
				KTfwd::extensions::discrete_mut_model && __m,
				KTfwd::extensions::discrete_rec_model && __recmap,
				const char * fitness)
  {
    const size_t simlen = Nvector_A_len;
    const double mu_tot = neutral + selected;

    auto x = std::max_element(Nvector_A,Nvector_A+Nvector_A_len);
    assert(x!=Nvector_A+Nvector_A_len);
    auto y = std::max_element(Nvector_B,Nvector_B+Nvector_B_len);
    assert(y!=Nvector_B+Nvector_B_len);
    //We reserve space for 2x b/c there are two demes
    reserve_space(mpop->gametes,mpop->mutations,2*std::max(*x,*y),2.*mu_tot);

    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
    KTfwd::extensions::discrete_mut_model m(std::move(__m));
    KTfwd::extensions::discrete_rec_model recmap(std::move(__recmap));
    const auto recpos = KTfwd::extensions::bind_drm(recmap,mpop->gametes,mpop->mutations,
						    rng,recrate);

    //The fitness model
    fwdpy::metapop_t::fitness_t dipfit = std::bind(KTfwd::multiplicative_diploid(),
						   std::placeholders::_1,
						   std::placeholders::_2,
						   std::placeholders::_3,
						   2.);

    if( std::string(fitness) == "additive" )
      {
     	dipfit = std::bind(KTfwd::additive_diploid(),
			   std::placeholders::_1,
			   std::placeholders::_2,
			   std::placeholders::_3,
			   2.);
      }

    //vector of fitness models
    std::vector<decltype(dipfit)> fmodels({dipfit,dipfit});

    //FIX later
    const double migrate = 0.;

    //Finally, we can evolve this thing
    for( size_t g = 0 ; g < simlen ; ++g, ++mpop->generation )
      {
	std::vector<unsigned> Ns_next({unsigned(Nvector_A[g]),unsigned(Nvector_A[g])});
	unsigned N = std::accumulate(Ns_next.begin(),Ns_next.end(),0u);
	std::vector<double> wbars = sample_diploid(rng,
						   mpop->gametes,
						   mpop->diploids,
						   mpop->mutations,
						   mpop->mcounts,
						   &mpop->Ns[0],
						   &Ns_next[0],
						   mu_tot,
						   KTfwd::extensions::bind_dmm(m,mpop->mutations,mpop->mut_lookup,rng,neutral,selected,mpop->generation),
						   recpos,
						   fmodels,
						   std::bind(migpop,std::placeholders::_1,rng,migrate),
						   mpop->neutral,
						   mpop->selected,
						   &fs[0]);

	mpop->Ns = Ns_next;
	//4*N b/c it needs to be fixed in the metapopulation
	KTfwd::update_mutations(mpop->mutations,mpop->fixations,mpop->fixation_times,mpop->mut_lookup,mpop->mcounts,mpop->generation,4*N);
      }
    gsl_rng_free(rng);
  }

  void split_and_evolve_t(GSLrng_t * rng,
			  std::vector<std::shared_ptr<metapop_t> > * mpops,
			  const unsigned * Nvector_A,
			  const size_t Nvector_A_len,
			  const unsigned * Nvector_B,
			  const size_t Nvector_B_len,
			  const double neutral,
			  const double selected,
			  const double recrate,
			  const std::vector<double> & fs,
			  const fwdpy::internal::region_manager * rm,
			  const char * fitness)
  {
    std::vector<std::thread> threads(mpops->size());
    for(unsigned i=0;i<mpops->size();++i)
      {
  	threads[i]=std::thread(split_and_evolve_details,
			       mpops->operator[](i).get(),
			       gsl_rng_get(rng->get()),
			       Nvector_A,
			       Nvector_A_len,
			       Nvector_B,
			       Nvector_B_len,
			       neutral,
			       selected,
			       recrate,
			       fs,
			       std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
			       std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)),
			       fitness);
      }
    for(unsigned i=0;i<threads.size();++i) threads[i].join();
  }
}

