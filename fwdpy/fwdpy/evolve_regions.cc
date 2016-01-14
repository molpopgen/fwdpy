#include <thread>
#include <future>
#include <functional>
#include <fwdpp/diploid.hh>
#include <fwdpp/extensions/regions.hpp>
#include <types.hpp>
#include <metapop.hpp>
#include <evolve_regions.hpp>

namespace fwdpy {

  void evolve_regions_details( fwdpy::singlepop_t * pop,
			       unsigned long seed,
			       const unsigned * Nvector,
			       const size_t Nvector_len,
			       const double neutral,
			       const double selected,
			       const double recrate,
			       const double f,
			       const int track,
			       const char * fitness,
			       KTfwd::extensions::discrete_mut_model && __m,
			       KTfwd::extensions::discrete_rec_model && __recmap)
  {    
    const size_t simlen = Nvector_len;
    
    const double mu_tot = neutral + selected;
    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
    KTfwd::extensions::discrete_mut_model m(std::move(__m));
    KTfwd::extensions::discrete_rec_model recmap(std::move(__recmap));
    fwdpy::singlepop_t pop2(std::move(*pop));
    //Recombination policy: more complex than the standard case...
    const auto recpos = KTfwd::extensions::bind_drm(recmap,pop2.gametes,pop2.mutations,
						    rng,recrate);

    /*
      The fitness model.

      Normally, we'd declare dipfit as "auto", but that won't work here b/c there is the chance 
      that we have to re-assign it using an additive model based on input from calling environment.

      The std::bind signature has a different type for the two models, and thus we must coerce it to
      the proper function signature, which is a member typedef provided by the fwdpp sugar type 
      from which fwdpy::singlepop_t inherits
    */
    fwdpy::singlepop_t::fitness_t dipfit = std::bind(KTfwd::multiplicative_diploid(),
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

    for( size_t g = 0 ; g < simlen ; ++g, ++pop2.generation )
      {
	const unsigned nextN = 	*(Nvector+g);
	KTfwd::sample_diploid(rng,
			      pop2.gametes,  
			      pop2.diploids, 
			      pop2.mutations,
			      pop2.mcounts,
			      pop2.N,
			      nextN,
			      mu_tot,
			      KTfwd::extensions::bind_dmm(m,pop2.mutations,pop2.mut_lookup,rng,neutral,selected,pop2.generation),
			      recpos,
			      dipfit,
			      pop->neutral,pop->selected,
			      f);
	if (track) pop2.updateTraj();
	pop2.N=nextN;
	KTfwd::update_mutations(pop2.mutations,pop2.fixations,pop2.fixation_times,pop2.mut_lookup,pop2.mcounts,pop2.generation,2*nextN);
	assert(KTfwd::check_sum(pop2.gametes,2*nextN));
      }
    //Update population's size variable to be the current pop size
    pop2.N = unsigned(pop2.diploids.size());
    //cleanup
    gsl_rng_free(rng);
    //restore data
    *pop=std::move(pop2);
  }

  void evolve_regions_t( GSLrng_t * rng, std::shared_ptr<singlepop_t> pop,
			 const unsigned * Nvector,
			 const size_t Nvector_len,
			 const double mu_neutral,
			 const double mu_selected,
			 const double littler,
			 const double f,
			 const int track,
			 const fwdpy::internal::region_manager * rm,
			 const char * fitness)
  {
    evolve_regions_details(pop.get(),gsl_rng_get(rng->get()),Nvector,Nvector_len,
			   mu_neutral,mu_selected,littler,f,track,fitness,
			   std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
			   std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)));
  }
      
  void evolve_regions_t( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
			 const unsigned * Nvector,
			 const size_t Nvector_len,
			 const double mu_neutral,
			 const double mu_selected,
			 const double littler,
			 const double f,
			 const int track,
			 const fwdpy::internal::region_manager * rm,
			 const char * fitness)
  {
    std::vector<std::thread> threads(pops->size());
    for(unsigned i=0;i<pops->size();++i)
      {
	threads[i]=std::thread(evolve_regions_details,pops->operator[](i).get(),gsl_rng_get(rng->get()),Nvector,Nvector_len,
			       mu_neutral,mu_selected,littler,f,track,fitness,
			       std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
			       std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)));    
      }
    for(unsigned i=0;i<threads.size();++i) threads[i].join();
  }

  std::shared_ptr<fwdpy::singlepop_t> evolve_regions_details_async(unsigned long seed,
								   const unsigned * Nvector,
								   const size_t Nvector_len,
								   const double neutral,
								   const double selected,
								   const double recrate,
								   const double f,
								   const int track,
								   const char * fitness,
								   const fwdpy::internal::region_manager * rm)
  {    
    const size_t simlen = Nvector_len;
    
    const double mu_tot = neutral + selected;
    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
    KTfwd::extensions::discrete_mut_model m(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks);//std::move(__m));
    KTfwd::extensions::discrete_rec_model recmap(rm->rb,rm->rw,rm->rw);//std::move(__recmap));
    fwdpy::singlepop_t pop2(Nvector[0]);
    //Recombination policy: more complex than the standard case...
    const auto recpos = KTfwd::extensions::bind_drm(recmap,pop2.gametes,pop2.mutations,
						    rng,recrate);

    //The fitness model
    fwdpy::singlepop_t::fitness_t dipfit = std::bind(KTfwd::multiplicative_diploid(),
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
    //This may not be the best thing, long-term, design-wise...

    for( size_t g = 0 ; g < simlen ; ++g, ++pop2.generation )
      {
	const unsigned nextN = 	*(Nvector+g);
	KTfwd::sample_diploid(rng,
			      pop2.gametes,  
			      pop2.diploids, 
			      pop2.mutations,
			      pop2.mcounts,
			      pop2.N,
			      nextN,
			      mu_tot,
			      KTfwd::extensions::bind_dmm(m,pop2.mutations,pop2.mut_lookup,rng,neutral,selected,pop2.generation),
			      recpos,
			      dipfit,
			      pop2.neutral,
			      pop2.selected,
			      f);
	if (track) pop2.updateTraj();
	pop2.N=nextN;
	KTfwd::update_mutations(pop2.mutations,pop2.fixations,pop2.fixation_times,pop2.mut_lookup,pop2.mcounts,pop2.generation,2*nextN);
	assert(KTfwd::check_sum(pop2.gametes,2*nextN));
      }
    //Update population's size variable to be the current pop size
    pop2.N = unsigned(pop2.diploids.size());
    //cleanup
    gsl_rng_free(rng);
    return std::make_shared<fwdpy::singlepop_t>(std::move(pop2));
  }
  
  std::vector<std::shared_ptr<singlepop_t> >  evolve_regions_async(const unsigned npops,
								   GSLrng_t * rng, 
								   const unsigned * Nvector,
								   const size_t Nvector_len,
								   const double mu_neutral,
								   const double mu_selected,
								   const double littler,
								   const double f,
								   const int track,
								   const fwdpy::internal::region_manager * rm,
								   const char * fitness)
  {
    std::vector<std::future<std::shared_ptr<singlepop_t> > > futures;
    for(unsigned i=0;i<npops;++i)
      {
	futures.emplace_back(std::async(std::launch::async,evolve_regions_details_async,
					gsl_rng_get(rng->get()),Nvector,Nvector_len,
					mu_neutral,mu_selected,littler,f,track,fitness,rm));
	//std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
	//				 std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw))));
      }
    std::vector<std::shared_ptr<singlepop_t> > rv;
    for_each(std::begin(futures),std::end(futures),[&rv](std::future<std::shared_ptr<singlepop_t> > & fut) {
	// 	 //fut.wait();
	rv.emplace_back(std::move(fut.get()));
      });
    return rv;
  }

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

    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
    KTfwd::extensions::discrete_mut_model m(std::move(__m));
    KTfwd::extensions::discrete_rec_model recmap(std::move(__recmap));
    metapop_t mpop2(std::move(*mpop));
    const auto recpos = KTfwd::extensions::bind_drm(recmap,mpop2.gametes,mpop2.mutations,
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
    for( size_t g = 0 ; g < simlen ; ++g, ++mpop2.generation )
      {
	std::vector<unsigned> Ns_next({unsigned(Nvector_A[g]),unsigned(Nvector_A[g])});
	unsigned N = std::accumulate(Ns_next.begin(),Ns_next.end(),0u);
	std::vector<double> wbars = sample_diploid(rng,
						   mpop2.gametes,
						   mpop2.diploids,
						   mpop2.mutations,
						   mpop2.mcounts,
						   &mpop2.Ns[0],
						   &Ns_next[0],
						   mu_tot,
						   KTfwd::extensions::bind_dmm(m,mpop2.mutations,mpop2.mut_lookup,rng,neutral,selected,mpop2.generation),
						   recpos,
						   fmodels,
						   std::bind(migpop,std::placeholders::_1,rng,migrate),
						   mpop2.neutral,
						   mpop2.selected,
						   &fs[0]);

	mpop2.Ns = Ns_next;
	//4*N b/c it needs to be fixed in the metapopulation
	KTfwd::update_mutations(mpop2.mutations,mpop2.fixations,mpop2.fixation_times,mpop2.mut_lookup,mpop2.mcounts,mpop2.generation,4*N);
      }
    gsl_rng_free(rng);
    *mpop=std::move(mpop2);
  }
  
  //Wow, that's too many args
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

