#include <thread>
#include <functional>
#include <fwdpp/diploid.hh>
#include <fwdpp/extensions/regions.hpp>
#include <types.hpp>
#include <metapop.hpp>
#include <evolve_regions.hpp>

namespace fwdpy {

  void evolve_regions_details( fwdpy::singlepop_t * pop,
			       gsl_rng * rng,
			       const unsigned * Nvector,
			       const size_t Nvector_len,
			       const double & neutral,
			       const double & selected,
			       const double & recrate,
			       const double & f,
			       const int track,
			       const char * fitness,
			       const KTfwd::extensions::discrete_mut_model & m,
			       const KTfwd::extensions::discrete_rec_model & recmap)
  {    
    const unsigned simlen = Nvector_len;
    
    const double mu_tot = neutral + selected;
    
    //Recombination policy: more complex than the standard case...
    std::function<double(void)> recpos = std::bind(&KTfwd::extensions::discrete_rec_model::operator(),&recmap,rng);

    //The fitness model
    std::function<double(const fwdpy::singlepop_t::dipvector_t::iterator &)> dipfit = std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,2.);
    if( std::string(fitness) == "additive" )
      {
     	dipfit = std::bind(KTfwd::additive_diploid(),std::placeholders::_1,2.);
      }
    for( unsigned g = 0 ; g < simlen ; ++g, ++pop->generation )
      {
	const unsigned nextN = 	*(Nvector+g);
	KTfwd::sample_diploid(rng,
			      &pop->gametes,  
			      &pop->diploids, 
			      &pop->mutations,
			      pop->N,
			      nextN,
			      mu_tot,
			      std::bind(&KTfwd::extensions::discrete_mut_model::make_mut<decltype(pop->mut_lookup)>,&m,rng,neutral,selected,pop->generation,&pop->mut_lookup),
			      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
					std::ref(pop->neutral),std::ref(pop->selected),
					&pop->gametes,
					recrate,
					rng,
					recpos),
			      std::bind(KTfwd::insert_at_end<fwdpy::singlepop_t::mutation_t,fwdpy::singlepop_t::mlist_t>,std::placeholders::_1,std::placeholders::_2),
			      std::bind(KTfwd::insert_at_end<fwdpy::singlepop_t::gamete_t,fwdpy::singlepop_t::glist_t>,std::placeholders::_1,std::placeholders::_2),
			      dipfit,
			      std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*nextN),
			      f);
	if (track) pop->updateTraj();
	pop->N=nextN;
	KTfwd::remove_fixed_lost(&pop->mutations,&pop->fixations,&pop->fixation_times,&pop->mut_lookup,pop->generation,2*nextN);
	assert(KTfwd::check_sum(pop->gametes,2*nextN));
      }
    //Update population's size variable to be the current pop size
    pop->N = pop->diploids.size();
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
			 // const std::vector<double> & nbegs,
			 // const std::vector<double> & nends,
			 // const std::vector<double> & nweights,
			 // const std::vector<double> & sbegs,
			 // const std::vector<double> & sends,
			 // const std::vector<double> & sweights,
			 // const std::vector<KTfwd::extensions::shmodel> * callbacks,
			 // const std::vector<double> & rbeg,
			 // const std::vector<double> & rend,
			 // const std::vector<double> & rweight,
			 const char * fitness)
  {
    const KTfwd::extensions::discrete_mut_model m(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks);
    auto recmap = KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw);
    std::vector<GSLrng_t> rngs;
    for(unsigned i=0;i<pops->size();++i)
      {
	//Give each thread a new RNG + seed
	rngs.emplace_back(GSLrng_t(gsl_rng_get(rng->get())) );
      }
    std::vector<std::thread> threads(pops->size());
    for(unsigned i=0;i<pops->size();++i)
      {
	threads[i]=std::thread(evolve_regions_details,pops->operator[](i).get(),rngs[i].get(),Nvector,Nvector_len,
			       mu_neutral,mu_selected,littler,f,track,fitness,std::cref(m),std::cref(recmap));
      }
    for(unsigned i=0;i<threads.size();++i) threads[i].join();
  }

  void split_and_evolve_details(metapop_t * mpop, 
				gsl_rng * rng,
				const unsigned * Nvector_A,
				const size_t Nvector_A_len,
				const unsigned * Nvector_B,
				const size_t Nvector_B_len,
				const double & neutral,
				const double & selected,
				const double & recrate,
				const std::vector<double> & fs,
				const KTfwd::extensions::discrete_mut_model & m,
				const KTfwd::extensions::discrete_rec_model & recmap,
				const char * fitness)
  {  
    const unsigned simlen = Nvector_A_len;
    const double mu_tot = neutral + selected;
    
    std::function<double(void)> recpos = std::bind(&KTfwd::extensions::discrete_rec_model::operator(),&recmap,rng);
    
    //The fitness model
    std::function<double(const fwdpy::singlepop_t::dipvector_t::iterator &)> dipfit = std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,2.);
    if( std::string(fitness) == "additive" )
      {
     	dipfit = std::bind(KTfwd::additive_diploid(),std::placeholders::_1,2.);
      }
    
    //vector of fitness models
    std::vector<decltype(dipfit)> fmodels({dipfit,dipfit});

    //FIX later
    const double migrate = 0.;
    
    //Finally, we can evolve this thing
    for( unsigned g = 0 ; g < simlen ; ++g, ++mpop->generation )
      {
	std::vector<unsigned> Ns_next({unsigned(Nvector_A[g]),unsigned(Nvector_A[g])});
	unsigned N = std::accumulate(Ns_next.begin(),Ns_next.end(),0u);
	std::vector<double> wbars = sample_diploid(rng,
						   &mpop->gametes,
						   &mpop->diploids,
						   &mpop->mutations,
						   &mpop->Ns[0],
						   &Ns_next[0],
						   mu_tot,
						   std::bind(&KTfwd::extensions::discrete_mut_model::make_mut<decltype(mpop->mut_lookup)>,&m,rng,neutral,selected,mpop->generation,&mpop->mut_lookup),
						   std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
							     std::ref(mpop->neutral),std::ref(mpop->selected),
							     &mpop->gametes,
							     recrate,
							     rng,
							     recpos),
						   std::bind(KTfwd::insert_at_end<metapop_t::mutation_t,metapop_t::mlist_t>,std::placeholders::_1,std::placeholders::_2),
						   std::bind(KTfwd::insert_at_end<metapop_t::gamete_t,metapop_t::glist_t>,std::placeholders::_1,std::placeholders::_2),
						   fmodels,
						   //4*N b/c it needs to be fixed in the metapopulation
						   std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,4*N),
						   std::bind(migpop,std::placeholders::_1,rng,migrate),
						   &fs[0]);
	mpop->Ns = Ns_next;
	//4*N b/c it needs to be fixed in the metapopulation
	KTfwd::remove_fixed_lost(&mpop->mutations,&mpop->fixations,&mpop->fixation_times,&mpop->mut_lookup,mpop->generation,4*N);
      }
  }
  
  //Wow, that's too many args
  void split_and_evolve_t(GSLrng_t * rng,
			  std::vector<std::shared_ptr<metapop_t> > * mpops,
			  const unsigned * Nvector_A,
			  const size_t Nvector_A_len,
			  const unsigned * Nvector_B,
			  const size_t Nvector_B_len,
			  const double & neutral,
			  const double & selected,
			  const double & recrate,
			  const std::vector<double> & fs,
			  const fwdpy::internal::region_manager * rm,
			  const char * fitness)
  {
    const KTfwd::extensions::discrete_mut_model m(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks);
    auto recmap = KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw);
    std::vector<GSLrng_t> rngs;
    for(unsigned i=0;i<mpops->size();++i)
      {
	//Give each thread a new RNG + seed
	rngs.emplace_back(GSLrng_t(gsl_rng_get(rng->get())) );
      }
    std::vector<std::thread> threads(mpops->size());
    for(unsigned i=0;i<mpops->size();++i)
      {
  	threads[i]=std::thread(split_and_evolve_details,
			       mpops->operator[](i).get(),
			       rngs[i].get(),
			       Nvector_A,
			       Nvector_A_len,
			       Nvector_B,
			       Nvector_B_len,
			       neutral,
			       selected,
			       recrate,
			       fs,
			       std::cref(m),std::cref(recmap),
			       fitness);
      }
    for(unsigned i=0;i<threads.size();++i) threads[i].join();
  }
}

