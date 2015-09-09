#include <fwdpp/diploid.hh>
#include <fwdpp/extensions/regions.hpp>
#include <types.hpp>
#include <thread>

#include <functional>

namespace fwdpy {

  void evolve_regions_details( fwdpy::singlepop_t * pop,
			       gsl_rng * rng,
			       const unsigned * Nvector,
			       const size_t Nvector_len,
			       const double & neutral,
			       const double & selected,
			       const double & recrate,
			       const double & f,
			       const char * fitness,
			       const KTfwd::extensions::discrete_mut_model & m,
			       const KTfwd::extensions::discrete_rec_model & recmap)
  {    
    const unsigned simlen = Nvector_len;
    
    const double mu_tot = neutral + selected;
    
    //Recombination policy: more complex than the standard case...
    std::function<double(void)> recpos = std::bind(&KTfwd::extensions::discrete_rec_model::operator(),&recmap,rng);
    
    //The fitness model
    std::function<double(const fwdpy::singlepop_t::glist_t::const_iterator &,
			 const fwdpy::singlepop_t::glist_t::const_iterator &)> dipfit = std::bind(KTfwd::additive_diploid(),std::placeholders::_1,std::placeholders::_2,2.);
    if( std::string(fitness) == "multiplicative" )
      {
	dipfit = std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.);
      }
    for( unsigned g = 0 ; g < simlen ; ++g, ++pop->generation )
      {
	const unsigned nextN = 	*(Nvector+g);
	KTfwd::sample_diploid(rng,
			      &pop->gametes,  
			      &pop->diploids, 
			      &pop->mutations,
			      pop->N,
			      pop->N,
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
			      std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*pop->N),
			      f);
	pop->N=nextN;
	KTfwd::remove_fixed_lost(&pop->mutations,&pop->fixations,&pop->fixation_times,&pop->mut_lookup,pop->generation,2*nextN);
	assert(KTfwd::check_sum(pop->gametes,2*pop->N));
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
			 const std::vector<double> & nbegs,
			 const std::vector<double> & nends,
			 const std::vector<double> & nweights,
			 const std::vector<double> & sbegs,
			 const std::vector<double> & sends,
			 const std::vector<double> & sweights,
			 const std::vector<KTfwd::extensions::shmodel> * callbacks,
			 const std::vector<double> & rbeg,
			 const std::vector<double> & rend,
			 const std::vector<double> & rweight,
			 const char * fitness)
  {
    const KTfwd::extensions::discrete_mut_model m(nbegs,nends,nweights,sbegs,sends,sweights,*callbacks);
    auto recmap = KTfwd::extensions::discrete_rec_model(rbeg,rend,rweight);
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
			       mu_neutral,mu_selected,littler,f,fitness,std::cref(m),std::cref(recmap));
      }
    for(unsigned i=0;i<threads.size();++i) threads[i].join();
  }
}

