//Eyre-Walker 2010, backwards
#include <ew_backwards_rules.hpp>
#include <fwdpp/experimental/sample_diploid.hpp>
/*
  Problems that we'll have:
  We need a way to track effects sizes on trait value.

  Currently, fwdpp's popgenmut doesn't have any extra available members.

  Some hacks are possible:
  1. Usurp the 'h' term in popgenmut for the mutation effect size.
  2. Pass the 'h' for fitness on to the rules struct

  Or:

  1. Keep a totally different data type relating mutation to trait size

  Both have pros and cons...
*/

namespace fwdpy
{
  namespace qtrait
  {
    void qtrait_sim_details_t( gsl_rng * rng,
			       fwdpy::singlepop_t * pop,
			       const unsigned * Nvector,
			       const size_t Nvector_len,
			       const double & neutral,
			       const double & selected,
			       const double & recrate,
			       const double & f,
			       const bool track,  //do we want to track the trajectories of all mutations?
			       const ew_backwards_rules & model_rules)   
  {
    const unsigned simlen = Nvector_len;
    
    const double mu_tot = neutral + selected;
    
    std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng);
    for( unsigned g = 0 ; g < simlen ; ++g, ++pop->generation )
      {
	const unsigned nextN = 	*(Nvector+g);
	KTfwd::experimental::sample_diploid(rng,
					    &pop->gametes,  
					    &pop->diploids, 
					    &pop->mutations,
					    pop->N,
					    nextN,
					    mu_tot,
					    //NEED TO FIX MUTATION MODEL -- this is copied from diploid_ind.cc as placeholder
					    std::bind(KTfwd::infsites(),rng,&pop->mut_lookup,pop->generation,
						      neutral,0.,[&rng](){return gsl_rng_uniform(rng);},[](){return 0.;},[](){return 0.;}),
					    //The recombination policy includes the uniform crossover rate
					    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
						      std::placeholders::_3,
						      //Pass as reference
						      std::ref(pop->neutral),std::ref(pop->selected),
						      &pop->gametes,
						      recrate,
						      rng,
						      recmap),
					    std::bind(KTfwd::insert_at_end<typename fwdpy::singlepop_t::mutation_t,typename fwdpy::singlepop_t::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					    std::bind(KTfwd::insert_at_end<typename fwdpy::singlepop_t::gamete_t,typename fwdpy::singlepop_t::glist_t>,std::placeholders::_1,std::placeholders::_2),
					    //We use an empty fitness fxn here b/c the rules policies keep track of it separately.
					    []( typename fwdpy::singlepop_t::dipvector_t::const_iterator & ) { return 0.; },
					    KTfwd::remove_nothing(),
					    f,
					    model_rules);
	KTfwd::remove_lost(&pop->mutations,&pop->mut_lookup);
	//This being put here ignores any mutation existing for only 1 generation
	if(track) pop->updateTraj();
	assert(KTfwd::check_sum(pop->gametes,2*nextN));
      }
    //Update population's size variable to be the current pop size
    pop->N = pop->diploids.size();
  }
  }
}
