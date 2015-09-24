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

  I think I like the following:

  vector<pair<double,double> >

  we'll keep it sorted and use log-time search
*/

using namespace std;
using esize_lookup = vector<pair<double,double> >;

namespace fwdpy
{
  namespace qtrait
  {
    inline KTfwd::popgenmut ewbw_mut_model(gsl_rng * r,
					   singlepop_t::lookup_table_t * lookup,
					   esize_lookup * esizes,
					   const unsigned & generation,
					   const double & neutral,
					   const double & selected,
					   const double & h,
					   const double & tau,
					   const double & sigma)
    {
      double pos = gsl_rng_uniform(r);
      while(lookup->find(pos)!=lookup->end())
	{
	  pos = gsl_rng_uniform(r);
	}
      lookup->insert(pos);

      if( gsl_rng_uniform(r) < neutral/(neutral+selected) )
	{
	  //no effect on trait value
	  return KTfwd::popgenmut(pos,0.,0.,generation,1);
	}
      //mutation affects trait value and has selection coefficient
      //effect size for mutation
      double z = gsl_ran_gaussian_ziggurat(r,1.);
      //update esizes
      auto itr = std::lower_bound(esizes->begin(),esizes->end(),pos,
				  [](const esize_lookup::value_type & p,
				     const double & __val)
				  {
				    return p.first < __val;
				  });
      esizes->emplace(itr,make_pair(pos,z));
      double epsilon = gsl_ran_gaussian_ziggurat(r,sigma);
      double temp = fabs(z/(1.+epsilon));
      // s = exp(log(temp)/tau)
      return KTfwd::popgenmut(pos,exp(log(temp)/tau),h,generation,1);
    }
    
    void qtrait_sim_details_t( gsl_rng * rng,
			       fwdpy::singlepop_t * pop,
			       const unsigned * Nvector,
			       const size_t Nvector_len,
			       const double & neutral,
			       const double & selected,
			       const double & recrate,
			       const double & tau,
			       const double & sigma,
			       const double & h,
			       const double & f,
			       const bool track,  //do we want to track the trajectories of all mutations?
			       const ew_backwards_rules & model_rules,
			       esize_lookup * esizes) //this is our "hack"
  {
    const unsigned simlen = Nvector_len;
    
    const double mu_tot = neutral + selected;
    
    function<double(void)> recmap = bind(gsl_rng_uniform,rng);
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
					    bind(ewbw_mut_model,rng,&pop->mut_lookup,esizes,pop->generation,neutral,selected,h,tau,sigma),
					    //The recombination policy includes the uniform crossover rate
					    bind(KTfwd::genetics101(),placeholders::_1,placeholders::_2,
						      placeholders::_3,
						      //Pass as reference
						      ref(pop->neutral),ref(pop->selected),
						      &pop->gametes,
						      recrate,
						      rng,
						      recmap),
					    bind(KTfwd::insert_at_end<typename fwdpy::singlepop_t::mutation_t,typename fwdpy::singlepop_t::mlist_t>,placeholders::_1,placeholders::_2),
					    bind(KTfwd::insert_at_end<typename fwdpy::singlepop_t::gamete_t,typename fwdpy::singlepop_t::glist_t>,placeholders::_1,placeholders::_2),
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
