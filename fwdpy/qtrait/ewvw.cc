/*
  Eyre-Walker 2010, but trait does not 
  contribute to 100% of variation in fitness.

  Will this really be any different from what is in qtrait_impl.cc?

  Is this not mathematically equivalent to increasing V(S), the 
  variance in the Gaussian fitness fxn?

  TODO: It isn't equivalent.  Need to make an ipynb on why not...
*/

#include <thread>
#include <qtrait/ewvw_rules.hpp>
#include <fwdpp/experimental/sample_diploid.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/extensions/regions.hpp>


using namespace std;
//using esize_lookup = vector<pair<double,double> >;

namespace fwdpy
{
  namespace qtrait
  {
    // inline KTfwd::popgenmut ewbw_mut_model(gsl_rng * r,
    // 					   singlepop_t::lookup_table_t * lookup,
    // 					   esize_lookup * esizes,
    // 					   const unsigned & generation,
    // 					   const double & neutral,
    // 					   const double & selected,
    // 					   const double & h,
    // 					   const double & tau,
    // 					   const double & sigma)
    // {
    //   double pos = gsl_rng_uniform(r);
    //   while(lookup->find(pos)!=lookup->end())
    // 	{
    // 	  pos = gsl_rng_uniform(r);
    // 	}
    //   lookup->insert(pos);

    //   if( gsl_rng_uniform(r) < neutral/(neutral+selected) )
    // 	{
    // 	  //no effect on trait value
    // 	  return KTfwd::popgenmut(pos,0.,0.,generation,1);
    // 	}
    //   //mutation affects trait value and has selection coefficient
    //   //effect size for mutation
    //   double z = gsl_ran_gaussian_ziggurat(r,1.);
    //   //update esizes
    //   auto itr = std::lower_bound(esizes->begin(),esizes->end(),pos,
    // 				  [](const esize_lookup::value_type & p,
    // 				     const double & __val)
    // 				  {
    // 				    return p.first < __val;
    // 				  });
    //   esizes->emplace(itr,make_pair(pos,z));
    //   double epsilon = gsl_ran_gaussian_ziggurat(r,sigma);
    //   double temp = fabs(z/(1.+epsilon));
    //   // s = exp(log(temp)/tau)
    //   return KTfwd::popgenmut(pos,exp(log(temp)/tau),h,generation,1);
    // }

    void ewvw_sim_details_t( const gsl_rng * rng,
			     fwdpy::singlepop_t * pop,
			     const unsigned * Nvector,
			     const size_t Nvector_len,
			     const double & neutral,
			     const double & selected,
			     const double & recrate,
			     const double & f,
			     const bool track,  //do we want to track the trajectories of all mutations?
			     const KTfwd::extensions::discrete_mut_model & m,
			     const KTfwd::extensions::discrete_rec_model & recmap,
			     ewvw_rules & model_rules)
    {
      const unsigned simlen = Nvector_len;
    
      const double mu_tot = neutral + selected;

      const auto recpos = KTfwd::extensions::bind_drm(recmap,pop->gametes,pop->mutations,
						      rng,recrate);
      
      //We use an empty fitness fxn here b/c the rules policies keep track of it separately.
      const auto ff = []( const fwdpy::singlepop_t::diploid_t &,
			  const fwdpy::singlepop_t::gcont_t &,
			  const fwdpy::singlepop_t::mcont_t ) noexcept { return 0.; };
      for( unsigned g = 0 ; g < simlen ; ++g, ++pop->generation )
	{
	  const unsigned nextN = *(Nvector+g);
	  KTfwd::experimental::sample_diploid(rng,
					      pop->gametes,  
					      pop->diploids, 
					      pop->mutations,
					      pop->mcounts,
					      pop->N,
					      nextN,
					      mu_tot,
					      KTfwd::extensions::bind_dmm(m,pop->mutations,pop->mut_lookup,rng,neutral,selected,pop->generation),
					      recpos,
					      ff,
					      pop->neutral,pop->selected,
					      f,
					      model_rules,
					      KTfwd::remove_nothing());
	  KTfwd::update_mutations(pop->mutations,pop->mut_lookup,pop->mcounts,2*nextN);
	  //This being put here ignores any mutation existing for only 1 generation
	  assert(KTfwd::check_sum(pop->gametes,2*nextN));
	}
      //Update population's size variable to be the current pop size
      pop->N = pop->diploids.size();
    }

    void evolve_ewvw_t( GSLrng_t * rng,
			std::vector<std::shared_ptr<singlepop_t> > * pops,
			const unsigned * Nvector,
			const size_t Nvector_length,
			const double mu_neutral,
			const double mu_selected,
			const double littler,
			const double f,
			const double sigmaE,
			const double VS_total,
			const double optimum,
			const int track,
			const std::vector<double> & nbegs,
			const std::vector<double> & nends,
			const std::vector<double> & nweights,
			const std::vector<double> & sbegs,
			const std::vector<double> & sends,
			const std::vector<double> & sweights,
			const std::vector<KTfwd::extensions::shmodel> * callbacks,
			const std::vector<double> & rbeg,
			const std::vector<double> & rend,
			const std::vector<double> & rweight)
    {
      const KTfwd::extensions::discrete_mut_model m(nbegs,nends,nweights,sbegs,sends,sweights,*callbacks);
      auto recmap = KTfwd::extensions::discrete_rec_model(rbeg,rend,rweight);
      std::vector<GSLrng_t> rngs;
      std::vector<ewvw_rules> rules;
      for(unsigned i=0;i<pops->size();++i)
	{
	  //Give each thread a new RNG + seed
	  rngs.emplace_back(GSLrng_t(gsl_rng_get(rng->get())) );
	  rules.emplace_back(ewvw_rules(VS_total,sigmaE,optimum,*std::max_element(Nvector,Nvector+Nvector_length)));
	}
      std::vector<std::thread> threads(pops->size());
      for(unsigned i=0;i<pops->size();++i)
	{
    	  threads[i]=std::thread(ewvw_sim_details_t,
	  			 rngs[i].get(),
	  			 pops->operator[](i).get(),
	  			 Nvector,Nvector_length,
	  			 mu_neutral,mu_selected,littler,f,track,
	  			 std::cref(m),std::cref(recmap),
	  			 std::ref(rules[i]));
	}
      for(unsigned i=0;i<threads.size();++i) threads[i].join();
    }
  }
}
