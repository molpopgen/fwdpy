#ifndef __FWDPY_QTRAIT_DETAILS_HPP__
#define __FWDPY_QTRAIT_DETAILS_HPP__

#include "types.hpp"
#include <fwdpp/extensions/regions.hpp>

namespace fwdpy {
namespace qtrait {
  template<typename rules>
  void qtrait_sim_details_t( unsigned long seed,
			     fwdpy::singlepop_t * pop,
			     const unsigned * Nvector,
			     const size_t Nvector_len,
			     const double neutral,
			     const double selected,
			     const double recrate,
			     const double f,
			     const double sigmaE,
			     const double optimum ,
			     const bool track,  //do we want to track the trajectories of all mutations?
			     KTfwd::extensions::discrete_mut_model && __m,
			     KTfwd::extensions::discrete_rec_model && __recmap,
			     rules && __model_rules)   
  {
    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
    const unsigned simlen = Nvector_len;
    const double mu_tot = neutral + selected;

    KTfwd::extensions::discrete_mut_model m(std::move(__m));
    KTfwd::extensions::discrete_rec_model recmap(std::move(__recmap));
    rules model_rules(std::move(__model_rules));
    std::function<double(void)> recpos = std::bind(&KTfwd::extensions::discrete_rec_model::operator(),&recmap,rng);

    fwdpy::singlepop_t pop2(std::move(*pop));
    for( unsigned g = 0 ; g < simlen ; ++g, ++pop->generation )
      {
	const unsigned nextN = 	*(Nvector+g);
	KTfwd::experimental::sample_diploid(rng,
					    &pop2.gametes,  
					    &pop2.diploids, 
					    &pop2.mutations,
					    pop2.N,
					    nextN,
					    mu_tot,
					    std::bind(&KTfwd::extensions::discrete_mut_model::make_mut<decltype(pop2.mut_lookup)>,&m,rng,neutral,selected,pop2.generation,&pop2.mut_lookup),
					    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
						      std::ref(pop2.neutral),std::ref(pop2.selected),
						      &pop2.gametes,
						      recrate,
						      rng,
						      recpos),
					    std::bind(KTfwd::insert_at_end<typename fwdpy::singlepop_t::mutation_t,typename fwdpy::singlepop_t::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					    std::bind(KTfwd::insert_at_end<typename fwdpy::singlepop_t::gamete_t,typename fwdpy::singlepop_t::glist_t>,std::placeholders::_1,std::placeholders::_2),
					    //We use an empty fitness fxn here b/c the rules policies keep track of it separately.
					    []( typename fwdpy::singlepop_t::dipvector_t::const_iterator & ) { return 0.; },
					    KTfwd::remove_nothing(),
					    f,
					    model_rules);
	KTfwd::remove_lost(&pop2.mutations,&pop2.mut_lookup);
	//This being put here ignores any mutation existing for only 1 generation
	if(track) pop2.updateTraj();
	assert(KTfwd::check_sum(pop2.gametes,2*nextN));
      }
    gsl_rng_free(rng);
    //Update population's size variable to be the current pop size
    pop2.N = pop2.diploids.size();
    *pop = std::move(pop2);
  }
}
} //namespace

#endif
