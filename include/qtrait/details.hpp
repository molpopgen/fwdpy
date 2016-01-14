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
    fwdpy::singlepop_t pop2(std::move(*pop));
    const auto recpos = KTfwd::extensions::bind_drm(recmap,pop2.gametes,pop2.mutations,
						    rng,recrate);

    //We use an empty fitness fxn here b/c the rules policies keep track of it separately.
    const auto ff = []( const fwdpy::singlepop_t::diploid_t &,
			const fwdpy::singlepop_t::gcont_t &,
			const fwdpy::singlepop_t::mcont_t ) noexcept { return 0.; };
    
    for( unsigned g = 0 ; g < simlen ; ++g, ++pop->generation )
      {
	const unsigned nextN = 	*(Nvector+g);
	KTfwd::experimental::sample_diploid(rng,
					    pop2.gametes,  
					    pop2.diploids, 
					    pop2.mutations,
					    pop2.mcounts,
					    pop2.N,
					    nextN,
					    mu_tot,
					    KTfwd::extensions::bind_dmm(m,pop2.mutations,pop2.mut_lookup,rng,neutral,selected,pop2.generation),
					    recpos,
					    ff,
					    pop2.neutral,pop2.selected,
					    f,
					    model_rules,
					    KTfwd::remove_nothing());
	KTfwd::update_mutations(pop2.mutations,pop2.mut_lookup,pop2.mcounts,2*nextN);
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
