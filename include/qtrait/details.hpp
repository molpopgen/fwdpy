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
			       const int track,  //do we want to track the trajectories of all mutations and how often?
			       const int trackStats,  //do we want to track VG, etc., and how often?
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
      const auto recpos = KTfwd::extensions::bind_drm(recmap,pop->gametes,pop->mutations,
						      rng,recrate);

      //We use an empty fitness fxn here b/c the rules policies keep track of it separately.
      const auto ff = []( const fwdpy::singlepop_t::diploid_t &,
			  const fwdpy::singlepop_t::gcont_t &,
			  const fwdpy::singlepop_t::mcont_t ) noexcept { return 0.; };

      for( unsigned g = 0 ; g < simlen ; ++g, ++pop->generation )
	{
	  if(track&&pop->generation&&pop->generation%track==0.) pop->updateTraj();
	  if(trackStats&&pop->generation&&pop->generation%trackStats==0) pop->updateStats();
	  const unsigned nextN = 	*(Nvector+g);
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
	  assert(KTfwd::check_sum(pop->gametes,2*nextN));
	}
      //make sure we update in the last generation if needed
      if(track&&pop->generation&&pop->generation%track==0.) pop->updateTraj();
      if(trackStats&&pop->generation&&pop->generation%trackStats==0) pop->updateStats();
      gsl_rng_free(rng);
      //Update population's size variable to be the current pop size
      pop->N = pop->diploids.size();
    }
  }
} //namespace

#endif
