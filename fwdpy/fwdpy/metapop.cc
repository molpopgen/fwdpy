#include <gsl/gsl_randist.h>
#include <types.hpp>
#include <fwdpp/util.hpp>
#include <memory>
namespace fwdpy {

  size_t migpop(const size_t & source_pop, const gsl_rng * r, const double & mig_prob)
  {
    if( gsl_rng_uniform(r) < mig_prob )
      {
	return ! source_pop;
      }
    return source_pop;
  }

  void re_init_mpop( metapop_t * mpop,
		     const singlepop_t * pop)
  {
    mpop->diploids.clear();
    mpop->gametes.clear();
    mpop->Ns.clear();

    mpop->fixations = pop->fixations;
    mpop->fixation_times = pop->fixation_times;
    mpop->mutations = pop->mutations;
    mpop->gametes = pop->gametes;
    mpop->diploids.emplace_back(pop->diploids);
    mpop->mut_lookup = pop->mut_lookup;
    mpop->Ns.push_back(unsigned(mpop->diploids[0].size()));
    mpop->mcounts = pop->mcounts;
    //Update other data
    mpop->generation = pop->generation;
  }

  //Make a copy of deme i
  void copy_deme( metapop_t * mpop,
		  const size_t i )
  {
    mpop->diploids.push_back(metapop_t::dipvector_t(mpop->diploids[i]));
    mpop->Ns.push_back(unsigned(mpop->diploids[i].size()));
    KTfwd::fwdpp_internal::process_gametes(mpop->gametes,mpop->mutations,mpop->mcounts);
  }

}
