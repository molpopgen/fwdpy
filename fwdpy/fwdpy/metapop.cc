#include <gsl/gsl_randist.h>
#include <types.hpp>
#include <fwdpp/util.hpp>
#include <memory>
namespace fwdpy {

  size_t migpop(const size_t & source_pop, gsl_rng * r, const double & mig_prob)
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
    singlepop_t spop(*pop);
    mpop->diploids.clear();
    mpop->gametes.clear();
    mpop->Ns.clear();
    // //Move construct from spop into mpop
    mpop->fixations = std::move(spop.fixations);
    mpop->fixation_times = std::move(spop.fixation_times);
    mpop->mutations = std::move(spop.mutations);
    mpop->gametes = std::move(spop.gametes);
    mpop->diploids.emplace_back(std::move(spop.diploids));
    mpop->mut_lookup = std::move(spop.mut_lookup);
    mpop->Ns.push_back(unsigned(mpop->diploids[0].size()));
    //Update other data
    mpop->generation = spop.generation;
  }

  //Make a copy of deme i
  void copy_deme( metapop_t * mpop,
		  const size_t i,
		  const int update_counts)
  {
    mpop->diploids.push_back(metapop_t::dipvector_t(mpop->diploids[i]));
    mpop->Ns.push_back(unsigned(mpop->diploids[i].size()));
    size_t newpop_idx = mpop->diploids.size()-1;
    if( update_counts ) //normally not needed--the next "sample_diploid" call will do the trick
      {
	for( auto itr = mpop->diploids[newpop_idx].begin() ;
	     itr != mpop->diploids[newpop_idx].end() ; ++itr )
	  {
	    itr->first->n++;
	    itr->second->n++;
	  }
	for( auto & m : mpop->mutations) m.checked=false;
	for(auto itr = mpop->gametes.begin() ;
	    itr != mpop->gametes.end() ; ++itr )
	  {
	    KTfwd::adjust_mutation_counts(itr,itr->n);
	  }
      }
  }

}
