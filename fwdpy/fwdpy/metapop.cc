#include <memory>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include <fwdpp/util.hpp>
#include <fwdpp/sugar/demography.hpp>
#include "types.hpp"
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
    //use copy-construction from singlepop.
    *mpop=metapop_t(*pop);
  }

  /*
    Generic demographic operations.

    Calls are made to fwdpp's sugar layer, which in turn 
    calls lower-level functions.  

    Any non-zero return value is treated as a run-time error
    and an exception is thrown.
   */
  void copy_deme( metapop_t * mpop,
		  const size_t i )
  {
    int rv = KTfwd::copy_pop(*mpop,i);
    if(rv) throw std::runtime_error("error copying deme");
  }

  void remove_deme( metapop_t * mpop, const std::size_t i )
  {
    int rv=KTfwd::remove_pop(*mpop,i);
    if(rv) throw std::runtime_error("error removing deme");
  }
  
  void merge_demes(metapop_t  * mpop, const std::size_t i, const std::size_t j)
  {
    int rv=KTfwd::merge_pops(*mpop,i,j);
    if(rv) throw std::runtime_error("error merging demes");
  }
  
  void split_deme(const gsl_rng * r, metapop_t * mpop, const std::size_t i, const unsigned N_new, const bool replacement )
  {
    int rv=KTfwd::split_pop(r,*mpop,i,N_new,replacement);
    if(rv) throw std::runtime_error("error splitting deme");
  }
  
  void admix_demes(const gsl_rng * r, metapop_t * mpop, const std::size_t i, const std::size_t j, const double prop_i, const bool replacement)
  {
    int rv=KTfwd::admix_pops(r,*mpop,i,j,prop_i,replacement);
    if(rv) throw std::runtime_error("error admixing demes");
  }
}
