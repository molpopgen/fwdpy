#include <gsl/gsl_randist.h>

namespace fwdpy {

//Should move to a more general .cpp file in future
  size_t migpop(const size_t & source_pop, gsl_rng * r, const double & mig_prob)
  {
    if( gsl_rng_uniform(r) < mig_prob )
      {
	return ! source_pop;
      }
    return source_pop;
  }

}
