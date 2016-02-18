#ifndef __FWDPY_METAPOP_HPP__
#define __FWDPY_METAPOP_HPP__

#include <gsl/gsl_randist.h>
#include "types.hpp"

namespace fwdpy
{
  size_t migpop(const size_t & source_pop, const gsl_rng * r, const double & mig_prob);

  void re_init_mpop( metapop_t * mpop, const singlepop_t * pop);

  void copy_deme( metapop_t * mpop,
		  const size_t i );
}

#endif
