#ifndef __FWDPY_METAPOP_HPP__
#define __FWDPY_METAPOP_HPP__

#include <gsl/gsl_randist.h>

namespace fwdpy
{
  size_t migpop(const size_t & source_pop, gsl_rng * r, const double & mig_prob);
}

#endif
