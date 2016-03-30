/*!
  \file metapop.hpp

  \brief Operations on metapopulations.  Wrappers around fwdpp's demographic functions.
*/

#ifndef __FWDPY_METAPOP_HPP__
#define __FWDPY_METAPOP_HPP__

#include <gsl/gsl_randist.h>
#include <fwdpp/sugar/metapop.hpp>
#include "types.hpp"

namespace fwdpy
{
  size_t migpop(const size_t & source_pop, const gsl_rng * r, const double & mig_prob);

  void re_init_mpop( metapop_t * mpop, const singlepop_t * pop);

  /*
    Generic demographic operations.
    These functions are exposed to Cython as except + b/c their implementations
    will throw std::runtime_error if the underying fwdpp function returns a non-zero value
   */
  void copy_deme( metapop_t * mpop, const std::size_t i );
  void remove_deme( metapop_t * mpop, const std::size_t i );
  void merge_demes(metapop_t  * mpop, const std::size_t i, const std::size_t j);
  void split_deme(const gsl_rng * r, metapop_t * mpop, const std::size_t i, const unsigned N_new, const bool replacement );
  void admix_demes(const gsl_rng * r, metapop_t * mpop, const std::size_t i, const std::size_t j, const double prop_i, const bool replacement);
}

#endif
