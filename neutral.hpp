#ifndef __FWDPP_CYTHON_SINGLEPOP_HPP__
#define __FWDPP_CYTHON_SINGLEPOP_HPP__

#include <fwdpp/diploid.hh>
#include <fwdpp/sugar.hpp> //lazy include of sugar library
#include <types.hpp>

namespace fwdpy {

  void evolve_pop(GSLrng_t * rng, popvector * pops,const unsigned & ngens,const double & theta, const double & rho);
}

#endif
