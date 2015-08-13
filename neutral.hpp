#ifndef __FWDPP_CYTHON_SINGLEPOP_HPP__
#define __FWDPP_CYTHON_SINGLEPOP_HPP__

#include <fwdpp/diploid.hh>
#include <fwdpp/sugar.hpp> //lazy include of sugar library
#include <types.hpp>

namespace fwdpy {

  void evolve_pop(GSLrng_t * rng, singlepop_t * pop, const unsigned & ngens,const double & theta, const double & rho);
  //std::vector<int> sfs_from_sample(GSLrng_t * rng,const singlepop_t * pop,const unsigned & nsam);
}

#endif
