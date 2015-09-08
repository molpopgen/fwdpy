#ifndef __FWDPP_CYTHON_SINGLEPOP_HPP__
#define __FWDPP_CYTHON_SINGLEPOP_HPP__

#include <fwdpp/diploid.hh>
#include <fwdpp/sugar.hpp> //lazy include of sugar library
#include <types.hpp>
#include <memory>

namespace fwdpy {

  void evolve_pop(GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,const std::vector<unsigned> & nlist,const double & theta, const double & rho);
}

#endif
