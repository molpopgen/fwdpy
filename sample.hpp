#ifndef __FWDPY_SAMPLE_HPP__
#define __FWDPY_SAMPLE_HPP__

#include <types.hpp>

namespace fwdpy {
  std::vector<std::pair<double,std::string> > take_sample_from_pop(GSLrng_t * rng,const singlepop_t * pop,const unsigned & nsam);
  double tajd( const std::vector<std::pair<double,std::string> > & __data );
}

#endif
