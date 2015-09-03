#ifndef __FWDPY_SAMPLE_HPP__
#define __FWDPY_SAMPLE_HPP__

#include <types.hpp>

namespace fwdpy {
  std::vector<std::vector<std::pair<double,std::string> >> take_sample_from_pop(GSLrng_t * rng,const popvector * pops,const unsigned & nsam);
  void get_sh( const std::vector< std::vector<std::pair<double,std::string> > >& samples,
	       const popvector * pops, const unsigned i,
	       std::vector<double> * s,
	       std::vector<double> * h,
	       std::vector<double> * p,
	       std::vector<double> * a);
  double tajd( const std::vector<std::pair<double,std::string> > & __data );
}

#endif
