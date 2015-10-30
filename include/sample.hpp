#ifndef __FWDPY_SAMPLE_HPP__
#define __FWDPY_SAMPLE_HPP__

#include <types.hpp>

namespace fwdpy {
    void get_sh( const std::vector<std::pair<double,std::string> > & ms_sample,
	       const singlepop_t * pop,
	       std::vector<double> * s,
	       std::vector<double> * h,
	       std::vector<double> * p,
	       std::vector<double> * a);
  void get_sh( const std::vector<std::pair<double,std::string> > & samples,
	       const metapop_t * pop,
	       std::vector<double> * s,
	       std::vector<double> * h,
	       std::vector<double> * p,
	       std::vector<double> * a);
}

#endif
