/*!
  \file sample.hpp
  \brief Functions related to taking samples from populations.
*/
#ifndef __FWDPY_SAMPLE_HPP__
#define __FWDPY_SAMPLE_HPP__

#include "types.hpp"

namespace fwdpy {
  /*!
    \brief Get detailed info about mutations in a sample taken from a fwdpy::singlepop_t
    \note Definition in fwdpy/fwdpy/sample.cc
  */
  void get_sh( const std::vector<std::pair<double,std::string> > & ms_sample,
	       const singlepop_t * pop,
	       std::vector<double> * s,
	       std::vector<double> * h,
	       std::vector<double> * p,
	       std::vector<double> * a);
  /*!
    \brief Get detailed info about mutations in a sample taken from a fwdpy::metapop_t
    \note Definition in fwdpy/fwdpy/sample.cc
  */
  void get_sh( const std::vector<std::pair<double,std::string> > & samples,
	       const metapop_t * pop,
	       std::vector<double> * s,
	       std::vector<double> * h,
	       std::vector<double> * p,
	       std::vector<double> * a);
}

#endif
