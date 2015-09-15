#ifndef __FWDPY_SAMPLE_HPP__
#define __FWDPY_SAMPLE_HPP__

#include <types.hpp>

namespace fwdpy {
  std::vector<std::pair<double,std::string> > take_sample_from_pop(GSLrng_t * rng,const singlepop_t * pop,const unsigned nsam,const int remove_fixed);
  std::pair<std::vector<std::pair<double,std::string> >,
			std::vector<std::pair<double,std::string> > >
  take_sample_from_pop_sep(GSLrng_t * rng,const singlepop_t * pop,const unsigned nsam, const int remove_fixed);
  std::pair<std::vector<std::pair<double,std::string> >,
	    std::vector<std::pair<double,std::string> > >
  take_sample_from_metapop_sep(GSLrng_t * rng,const metapop_t * mpop,const unsigned & nsam, const int remove_fixed, const int deme);
  std::pair< std::vector<std::pair<double,std::string> >,
	     std::vector<std::pair<double,std::string> > >
  sample_specific_diploids(const singlepop_t * pop, const std::vector<unsigned> & indlist, const int remove_fixed);
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
  double tajd( const std::vector<std::pair<double,std::string> > & __data );
}

#endif
