#ifndef __FWDPY_LIBSEQ_HPP__
#define __FWDPY_LIBSEQ_HPP__

#include <vector>
#include <utility>
#include <string>
#include <map>
namespace fwdpy
{
  namespace libseq
  {
    std::map<std::string,double> libseq_basic_stats( const std::vector<std::pair<double,std::string> > & __data );
    std::map<std::string,double> libseq_extra_ld_stats( const std::vector<std::pair<double,std::string> > & __data,
							const double minfreq,
							const double binsize,
							const std::vector<double> & gmap);
    std::map<std::string,double> libseq_extra_ld_stats( const std::vector<std::pair<double,std::string> > & __data,
							const double minfreq,
							const double binsize );
    std::vector<double> lhaf( const std::vector<std::pair<double,std::string> > & __data, const double l );
    std::vector< std::vector<std::pair<double,std::string> > >
    sliding_windows_cpp( const std::vector<std::pair<double,std::string> > & sample,
			 const double window_size,
			 const double steplen,
			 const double starting_pos,
			 const double ending_pos);
  }
}

#endif
