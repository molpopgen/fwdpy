#ifndef __FWDPY_SUMMSTATS_HPP__
#define __FWDPY_SUMMSTATS_HPP__

#include <vector>
#include <utility>
#include <string>
#include <map>
namespace fwdpy
{
  namespace libseq
  {
    std::map<std::string,double> libseq_basic_stats( const std::vector<std::pair<double,std::string> > & __data );
  }
}

#endif
