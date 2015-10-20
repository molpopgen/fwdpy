#include <map>
#include <limits>
#include <vector>
#include <functional>
#include <algorithm>
#include <stdexcept>
#include <Sequence/SimData.hpp>
#include <Sequence/SummStats/nSL.hpp>
#include <Sequence/SummStats/Garud.hpp>
#include <Sequence/SummStats/lHaf.hpp>

using namespace std;

namespace fwdpy
{
  namespace libseq
  {
    map<string,double> libseq_extra_ld_stats( const vector<pair<double,string> > & __data,
					      const double minfreq,
					      const double binsize,
					      const std::vector<double> & gmap)
    {
      if ( !is_sorted(__data.begin(),__data.end(),
		     []( const pair<double,string> & a,
			 const pair<double,string> & b ) {
		       return a.first<b.first;
		     } ))
	{
	  throw runtime_error("libseq_extra_ld_stats: input data not sorted in ascending order of variant position");
	}
      Sequence::SimData d(__data.begin(),__data.end());
      Sequence::PolySIM ad(&d);
      map<string,double> rv;

      auto gstats = Sequence::H1H12(d);
      auto nsl = Sequence::snSL(d,minfreq,binsize,gmap.empty() ? nullptr : &gmap[0]);
      rv["H1"]=gstats.H1;
      rv["H12"]=gstats.H12;
      rv["H2H1"]=gstats.H2H1;
      rv["mnsSL"]=nsl.first;
      rv["mnsiHs"]=nsl.second;
      return rv;
    }

    map<string,double> libseq_extra_ld_stats( const vector<pair<double,string> > & __data,
					      const double minfreq,
					      const double binsize)
    {
      return libseq_extra_ld_stats(__data,minfreq,binsize,vector<double>());
    }

    std::vector<double> lhaf( const std::vector<std::pair<double,std::string> > & __data, const double l )
    {
      Sequence::SimData d(__data.begin(),__data.end());
      return Sequence::lHaf(d,l);
    }
  } //ns libseq
} //ns fwdpy
