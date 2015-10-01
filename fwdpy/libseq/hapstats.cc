#include <map>
#include <limits>
#include <vector>
#include <functional>
#include <Sequence/SimData.hpp>
#include <Sequence/SummStats/nSL.hpp>
#include <Sequence/SummStats/Garud.hpp>
using namespace std;

namespace fwdpy
{
  namespace libseq
  {
    map<string,double> libseq_extra_ld_stats( const vector<pair<double,string> > & __data,
					      const double minfreq,
					      const double binsize,
					      const double * gmap)
    {
      Sequence::SimData d(__data.begin(),__data.end());
      Sequence::PolySIM ad(&d);
      map<string,double> rv;

      auto gstats = Sequence::H1H12(d);
      auto nsl = Sequence::snSL(d,minfreq,binsize,gmap);
      rv["H1"]=gstats.H1;
      rv["H12"]=gstats.H12;
      rv["H2H1"]=gstats.H2H1;
      rv["mnsSL"]=nsl.first;
      rv["mnsiHs"]=nsl.second;
      return rv;
    }
  } //ns libseq
} //ns fwdpy
