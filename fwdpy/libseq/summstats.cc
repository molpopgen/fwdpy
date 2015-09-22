#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <map>
#include <string>
#include <algorithm>
#include <stdexcept>
using namespace std;

namespace fwdpy
{
  namespace libseq
  {
    map<string,double> libseq_basic_stats( const vector<pair<double,string> > & __data )
    {
      if(! is_sorted(__data.cbegin(),__data.cend(),
		     []( const pair<double,string> & a,const pair<double,string> & b ) {
		       return a.first < b.first;
		     }))
	{
	  throw runtime_error("data not sorted");
	}
      Sequence::SimData d(__data.begin(),__data.end());
      Sequence::PolySIM ad(&d);
      map<string,double> rv;

      rv["S"] = ad.NumPoly();
      rv["singletons"] = ad.NumSingletons();
      rv["dsingletons"] = ad.NumExternalMutations();
      rv["thetaw"] = ad.ThetaW();
      rv["thetapi"] = ad.ThetaPi();
      rv["tajd"] = ad.TajimasD();
      rv["thetah"] = ad.ThetaH();
      rv["hprime"] = ad.Hprime();
      rv["tajd"] = ad.TajimasD();
      return rv;
    }
  }
}
