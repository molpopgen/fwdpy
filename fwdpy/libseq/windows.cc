#include <Sequence/SimData.hpp>
#include <Sequence/PolyTableSlice.hpp>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>

using namespace std;

namespace fwdpy
{
  namespace libseq
  {
    std::vector< std::vector<std::pair<double,std::string> > >
    sliding_windows_cpp( const std::vector<std::pair<double,std::string> > & sample,
			 const double window_size,
			 const double steplen,
			 const double starting_pos,
			 const double ending_pos)
    {
      std::vector< std::vector<std::pair<double,std::string> > > rv;
      if(sample.empty()) return rv;
      
      Sequence::SimData d(sample.begin(),sample.end());
      Sequence::PolyTableSlice<Sequence::SimData> s(d.sbegin(),d.send(),window_size,steplen,starting_pos,ending_pos);
      std::for_each(s.cbegin(),s.cend(),[&rv](const Sequence::PolyTableSlice<Sequence::SimData>::const_iterator::value_type & v) {
	  rv.push_back( std::vector<std::pair<double,std::string> >(v.first,v.second) );
	});
      return rv;
    }
  }
}
