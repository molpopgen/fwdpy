#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>

namespace fwdpy {

  namespace libseq {
    double tajd( const std::vector<std::pair<double,std::string> > & __data )
    {
      Sequence::SimData d(__data.begin(),__data.end());
      Sequence::PolySIM ad(&d);
      return ad.TajimasD();
    }
  }
}
