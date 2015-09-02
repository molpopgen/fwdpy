#include <sample.hpp>

#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <fwdpp/diploid.hh>

namespace fwdpy {
  
  std::vector<std::vector<std::pair<double,std::string> >> take_sample_from_pop(GSLrng_t * rng,const popvector * pops,const unsigned & nsam)
  {
    std::vector<std::vector<std::pair<double,std::string> > > rv;
    for(unsigned i=0;i<pops->pops.size();++i)
      {
	rv.emplace_back( KTfwd::ms_sample(rng->get(),&(pops->pops[i].get()->diploids),nsam,true));
      }
    return rv;
  }

  double tajd( const std::vector<std::pair<double,std::string> > & __data )
  {
    Sequence::SimData d(__data.begin(),__data.end());
    Sequence::PolySIM ad(&d);
    return ad.TajimasD();
  }

}
