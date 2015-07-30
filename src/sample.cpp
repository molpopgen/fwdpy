#include <sample.hpp>

#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <fwdpp/diploid.hh>

namespace fwdpy {
  
  std::vector<std::pair<double,std::string> > take_sample_from_pop(GSLrng_t * rng,const singlepop_t * pop,const unsigned & nsam)
  {
    return KTfwd::ms_sample(rng->get(),&pop->diploids,nsam,true);
  }

  double tajd( const std::vector<std::pair<double,std::string> > & __data )
  {
    Sequence::SimData d(__data.begin(),__data.end());
    Sequence::PolySIM ad(&d);
    return ad.TajimasD();
  }

}
