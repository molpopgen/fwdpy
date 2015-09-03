#include <sample.hpp>

#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <fwdpp/diploid.hh>
#include <algorithm>

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

  void get_sh_details( const std::vector<std::pair<double,std::string> > & sample,
		       const singlepop_t::mlist_t & mutations,
		       std::vector<double> * s,
		       std::vector<double> * h)
  {
    std::for_each(sample.begin(),sample.end(),[&mutations,&s,&h](const std::pair<double,std::string> & p) {
	auto mitr = std::find_if(mutations.begin(),mutations.end(),[&p]( const singlepop_t::mutation_t & m ) {
	    return p.first==m.pos;
	  });
	s->push_back(mitr->s);
	h->push_back(mitr->h);
      });
  }
  
  void get_sh( const std::vector< std::vector<std::pair<double,std::string> > > & samples,
	       const popvector * pops, const unsigned i,
	       std::vector<double> * s,
	       std::vector<double> * h)
  {
    get_sh_details(samples[i],pops->pops[i].get()->mutations,s,h);
  }
}
