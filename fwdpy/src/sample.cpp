#include <sample.hpp>

#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <fwdpp/diploid.hh>
#include <algorithm>

namespace fwdpy {
  
  std::vector<std::pair<double,std::string>> take_sample_from_pop(GSLrng_t * rng,const singlepop_t * pop,const unsigned & nsam)
  {
    return KTfwd::ms_sample(rng->get(),&(pop->diploids),nsam,true);
  }

  double tajd( const std::vector<std::pair<double,std::string> > & __data )
  {
    Sequence::SimData d(__data.begin(),__data.end());
    Sequence::PolySIM ad(&d);
    return ad.TajimasD();
  }

  void get_sh_details( const std::vector<std::pair<double,std::string> > & sample,
		       const singlepop_t::mlist_t & mutations,
		       const unsigned & twoN, const unsigned & gen,
		       std::vector<double> * s,
		       std::vector<double> * h,
		       std::vector<double> * p,
		       std::vector<double> * a)
  {
    std::for_each(sample.begin(),sample.end(),[&mutations,&s,&h,&p,&a,&twoN,&gen](const std::pair<double,std::string> & __pair) {
	auto mitr = std::find_if(mutations.begin(),mutations.end(),[&__pair]( const singlepop_t::mutation_t & m ) {
	    return __pair.first==m.pos;
	  });
	s->push_back(mitr->s);
	h->push_back(mitr->h);
	p->push_back(double(mitr->n)/double(twoN));
	a->push_back(double(gen-mitr->g)+1.); //mutation age
      });
  }
  
void get_sh( const std::vector<std::pair<double,std::string> > & samples,
	     const singlepop_t * pop,
	     std::vector<double> * s,
	     std::vector<double> * h,
	     std::vector<double> * p,
	     std::vector<double> * a)
  {
    get_sh_details(samples,
		   pop->mutations,
		   2*pop->diploids.size(),
		   pop->generation,
		   s,h,p,a);
  }
}
