#include <sample.hpp>
#include <fwdpp/diploid.hh>
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace fwdpy {
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
	a->push_back(double(gen-mitr->g)); //mutation age--this is correct b/c of def'n of 'gen' in the pop objects!
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

  void get_sh( const std::vector<std::pair<double,std::string> > & samples,
	       const metapop_t * pop,
	       std::vector<double> * s,
	       std::vector<double> * h,
	       std::vector<double> * p,
	       std::vector<double> * a)
  {
    unsigned ttlN=0;
    for(auto itr = pop->diploids.begin();itr!=pop->diploids.end();++itr)
      {
	ttlN+=itr->size();
      }
    get_sh_details(samples,
		   pop->mutations,
		   ttlN,
		   pop->generation,
		   s,h,p,a);
  }
}

