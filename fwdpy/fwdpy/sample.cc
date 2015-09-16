#include <sample.hpp>
#include <fwdpp/diploid.hh>
#include <algorithm>

namespace {
  enum class treat_neutral {ALL,NEUTRAL,SELECTED};
  void add_fixations( std::vector<std::pair<double,std::string>> * sample,
		      const fwdpy::singlepop_t::mvector_t & fixations,
		      const unsigned nsam,
		      const treat_neutral treat )
  {
    for( const auto & f : fixations)
      {
	if( treat == treat_neutral::ALL )
	  {
	    sample->emplace_back( std::make_pair(f.pos,std::string(nsam,'1')) );
	  }
	else if (treat == treat_neutral::NEUTRAL && f.neutral ) //only add neutral mutations
	  {
	    sample->emplace_back( std::make_pair(f.pos,std::string(nsam,'1')) );
	  }
	else if (treat == treat_neutral::SELECTED && !f.neutral ) //only add selected mutations
	  {
	    sample->emplace_back( std::make_pair(f.pos,std::string(nsam,'1')) );
	  }
      }
  }
}

namespace fwdpy {
  
  std::vector<std::pair<double,std::string>> take_sample_from_pop(GSLrng_t * rng,const singlepop_t * pop,const unsigned nsam, const int remove_fixed)
  {
    auto rv = KTfwd::ms_sample(rng->get(),&(pop->diploids),nsam,remove_fixed);
    if(! remove_fixed)
      {
	add_fixations(&rv,pop->fixations,nsam,treat_neutral::ALL);
      }
    return rv;
  }

  std::pair<std::vector<std::pair<double,std::string> >,
	    std::vector<std::pair<double,std::string> > >
  take_sample_from_pop_sep(GSLrng_t * rng,const singlepop_t * pop,const unsigned nsam, const int remove_fixed)
  {
    auto rv = KTfwd::ms_sample_separate(rng->get(),&(pop->diploids),nsam,remove_fixed);
    if(! remove_fixed)
      {
	add_fixations(&rv.first,pop->fixations,nsam,treat_neutral::NEUTRAL);
	add_fixations(&rv.second,pop->fixations,nsam,treat_neutral::SELECTED);
      }
    return rv;
  }

  std::pair<std::vector<std::pair<double,std::string> >,
	    std::vector<std::pair<double,std::string> > >
  take_sample_from_metapop_sep(GSLrng_t * rng,const metapop_t * mpop,const unsigned & nsam, const int remove_fixed, const int deme)
  {
    auto temp = KTfwd::ms_sample_separate(rng->get(),&(mpop->diploids[deme]),nsam,remove_fixed);
    if(! remove_fixed)
      {
	add_fixations(&temp.first,mpop->fixations,nsam,treat_neutral::NEUTRAL);
	add_fixations(&temp.second,mpop->fixations,nsam,treat_neutral::SELECTED);
      }
    return temp;
  }

  std::pair< std::vector<std::pair<double,std::string > >,
	     std::vector<std::pair<double,std::string> > >
  sample_specific_diploids(const singlepop_t * pop, const std::vector<unsigned> & indlist, const int remove_fixed)
  {
    auto rv = KTfwd::fwdpp_internal::ms_sample_separate_single_deme(&(pop->diploids),indlist,2*indlist.size(),remove_fixed);
    if(!remove_fixed)
      {
	add_fixations(&rv.first,pop->fixations,2*indlist.size(),treat_neutral::NEUTRAL);
	add_fixations(&rv.first,pop->fixations,2*indlist.size(),treat_neutral::SELECTED);
      }
    return rv;
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

