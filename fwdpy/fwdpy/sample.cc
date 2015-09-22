#include <sample.hpp>
#include <fwdpp/diploid.hh>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <Sequence/PolyTableSlice.hpp>

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

  void fill_diploid_view_details( std::map<std::string, std::vector<double> > & rv,
				  //that's a mouthful!
				  const fwdpy::singlepop_t::gamete_t::mutation_container & muts,
				  const size_t ind,
				  const unsigned twoN,
				  const int remove_fixed,
				  const int chrom)
  {
    for( const auto itr : muts )
      {
	//is mutation fixed or not?
	double p = double(itr->n)/double(twoN);
	if( p < 1. || !remove_fixed )
	  {
	    rv["ind"].push_back(double(ind));
	    rv["pos"].push_back(itr->pos);
	    rv["esize"].push_back(itr->s);
	    rv["h"].push_back(itr->h);
	    rv["p"].push_back(p);
	    rv["hap"].push_back(double(chrom));
	    rv["origin"].push_back(double(itr->g));
	  }
      }
  }
  /*
    This will work for singlepop_t and metapop_t dipvectors, as they instantiate to the same type.
    Can throw std::out_of_range
  */
  std::map<std::string, std::vector<double> > diploid_view_details(const fwdpy::singlepop_t::dipvector_t & diploids,
								   const size_t ind,
								   const int remove_fixed)
  {
    if(ind >= diploids.size())
      {
	throw std::out_of_range("diploid_view_details: individual index out of range");
      }
    const unsigned twoN = 2*diploids.size();

    using vd = std::vector<double>;
    std::map<std::string, vd > rv;
    std::vector<std::string> rvlabels = {"ind","pos","esize","h","p","hap","origin"};
    for (const auto & l : rvlabels)
      {
	rv[l]=vd();
      }
    if( diploids[ind].first->mutations.empty() &&
	diploids[ind].first->smutations.empty() &&
	diploids[ind].second->mutations.empty() &&
	diploids[ind].second->smutations.empty() )
      {
	//diploid has no mutations!
	rv["ind"].push_back(double(ind));
	for (const auto & l : rvlabels)
	  {
	    if(l!="ind")
	      {
		rv[l].push_back(std::numeric_limits<double>::quiet_NaN());
	      }
	  }
      }
    else
      {
	fill_diploid_view_details(rv,diploids[ind].first->mutations,ind,twoN,remove_fixed,0);
	fill_diploid_view_details(rv,diploids[ind].first->smutations,ind,twoN,remove_fixed,0);
	fill_diploid_view_details(rv,diploids[ind].second->mutations,ind,twoN,remove_fixed,1);
	fill_diploid_view_details(rv,diploids[ind].second->smutations,ind,twoN,remove_fixed,1);
      }
    return rv;
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

  std::map<std::string, std::vector<double> > diploid_view_cpp(const singlepop_t *pop,
							       const size_t ind,
							       const int remove_fixed)
  {
    return diploid_view_details(pop->diploids,ind,remove_fixed);
  }

  std::map<std::string, std::vector<double> > diploid_view_cpp(const metapop_t * pop,
							       const size_t ind,
							       const int remove_fixed,
							       const int deme)
  {
    if (unsigned(deme) >= pop->diploids.size())
      {
	throw std::out_of_range("diploid_view: deme index out of range");
      }
    return diploid_view_details(pop->diploids[deme],ind,remove_fixed);
  }

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

