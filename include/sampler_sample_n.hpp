#ifndef FWDPY_SAMPLE_N_HPP
#define FWDPY_SAMPLE_N_HPP

#include <vector>
#include <utility>
#include <algorithm>
#include <type_traits>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/sampling.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>
#include "types.hpp"

namespace fwdpy
{
  struct detailed_deme_sample
  {
    KTfwd::sep_sample_t genotypes;
    std::vector<std::pair<double,double> > sh;
    template<typename T1,typename T2>
    detailed_deme_sample(T1 && t1, T2 && t2) : genotypes(std::forward<T1>(t1)),
					       sh(std::forward<T2>(t2))
    {
    }
  };


  template<typename pop_t = singlepop_t>
  class sample_n //take a sample of size n from a population
  /*
    \brief A "sampler" that takes a sample of n gametes from a population
    \ingroup samplers
  */
  {
  public:
    using singlepop_type_return_t = std::vector<std::pair<unsigned,detailed_deme_sample> >;
    using multilocus_type_return_t = std::vector<std::pair<unsigned,std::vector<detailed_deme_sample> > >;
    using final_t = typename std::conditional < std::is_same<typename pop_t::popmodel_t,KTfwd::sugar::SINGLEPOP_TAG>::value,
						singlepop_type_return_t,
						multilocus_type_return_t >::type;
    inline void operator()(const pop_t * pop,
			   const unsigned generation)
    {
      auto s = KTfwd::sample_separate(r.get(),*pop,nsam,true);
      std::vector< std::pair<double,double> > sh;
      for( const auto & i : s.second )
	{
	  auto itr = std::find_if(pop->mutations.begin(),pop->mutations.end(),[&i](const singlepop_t::mutation_t & m) noexcept
				  {
				    return m.pos == i.first;
				  });
	  sh.emplace_back(itr->s,itr->h);
	}
      rv.emplace_back(generation,detailed_deme_sample(std::move(s),std::move(sh)));
    }

    final_t final() const
    {
      return rv;
    }
    explicit sample_n(unsigned nsam_, const gsl_rng * r_) :
      rv(final_t()),nsam(nsam_),r(GSLrng_t(gsl_rng_get(r_)))
      /*!
	Note the implementation of this constructor!!

	By taking a gsl_rng * from outside, we are able to guarantee
	that this object is reproducibly seeded to the extent that
	this constructor is called in a reproducible order.
      */
    {
    }
  private:
    final_t rv;
    const unsigned nsam;
    GSLrng_t r;
  };

  template<>
  inline void sample_n<multilocus_t>::operator()(const multilocus_t * pop,
						 const unsigned generation)
  {
    auto s = KTfwd::sample_separate(r.get(),*pop,nsam,true);
    std::vector<detailed_deme_sample> vds;
    for(unsigned i=0;i<s.size();++i)
      {	
	std::vector< std::pair<double,double> > sh;
	for( const auto & si : s[i].second)
	  {
	    auto itr = std::find_if(pop->mutations.begin(),pop->mutations.end(),[&si](const singlepop_t::mutation_t & m) noexcept
				    {
				      return m.pos == si.first;
				    });
	    sh.emplace_back(itr->s,itr->h);
	  }
	vds.emplace_back(std::move(s[i]),std::move(sh));
      }
    rv.emplace_back(generation,std::move(vds));
  }
}

#endif
