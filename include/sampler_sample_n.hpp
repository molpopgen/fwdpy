#ifndef FWDPY_SAMPLE_N_HPP
#define FWDPY_SAMPLE_N_HPP

#include <vector>
#include <utility>
#include <algorithm>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/sampling.hpp>
#include "types.hpp"

//Namespace pollution!!
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

namespace fwdpy
{
  class sample_n //take a sample of size n from a population
  /*
    \brief A "sampler" that takes a sample of n gametes from a population
    \ingroup samplers
   */
  {
  public:
    using final_t = std::vector<std::pair<unsigned,detailed_deme_sample> >;
    inline void operator()(const singlepop_t * pop,gsl_rng * r,
			   const unsigned generation)
    {
      auto s = KTfwd::sample_separate(r,*pop,nsam,true);
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
    explicit sample_n(unsigned nsam_) : rv(final_t()),nsam(nsam_)
    {
    }
  private:
    final_t rv;
    const unsigned nsam;
  };
}

#endif
