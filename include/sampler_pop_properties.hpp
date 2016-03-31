#ifndef FWDPY_POP_PROPERTIES_HPP
#define FWDPY_POP_PROPERTIES_HPP

#include <vector>
#include <array>
#include <string>
#include "types.hpp"

/*
  NAMESPACE POLLUTION!!!!!

  Types defined here are in the global namespace.

  This is bad, BUT, it allows for auto-conversion of
  struct to dict via Cython.

  Currently, Cython will fail to compile auto-conversion
  code for structs declared inside a C++ namespace.
*/
struct qtrait_stats_cython
{
  std::string stat;
  double value;
  unsigned generation;
  qtrait_stats_cython(std::string _stat,
		      double _v, unsigned _g) : stat(std::move(_stat)),value(_v),generation(_g)
  {
  }
};

namespace fwdpy {
  /*!
    We use a vector of std::array here for compactness and speed
    We could use a std::map, or a vector of final_t (see below),
    but that gets slow and RAM-intensive.
  */
  using qtrait_stats_t = std::vector<std::array<double,10>>;

  class pop_properties
  /*!
    \brief A "sampler" that records "quantitative genetics" kinda stuff.
    \ingroup samplers
   */
  {
  public:
    using final_t = std::vector<qtrait_stats_cython>;

    inline void operator()(const singlepop_t * pop, gsl_rng *,
			   const unsigned generation)
    {
      std::vector<double> VG,VE,wbar,trait;
      VG.reserve(pop->diploids.size());
      VE.reserve(pop->diploids.size());
      trait.reserve(pop->diploids.size());
      wbar.reserve(pop->diploids.size());

      for(const auto & dip : pop->diploids)
	{
	  VG.push_back(dip.g);
	  VE.push_back(dip.e);
	  trait.push_back(dip.g+dip.e);
	  wbar.push_back(dip.w);
	}

      double twoN=2.*double(pop->diploids.size());
      double mvexpl = 0.,
        leading_e=std::numeric_limits<double>::quiet_NaN(),
        leading_f=std::numeric_limits<double>::quiet_NaN();
      double sum_e = 0.;
      unsigned nm=0;
      for(std::size_t i = 0 ; i < pop->mcounts.size() ; ++i )
        {
          if(pop->mcounts[i] && pop->mcounts[i]<twoN&&!pop->mutations[i].neutral)
            {
              auto n = pop->mcounts[i];
              double p=double(n)/twoN,q=1.-p;
	      double s = pop->mutations[i].s;
	      double temp = 2.*p*q*std::pow(s,2.0);
              if (temp > mvexpl)
                {
                  mvexpl=temp;
                  leading_e = s;
                  leading_f = p;
                }
	      sum_e += pop->mutations[i].s;
	      ++nm;
            }
        }

      //Calcate V(G) here b/c we're going to mess
      //around with this container below when
      //calculating V_{s,t}
      auto VG_ = gsl_stats_variance(VG.data(),1,VG.size());

      //Eq'n 5 from Zhang et al (2004) Genetics 166: 597,
      //but we calculate V_{G2} "manually"
      auto meanTrait = gsl_stats_mean(trait.data(),1,trait.size());
      auto meanG = gsl_stats_mean(VG.data(),1,VG.size());
      /*
	Transform arrays so that they represent squared deviations from mean.

	For "trait", this should be the optimum, but the API isn't allowing that right now...
      */
      std::transform(VG.begin(),VG.end(),VG.begin(),[meanG](const double d){ return std::pow(d-meanG,2.0);});
      std::transform(trait.begin(),trait.end(),trait.begin(),[this](const double d){ return std::pow(d-this->optimum,2.0);});
      //This is the apparent strenght of selection on the trait,
      //which is a regression of fitness onto trait value.
      double vst = -gsl_stats_variance(VG.data(),1,VG.size())/(2.0*gsl_stats_covariance(wbar.data(),1,trait.data(),1,wbar.size()));
      qstats.emplace_back(qtrait_stats_t::value_type{{double(generation),
	      VG_,
	      gsl_stats_variance(VE.data(),1,VE.size()),
	      leading_f,
	      leading_e,
	      mvexpl,
	      sum_e/double(nm),
	      gsl_stats_mean(wbar.data(),1,wbar.size()),
	      meanTrait,
	      vst}});
    }

    final_t final() const
    /*!
      The implementation details of this are messy,
      and the enum is used to guard against indexing errors.
    */
    {
      final_t rv;
      for( const auto & q : qstats)
	{
	  rv.emplace_back( "VG", q[static_cast<size_t>(qtrait_stat_list::VG)], unsigned(q[static_cast<size_t>(qtrait_stat_list::GEN)]) );
	  rv.emplace_back( "VE", q[static_cast<size_t>(qtrait_stat_list::VE)], unsigned(q[static_cast<size_t>(qtrait_stat_list::GEN)]) );
	  rv.emplace_back( "leading_q", q[static_cast<size_t>(qtrait_stat_list::PLF)], unsigned(q[static_cast<size_t>(qtrait_stat_list::GEN)]) );
	  rv.emplace_back( "leading_e", q[static_cast<size_t>(qtrait_stat_list::LE)], unsigned(q[static_cast<size_t>(qtrait_stat_list::GEN)]) );
	  rv.emplace_back( "max_expl", q[static_cast<size_t>(qtrait_stat_list::MAXEXP)], unsigned(q[static_cast<size_t>(qtrait_stat_list::GEN)]) );
	  rv.emplace_back( "ebar", q[static_cast<size_t>(qtrait_stat_list::EBAR)], unsigned(q[static_cast<size_t>(qtrait_stat_list::GEN)]) );
	  rv.emplace_back( "wbar", q[static_cast<size_t>(qtrait_stat_list::WBAR)], unsigned(q[static_cast<size_t>(qtrait_stat_list::GEN)]) );
	  rv.emplace_back( "tbar", q[static_cast<size_t>(qtrait_stat_list::TBAR)], unsigned(q[static_cast<size_t>(qtrait_stat_list::GEN)]) );
	  rv.emplace_back( "Vst", q[static_cast<size_t>(qtrait_stat_list::VST)], unsigned(q[static_cast<size_t>(qtrait_stat_list::GEN)]) );
	}
      return rv;
    }

    explicit pop_properties(double optimum_) : optimum(optimum_)
    {
    }
  private:
    qtrait_stats_t qstats;
    double optimum;
    enum class qtrait_stat_list : std::size_t { GEN,VG,VE,PLF,LE,MAXEXP,EBAR,WBAR,TBAR,VST };
  };
}
#endif
