#ifndef FWDPY_GET_SELECTED_MUT_DATA_HPP
#define FWDPY_GET_SELECTED_MUT_DATA_HPP

#include <limits>
#include "types.hpp"

struct selected_mut_data
{
  double pos,esize;
  unsigned origin;
  selected_mut_data(unsigned g,double p,double e) :
    pos(p),esize(e),origin(g)
  {
  }
  selected_mut_data() : pos(std::numeric_limits<double>::quiet_NaN()),
			esize(std::numeric_limits<double>::quiet_NaN()),
			origin(std::numeric_limits<unsigned>::max())
  {
  }
  inline bool operator==(const selected_mut_data & rhs) const noexcept
  {
    return this->origin == rhs.origin &&
      this->pos == rhs.pos &&
      this->esize == rhs.esize;
  }
};

namespace fwdpy
{
  //! Used internally to convert C++11 types to something Cython will understand
  enum class traj_key_values : std::size_t { deme,origin,pos,esize };

  /*!
    \brief Unique key for a mutation.  Used when tracking mutation frequencies.

    Values are: deme, generation of mutation origin, position, effect size.

    \note Used in fwdpy::selected_mut_tracker
  */
  using trajectories_key_t = std::tuple<unsigned,unsigned,double,double>;
  /*!
    \brief Internal representation of mutation frequencies during a simulation

    \note Used in fwdpy::selected_mut_tracker
  */
  using trajectories_t = std::map< trajectories_key_t , std::vector<double> >;
  
  class selected_mut_tracker
  /*!
    \brief A "sampler" for recording frequency trajectories of selected mutations.
    \ingroup samplers
  */
  {
  public:
    using final_t = std::vector< std::pair<selected_mut_data, std::vector<double> > >;
    inline void operator()(const singlepop_t * pop,gsl_rng * ,
			   const unsigned)
    {
      for(std::size_t i = 0 ; i < pop->mcounts.size() ; ++i )
      	{
	  if(pop->mcounts[i]) //if mutation is not extinct
	    {
	      const auto & __m = pop->mutations[i];
	      if( !__m.neutral )
		{
		  const auto freq = double(pop->mcounts[i])/double(2*pop->diploids.size());
		  auto __p = std::make_tuple(0u,__m.g,__m.pos,__m.s);
		  auto __itr = trajectories.find(__p);
		  if(__itr == trajectories.end())
		    {
		      trajectories[__p] = std::vector<double>(1,freq);
		    }
		  else
		    {
		      //Don't keep updating for fixed variants
		      if( __itr->second.back() < 1.)
			{
			  __itr->second.push_back(freq);
			}
		    }
		}
	    }
      	}
    }
    
    final_t final() const
    {
      final_t rv;
      for( const auto & i : trajectories )
	{
	  rv.emplace_back( std::make_pair( selected_mut_data(std::get<static_cast<std::size_t>(traj_key_values::origin)>(i.first),
							     std::get<static_cast<std::size_t>(traj_key_values::pos)>(i.first),
							     std::get<static_cast<std::size_t>(traj_key_values::esize)>(i.first)),
					   i.second) );
	}
      return rv;
    }
    
    explicit selected_mut_tracker() : trajectories(trajectories_t())
    {
    }
  private:
    trajectories_t trajectories;
  };
}

#endif
