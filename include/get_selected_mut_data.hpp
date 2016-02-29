#ifndef FWDPY_GET_SELECTED_MUT_DATA_HPP
#define FWDPY_GET_SELECTED_MUT_DATA_HPP

#include <types.hpp>

struct selected_mut_data
{
  double pos,esize;
  unsigned origin;
  selected_mut_data(unsigned g,double p,double f,double e) :
    origin(g),pos(p),esize(e)
  {
  }
};

namespace fwdpy
{
  class get_selected_mut_data //record info on selected mutations in population, including fixations
  {
  public:
    using final_t = std::map<std::string,std::vector<double> >;
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
      std::vector<double> pos,freq,s;
      std::vector<double> generations;
      for( auto itr = this->trajectories.cbegin() ;
	   itr != this->trajectories.cend() ; ++itr )
	{
	  std::vector<unsigned> times(itr->second.size());
	  unsigned itime = std::get<static_cast<std::size_t>(traj_key_values::origin)>(itr->first);
	  generate(times.begin(),times.end(),[&itime]{ return itime++; });
	  generations.insert(generations.end(),times.begin(),times.end());
	  fill_n(std::back_inserter(pos),itr->second.size(),std::get<static_cast<std::size_t>(traj_key_values::pos)>(itr->first));
	  fill_n(std::back_inserter(s),itr->second.size(),std::get<static_cast<std::size_t>(traj_key_values::esize)>(itr->first));
	  std::copy(itr->second.begin(),itr->second.end(),back_inserter(freq));
	}
      return final_t{
	{"pos",std::move(pos)},
	  {"freq",std::move(freq)},
	    {"generation",std::move(generations)},
	      {"esize",std::move(s)}
      };
    }
    
    explicit get_selected_mut_data() : trajectories(trajectories_t())
    {
    }
  private:
    trajectories_t trajectories;
  };
}

#endif
