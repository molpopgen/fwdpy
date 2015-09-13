#include <types.hpp>

namespace fwdpy {
void singlepop_t::updateTraj() {
  for( const auto & __m : mutations )
    {
      if( !__m.neutral )
        {
          auto __p = std::make_pair(__m.g,std::make_pair(__m.pos,__m.s));
          auto __itr = trajectories.find(__p);
          if(__itr == trajectories.end())
            {
              trajectories[__p] = std::vector<double>(1,double(__m.n)/double(2*diploids.size()));
            }
          else
            {
              //Don't keep updating for fixed variants
              if( *(__itr->second.end()-1) < 1.) //2*diploids.size() )
                {
                  __itr->second.push_back(double(__m.n)/double(2*diploids.size()));
                }
            }
        }
    }
}
}
