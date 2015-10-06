#ifndef __FWDPY_TRAJECTORIES_HPP__
#define __FWDPY_TRAJECTORIES_HPP__

#include <map>
#include <string>
#include <vector>

#include <types.hpp>
namespace fwdpy {
  std::map<std::string,std::vector<double> > get_singlepop_traj(const singlepop_t *pop,const unsigned minsojourn,const double minfreq);
}

#endif
