#ifndef __FWDPY_SERIALIZE_HPP__
#define __FWDPY_SERIALIZE_HPP__

#include <memory>
#include <vector>
#include <string>
#include <types.hpp>
#include <fwdpp/sugar/serialization.hpp>

namespace fwdpy
{
  namespace serialize
  {
    std::string serialize_singlepop(const singlepop_t * spop);
    void deserialize_singlepop(const std::vector<std::string> & strings,
			       std::vector<std::shared_ptr<singlepop_t> > * pops);
  }
}
  

#endif
