#ifndef __FWDPY_SERIALIZE_HPP__
#define __FWDPY_SERIALIZE_HPP__

#include <memory>
#include <vector>
#include <string>
#include <fwdpp/sugar/serialization.hpp>
#include "types.hpp"

namespace fwdpy
{
  namespace serialize
  {
    std::string serialize_singlepop(const singlepop_t * spop);
    std::vector<std::shared_ptr<singlepop_t> > deserialize_singlepop(const std::vector<std::string> & strings);
    std::string serialize_metapop(const fwdpy::metapop_t * pop);
    std::vector<std::shared_ptr<metapop_t> > deserialize_metapop(const std::vector<std::string> & strings);
  }
}
  

#endif
