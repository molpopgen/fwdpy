#ifndef __FWDPY_INTERNAL_HPP__
#define __FWDPY_INTERNAL_HPP__

#include <vector>
#include <callbacks.hpp>
#include <fwdpp/extensions/callbacks.hpp>

namespace fwdpy
{
  namespace internal
  {
    struct region_manager
    {
      std::vector<KTfwd::extensions::shmodel> callbacks;
      std::vector<double> nb,ne,nw,sb,se,sw,rb,re,rw;
    };
  }
}

#endif
