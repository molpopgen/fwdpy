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
      region_manager() : callbacks(std::vector<KTfwd::extensions::shmodel>()),
			 nb(std::vector<double>()),
			 ne(std::vector<double>()),
			 nw(std::vector<double>()),
			 sb(std::vector<double>()),
			 se(std::vector<double>()),
			 sw(std::vector<double>()),
			 rb(std::vector<double>()),
			 re(std::vector<double>()),
			 rw(std::vector<double>())
      {
      }
    };
  }
}

#endif
