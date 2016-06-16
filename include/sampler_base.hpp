#ifndef FWDPY_SAMPLER_BASE_HPP
#define FWDPY_SAMPLER_BASE_HPP

#include "types.hpp"
#include <stdexcept>
namespace fwdpy
{
  struct sampler_base
  {
    virtual void operator()(const singlepop_t *,const unsigned)
    {
      throw std::runtime_error("sampler type not implemented for single deme simulations");
    };
    virtual void operator()(const multilocus_t *,const unsigned)
    {
      throw std::runtime_error("sampler type not implemented for multi-locus simulations");
    };
    virtual void operator()(const metapop_t *,const unsigned)
    {
      throw std::runtime_error("sampler type not implemented for metapopulation simulations");
    };
    virtual ~sampler_base(){}
  };
}

#endif
