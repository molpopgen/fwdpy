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

  template<typename T>
  inline void apply_sampler_cpp( const std::vector<std::shared_ptr<T> > & popvec,
				 const std::vector<std::unique_ptr<sampler_base> > & samplers )
  {
    if(popvec.size()!=samplers.size()) throw std::runtime_error("Containers of populations and samplers must be equal in length");
    for(std::size_t i=0;i<popvec.size();++i)
      {
	samplers[i]->operator()(popvec[i].get(),popvec[i]->generation);
      }
  }
}

#endif
