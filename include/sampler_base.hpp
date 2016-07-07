#ifndef FWDPY_SAMPLER_BASE_HPP
#define FWDPY_SAMPLER_BASE_HPP

#include "types.hpp"
#include <stdexcept>
#include <thread>
#include <vector>

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
    virtual void cleanup()
    /*!
      Evolve functions will call this before returning.
      This function gives samplers a chance to do things like clean
      up any RAM that they may have allocated during the course of 
      a simulation. 

      \note: the stuff to return to Python should not be cleaned up by this function!
     */
    {
    }
    virtual ~sampler_base(){}
  };

  template<typename pop_t>
  inline void apply_sampler_wrapper( sampler_base * s, const pop_t * pop )
  /*!
    Thin wrapper function for apply_sampler_cpp.  This funcion is needed to 
    avoid slicing the pointer to a sampler down to a pointer to the base class.
   */
  {
    s->operator()(pop,pop->generation);
  }
  
  template<typename T>
  inline void apply_sampler_cpp( const std::vector<std::shared_ptr<T> > & popvec,
				 const std::vector<std::unique_ptr<sampler_base> > & samplers )
  /*!
    Apply the i-th ampler to the i-th pop using std::thread.

    Throws runtime_error if popvec.size()!=samplers.size().
   */
  {
    if(popvec.size()!=samplers.size()) throw std::runtime_error("Containers of populations and samplers must be equal in length");
    std::vector<std::thread> threads;
    for(std::size_t i=0;i<popvec.size();++i)
      {
	threads.emplace_back(apply_sampler_wrapper<T>,samplers[i].get(),popvec[i].get());
      }
    for(auto & t : threads) t.join();
  }

  template<typename final_t>
  struct custom_sampler : public sampler_base
  {
    using singlepop_call_operator = void(*)(const singlepop_t *, const unsigned, final_t & );
    using multilocus_call_operator = void(*)(const multilocus_t *, const unsigned, final_t & );
    using metapop_call_operator = void(*)(const metapop_t *, const unsigned, final_t & );    
    final_t f;
    singlepop_call_operator scall;
    multilocus_call_operator mcall;
    metapop_call_operator metacall;

    custom_sampler( singlepop_call_operator s ) : scall(s),mcall(nullptr),metacall(nullptr)
    {
    }

    void operator()(const singlepop_t * pop ,const unsigned generation)
    {
      if(scall != nullptr)
	{
	  scall(pop,generation,f);
	}
      else sampler_base::operator()(pop,generation);
    };
    void operator()(const multilocus_t * pop ,const unsigned generation)
    {
      if(mcall != nullptr)
	{
	  mcall(pop,generation,f);
	}
      else sampler_base::operator()(pop,generation);
    };
    void operator()(const metapop_t * pop,const unsigned generation)
    {
      if(metacall != nullptr)
	{
	  metacall(pop,generation,f);
	}
      else sampler_base::operator()(pop,generation);
    };    

    final_t final()
    {
      return f;
    }
  };
}

#endif
