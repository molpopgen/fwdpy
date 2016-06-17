#ifndef FWDPY_FITNESS_MODELS_HPP
#define FWDPY_FITNESS_MODELS_HPP

#include "types.hpp"
#include <functional>
namespace fwdpy
{
  struct singlepop_fitness
  /*!
    Base class for fitness schemes for single-deme simulations
   */
  {
    //! This is the form of a fitness function for a single-deme simulation
    using fitness_fxn_t = std::function<double(const diploid_t &,const gcont_t &, const mcont_t &)>;

    //! The fitness function itself
    fitness_fxn_t fitness_function;

    /*!
      Placeholder for future functionality
     */
    virtual void update(const singlepop_t *) {}

    //! Allows us to allocate on stack in Cython
    singlepop_fitness() : fitness_function(fitness_fxn_t()) {}
    //! Constructor is a sink for a fitness_fxn_t 
    singlepop_fitness(fitness_fxn_t ff) : fitness_function(std::move(ff)) {}
  };
    
  inline singlepop_fitness make_additive_fitness(double scaling = 2.0)
  {
    return singlepop_fitness(std::bind(KTfwd::additive_diploid(),
				       std::placeholders::_1,
				       std::placeholders::_2,
				       std::placeholders::_3,
				       scaling));
  }
  inline singlepop_fitness make_multiplicative_fitness(double scaling = 2.0)
  {
    return singlepop_fitness(std::bind(KTfwd::multiplicative_diploid(),
				       std::placeholders::_1,
				       std::placeholders::_2,
				       std::placeholders::_3,
				       scaling));
  }
}


#endif
