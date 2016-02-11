#ifndef FWDPY_EVOLVE_REGIONS_SAMPLER_HPP
#define FWDPY_EVOLVE_REGIONS_SAMPLER_HPP

#include <vector>
#include <type_traits>
#include <types.hpp>
#include <internal/internal.hpp>

/*
  This is the generic function.

  The type "sampler" must define an operator() taking a const singlepop_t * and a gsl_rng * as arguments.
  
  Other arguments may be bound by the caller, etc.

  The object s will be applied every "interval" generations.

  If result_t is the return type of s(const singlepop_t *, gsl_rng *), then the return type
  of this function is vector<pair<unsigned,result_t> >, where the unsigned values refer
  to the generation in which a specific result_t was obtained.

  Design issues to work out:

  1. The result_t is an object is considered below.  What happens when it is 
  a vector of objects?

  For example, a sample is a single object.  The selected mutation frequencies are a vector
  of data on mutations.

  Perhaps that is just ok? It isn't super-friendly, and would need conversion to pandas.DataFrame,
  and I guess we'd have to supply those functions.

  The problem is one of space.  Doing frequency trajectories this way will be much more RAM-intensive
  Than the current implementation.  The "qtrait stats" is probably fine, with some overhead from strings.
*/
template<typename sampler>
inline auto
evolve_regions_sampler_async(GSLrng_t * rng,
			     std::vector<std::shared_ptr<singlepop_t> > * pops,
			     const unsigned * Nvector,
			     const size_t Nvector_len,
			     const double mu_neutral,
			     const double mu_selected,
			     const double littler,
			     const double f,
			     const fwdpy::internal::region_manager * rm,
			     const char * fitness,
			     const sampler & s,
			     const int interval) -> std::vector< std::pair<unsigned,typename std::result_of<sampler(const singlepop_t *, gsl_rng *)>::type> >
{
  std::vector< std::pair<unsigned,typename std::result_of<sampler(const singlepop_t *, gsl_rng *)>::type> > rv;

  return rv;
}

#endif
