#ifndef FWDPY_EVOLVE_REGIONS_SAMPLER_HPP
#define FWDPY_EVOLVE_REGIONS_SAMPLER_HPP

#include <vector>
#include <type_traits>
#include <types.hpp>
#include <internal/internal.hpp>

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
			     const int track,
			     const fwdpy::internal::region_manager * rm,
			     const char * fitness,
			     const sampler & s,
			     const int interval) -> std::vector< std::pair<unsigned,typename std::result_of<sampler(const singlepop_t *, gsl_rng *)>::type> >
{
  std::vector< std::pair<unsigned,typename std::result_of<sampler(const singlepop_t *, gsl_rng *)>::type> > rv;

  return rv;
}

#endif
