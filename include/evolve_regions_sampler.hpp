#ifndef FWDPY_EVOLVE_REGIONS_SAMPLER_HPP
#define FWDPY_EVOLVE_REGIONS_SAMPLER_HPP
#include <vector>
#include <memory>
#include "types.hpp"
#include "sampler_base.hpp"
#include "internal_region_manager.hpp"
#include "fwdpy_fitness.hpp"
namespace fwdpy
{
  void evolve_regions_sampler_cpp( GSLrng_t * rng,
				   std::vector<std::shared_ptr<singlepop_t> > * pops,
				   std::vector<std::unique_ptr<sampler_base> > & samplers,
				   const unsigned * Nvector,
				   const size_t Nvector_length,
				   const double mu_neutral,
				   const double mu_selected,
				   const double littler,
				   const double f,
				   const int sample,
				   const internal::region_manager * rm,
				   const singlepop_fitness & fitness);
} //ns fwdpy
#endif
