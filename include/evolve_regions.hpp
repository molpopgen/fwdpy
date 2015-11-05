#ifndef __FWDPY_EVOLVE_REGIONS_HPP__
#define __FWDPY_EVOLVE_REGIONS_HPP__

#include <types.hpp>
#include <vector>
#include <internal/internal.hpp>

namespace fwdpy
{
  void evolve_regions_t( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
			 const unsigned * Nvector,
			 const size_t Nvector_length,
			 const double mu_neutral,
			 const double mu_selected,
			 const double littler,
			 const double f,
			 const int track,
			 const internal::region_manager * rm,
			 const char * fitness);

  void split_and_evolve_t(GSLrng_t * rng,
			  std::vector<std::shared_ptr<metapop_t> > * mpops,
			  const unsigned * Nvector_A,
			  const size_t Nvector_A_len,
			  const unsigned * Nvector_B,
			  const size_t Nvector_B_len,
			  const double neutral,
			  const double selected,
			  const double recrate,
			  const std::vector<double> & fs,
			  const internal::region_manager * rm,
			  const char * fitness);
}

#endif
