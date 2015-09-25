#ifndef __FWDPY_EWVW_HPP__
#define __FWDPY_EWVW_HPP__

#include <vector>
#include <memory>
#include <fwdpp/extensions/regions.hpp>
#include <types.hpp>

namespace fwdpy
{
  namespace qtrait {
    void evolve_ewvw_t( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
			const unsigned * Nvector,
			const size_t Nvector_length,
			const double mu_neutral,
			const double mu_selected,
			const double littler,
			const double f,
			const double prop_vw,
			const double optimum,
			const int track,
			const std::vector<double> & nbegs,
			const std::vector<double> & nends,
			const std::vector<double> & nweights,
			const std::vector<double> & sbegs,
			const std::vector<double> & sends,
			const std::vector<double> & sweights,
			const std::vector<KTfwd::extensions::shmodel> * callbacks,
			const std::vector<double> & rbeg,
			const std::vector<double> & rend,
			const std::vector<double> & rweight);
  }
}

#endif
