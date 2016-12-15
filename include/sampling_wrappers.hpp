/*!
  \file sampling_wrappers.hpp
  \brief Wrappers simplifying access to fwdpp's "sugar" fxns for taking samples
  from populations

  fwdpp's templates in fwdpp/sugar/sampling.hpp have an 'auto' return type,
  which Cython can't
  understand.  Hence, these thin wrappers "force" fwdpp to return the right
  thing.
*/

#ifndef FWDPY_SAMPLING_WRAPPERS_HPP
#define FWDPY_SAMPLING_WRAPPERS_HPP

#include "types.hpp"
#include <fwdpp/sugar/sampling.hpp>

namespace fwdpy
{
    template <typename poptype>
    KTfwd::sample_t
    sample_single(const gsl_rng *r, const poptype &p, const unsigned nsam,
                  const bool removeFixed)
    //! Single-deme, single region
    {
        return KTfwd::sample(r, p, nsam, removeFixed);
    }

    template <typename poptype>
    KTfwd::sep_sample_t
    sample_sep_single(const gsl_rng *r, const poptype &p, const unsigned nsam,
                      const bool removeFixed)
    //! Single deme, single region
    {
        return KTfwd::sample_separate(r, p, nsam, removeFixed);
    }

    template <typename poptype>
    std::vector<KTfwd::sample_t>
    sample_single_mloc(
        const gsl_rng *r, const poptype &p, const unsigned nsam,
        const bool removeFixed,
        const std::vector<std::pair<double, double>> &locus_boundaries)
    //! Single deme, multi-region
    {
        return KTfwd::sample(r, p, nsam, removeFixed, locus_boundaries);
    }

    template <typename poptype>
    std::vector<KTfwd::sep_sample_t>
    sample_sep_single_mloc(
        const gsl_rng *r, const poptype &p, const unsigned nsam,
        const bool removeFixed,
        const std::vector<std::pair<double, double>> &locus_boundaries)

    //! Single deme, multi-region
    {
        return KTfwd::sample_separate(r, p, nsam, removeFixed,
                                      locus_boundaries);
    }
}

#endif
