/*!
  \file demography_migrates.hpp

  \brief Defines object for choosing parents in models with migration.
*/
#ifndef FWDPY_DEMOGRAPHY_MIGRATES_HPP
#define FWDPY_DEMOGRAPHY_MIGRATES_HPP

#include <fwdpp/internal/gsl_discrete.hpp>
#include <stdexcept>
#include <vector>

namespace fwdpy
{
    namespace demography
    {
        struct migrates
        /*!
          Structure representing weights on migration probabilities.
         */
        {
            using lookup_t = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr;
            using vec_lookup_t = std::vector<lookup_t>;
            vec_lookup_t lookups;
            migrates() : lookups(vec_lookup_t())
            /*!
              Satisfies Cython's requirements for stack allocation
            */
            {
            }
            migrates(const std::vector<std::vector<double>> &weights)
                : lookups(vec_lookup_t())
            /*!
              Constructor

              All elements in weights must be equal-length.

              \note In Cython/Python, it must also be enforced that
              weights.size() == number of demes
              in the simulation.
             */
            {
                lookups.reserve(weights.size());
                if (weights.empty())
                    return;
                const auto s = weights[0].size();
                for (auto &i : weights)
                    {
                        if (i.size() != s)
                            throw std::runtime_error("all vectors of "
                                                     "migration weights must "
                                                     "be same length");
                        lookups.emplace_back(
                            lookup_t(gsl_ran_discrete_preproc(s, i.data())));
                    }
            }

            inline std::size_t
            operator()(const size_t deme, const gsl_rng *r) const
            /*!
              Call operator conforms to fwdpp's requirements
            */
            {
                return gsl_ran_discrete(r, lookups[deme].get());
            }
        };

        migrates
        make_migrates(const std::vector<std::vector<double>> &weights)
        /*!
          Convenience function makes life easier in Cython via
          cdef migrates x = make_migrates(y)
        */
        {
            return migrates(weights);
        }
    }
}

#endif
