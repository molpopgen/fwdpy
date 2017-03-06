#ifndef __FWDPY_HOCRULES_HPP__
#define __FWDPY_HOCRULES_HPP__

#include "rules_base.hpp"
#include <cmath>
#include <gsl/gsl_sf_pow_int.h>
/*
  Custom "rules" policy for single-region House-of-Cards simulations.

  This file is partly KT having fun, but also an effort to stress-test fwdpp's
  "experimental" API.

  The advantage of this struct is:
  1. An offspring has its G and E values automatically assigned.  This allows
  us to record the
  exact fitnesses/heritabilities used in the simulation over time.
*/

namespace fwdpy
{
    namespace qtrait
    {
        struct qtrait_model_rules : public fwdpy::single_region_rules_base
        {
            using base_t = fwdpy::single_region_rules_base;
            const double sigE, optimum, VS;
            qtrait_model_rules(const double &sigE_, const double &optimum_,
                               const double &VS_,
                               const unsigned maxN_ = 100000,
                               const int power_ = 2) noexcept(false)
                : base_t(), sigE(sigE_), optimum(optimum_), VS(VS_)
            /*!
              Constructor throws std::runtime_error if params are not valid.
            */
            {
                if (sigE < 0.)
                    throw std::runtime_error(
                        "Environmental noise term must be >= 0.");
                if (VS <= 0.)
                    throw std::runtime_error("VS must be > 0.");
            }

            qtrait_model_rules(qtrait_model_rules &&) = default;

            qtrait_model_rules(const qtrait_model_rules &rhs)
                : base_t(rhs), sigE(rhs.sigE), optimum(rhs.optimum), VS(rhs.VS)
            {
            }

            virtual void
            w(const dipvector_t &diploids, gcont_t &gametes,
              const mcont_t &mutations)
            {
                auto N_curr = diploids.size();
                if (fitnesses.size() < N_curr)
                    fitnesses.resize(N_curr);
                wbar = 0.;
                for (size_t i = 0; i < N_curr; ++i)
                    {
                        gametes[diploids[i].first].n
                            = gametes[diploids[i].second].n = 0;
                        fitnesses[i] = diploids[i].w;
                        wbar += diploids[i].w;
                    }
                wbar /= double(N_curr);
                lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(
                    gsl_ran_discrete_preproc(N_curr, &fitnesses[0]));
            }

            //! \brief Update some property of the offspring based on
            //! properties of the parents
            virtual void
            update(const gsl_rng *r, diploid_t &offspring, const diploid_t &,
                   const diploid_t &, const gcont_t &gametes,
                   const mcont_t &mutations,
                   const single_region_fitness_fxn &ff) noexcept
            {
                offspring.g = ff(offspring, gametes, mutations);
                offspring.e = gsl_ran_gaussian_ziggurat(r, sigE);
                double dev = (offspring.g + offspring.e - optimum);
                offspring.w = std::exp(-(dev * dev) / (2. * VS));
                assert(std::isfinite(offspring.w));
                return;
            }
        };
    } // namespace qtrait
} // namespace fwdpy
#endif
