#ifndef FWDPY_QTRAIT_EVOLVE_MLOCUS_HPP
#define FWDPY_QTRAIT_EVOLVE_MLOCUS_HPP

#include "fwdpp_features.hpp"
#include "fwdpy_fitness.hpp"
#include "internal_region_manager.hpp"
#include "reserve.hpp"
#include "sampler_base.hpp"
#include "sampler_base.hpp"
#include "types.hpp"
#include <algorithm>
#include <exception>
#include <functional>
#include <future>
#include <fwdpp/diploid.hh>
#include <fwdpp/experimental/sample_diploid_mloc.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <set>
#include <type_traits>
#include <vector>

namespace fwdpy
{
    namespace qtrait
    {
        template <typename mutation_policies, typename recombination_policies,
                  typename rules_type>
        inline void
        evolve_qtrait_mloc_details_common(
            multilocus_t *pop, gsl_rng const *rng,
            const KTfwd::uint_t *Nvector, const std::size_t Nvector_len,
            const mutation_policies &mmodels,
            const recombination_policies &recpols,
            const std::vector<double> &tmu,
            const std::vector<double> &between_region_rec_rates,
            std::unique_ptr<multilocus_fitness> &fitness, sampler_base &s,
            const unsigned interval, const double f, rules_type &&rules)
        /*!
         * Common loop shared by the two functions defined
         * below
         */
        {
            auto rules_local(std::forward<rules_type>(rules));
            // evolve...
            const unsigned simlen = unsigned(Nvector_len);
            // fitness->update(pop);
            for (unsigned g = 0; g < simlen; ++g, ++pop->generation)
                {
                    const unsigned nextN = *(Nvector + g);
                    if (interval && pop->generation
                        && pop->generation % interval == 0.)
                        {
                            s(pop, pop->generation);
                        }
                    KTfwd::experimental::sample_diploid(
                        rng, pop->gametes, pop->diploids, pop->mutations,
                        pop->mcounts, pop->N, nextN, tmu.data(), mmodels,
                        recpols, between_region_rec_rates.data(),
                        // rec b/w loci is interpreted as cM!!!!!
                        [](const gsl_rng *__r, const double __d) {
                            return gsl_ran_bernoulli(__r, __d);
                        },
                        fitness->fitness_function, pop->neutral, pop->selected,
                        f, rules_local, KTfwd::remove_neutral());
                    fwdpy::update_mutations_n(pop->mutations, pop->fixations,
                                              pop->fixation_times,
                                              pop->mut_lookup, pop->mcounts,
                                              pop->generation, 2 * nextN);
                    pop->N = nextN;
                    // fitness->update(pop);
                }
            if (interval && pop->generation
                && pop->generation % interval == 0.)
                {
                    s(pop, pop->generation);
                }
        }

        template <typename rules_type>
        inline void
        evolve_qtrait_mloc_regions_cpp_details(
            fwdpy::multilocus_t *pop,
            std::unique_ptr<multilocus_fitness> &fitness, sampler_base &s,
            const unsigned long seed, const unsigned *Nvector,
            const size_t Nvector_len, const internal::region_manager *rm,
            const std::vector<double> &between_region_rec_rates,
            const double f, const int interval, rules_type &&rules)
        /*!
         * Evolve a multilocus model with support for "regions".
         * Current region support is limited: 1 neutral, 1 selected,
         * and 1 rec region per "locus".
         */
        {
            // We need to set up mutation and recombination
            // models
            // Get local rng 4 this thread
            gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(rng, seed);
            std::vector<double> tmu;
            std::vector<std::function<std::vector<double>(
                const multilocus_t::gamete_t &, const multilocus_t::gamete_t &,
                const multilocus_t::mcont_t &)>>
                recpols;
            std::vector<std::function<std::size_t(std::queue<std::size_t> &,
                                                  multilocus_t::mcont_t &)>>
                mmodels;

            for (std::size_t i = 0; i < rm->nb.size(); ++i)
                {
                    tmu.push_back(rm->nw[i] + rm->sw[i]);
                    mmodels.emplace_back(std::bind(
                        KTfwd::infsites(), std::placeholders::_1,
                        std::placeholders::_2, rng, std::ref(pop->mut_lookup),
                        &pop->generation,
                        rm->nw[i], // mutation rate, neutral
                        rm->sw[i], // mutation rate, selected
                        [&rng, &rm, i]() {
                            return gsl_ran_flat(rng, rm->nb[i], rm->ne[i]);
                        }, // neutral mutation pos'n
                        [&rng, &rm, i]() {
                            return gsl_ran_flat(rng, rm->sb[i], rm->se[i]);
                        }, // selected mutation pos'n
                        [&rm, i, rng]() { return rm->callbacks[i].s(rng); },
                        [&rm, i, rng]() { return rm->callbacks[i].h(rng); }));
                    recpols.emplace_back(std::bind(
                        KTfwd::poisson_xover(), rng, rm->rw[i], rm->rb[i],
                        rm->re[i], std::placeholders::_1,
                        std::placeholders::_2, std::placeholders::_3));
                }
            auto x = std::max_element(Nvector, Nvector + Nvector_len);
            reserve_space(pop->gametes, pop->mutations, *x,
                          std::accumulate(tmu.begin(), tmu.end(), 0.));
            evolve_qtrait_mloc_details_common(
                pop, rng, Nvector, Nvector_len, mmodels, recpols, tmu,
                between_region_rec_rates, fitness, s, interval, f, rules);
            s.cleanup();
            gsl_rng_free(rng);
        }

        template <typename rules_type>
        inline void
        evolve_qtrait_mloc_cpp_details(
            fwdpy::multilocus_t *pop,
            std::unique_ptr<multilocus_fitness> &fitness, sampler_base &s,
            const unsigned long seed, const unsigned *Nvector,
            const size_t Nvector_len,
            const std::vector<double> &neutral_mutation_rates,
            const std::vector<double> &selected_mutation_rates,
            const std::vector<KTfwd::extensions::shmodel> &effects_dominance,
            const std::vector<double> &within_region_rec_rates,
            const std::vector<double> &between_region_rec_rates,
            const double f, const int interval, rules_type &&rules)
        /*!
         * \deprecated
         * Simplistic evolution of multi-locus quant-trait model.
         */
        {
            // TODO: check that vector input sizes ok.  Ideally, this needs to
            // be done b4 this point, else we
            // have threads throwing exceptions...

            // Reserve space
            auto x = std::max_element(Nvector, Nvector + Nvector_len);
            assert(x != Nvector + Nvector_len);
            reserve_space(
                pop->gametes, pop->mutations, *x,
                std::accumulate(neutral_mutation_rates.begin(),
                                neutral_mutation_rates.end(), 0.)
                    + std::accumulate(selected_mutation_rates.begin(),
                                      selected_mutation_rates.end(), 0.));

            // Get local rng 4 this thread
            gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(rng, seed);

            // Establish recombination maps--uniform w/in each locus
            std::vector<std::function<std::vector<double>(
                const multilocus_t::gamete_t &, const multilocus_t::gamete_t &,
                const multilocus_t::mcont_t &)>>
                recpols;
            unsigned i = 0;
            for (auto ri : within_region_rec_rates)
                {
                    recpols.emplace_back(std::bind(
                        KTfwd::poisson_xover(), rng, ri, double(i),
                        double(i + 1), std::placeholders::_1,
                        std::placeholders::_2, std::placeholders::_3));
                    ++i;
                }

            // Establish mutation models--uniform process w/in each locus,
            // w/mutations affecting trait value having DFE N(0,sigma_mus[i])
            // at the i-th locus
            i = 0;
            std::vector<std::function<std::size_t(std::queue<std::size_t> &,
                                                  multilocus_t::mcont_t &)>>
                mmodels;
            for (auto mi : neutral_mutation_rates)
                {
                    mmodels.emplace_back(std::bind(
                        KTfwd::infsites(), std::placeholders::_1,
                        std::placeholders::_2, rng, std::ref(pop->mut_lookup),
                        &pop->generation,
                        mi,                         // mutation rate
                        selected_mutation_rates[i], // mutation rate
                        [&rng, i]() {
                            return gsl_ran_flat(rng, double(i), double(i + 1));
                        }, // mutation pos'n
                        [&effects_dominance, i, rng]() {
                            return effects_dominance[i].s(rng);
                        },
                        [&effects_dominance, i, rng]() {
                            return effects_dominance[i].h(rng);
                        }));
                    ++i;
                }

            // Total mutation rates
            auto tmu(neutral_mutation_rates);
            std::transform(tmu.begin(), tmu.end(),
                           selected_mutation_rates.begin(), tmu.begin(),
                           [](double a, double b) { return a + b; });

            evolve_qtrait_mloc_details_common(
                pop, rng, Nvector, Nvector_len, mmodels, recpols, tmu,
                between_region_rec_rates, fitness, s, interval, f, rules);
            // auto rules_local(std::forward<rules_type>(rules));
            // evolve...
            // const unsigned simlen = unsigned(Nvector_len);
            // fitness->update(pop);
            /*
                        for (unsigned g = 0; g < simlen; ++g,
            ++pop->generation)
                {
                    const unsigned nextN = *(Nvector + g);
                    if (interval && pop->generation
                        && pop->generation % interval == 0.)
                        {
                            s(pop, pop->generation);
                        }
                    KTfwd::experimental::sample_diploid(
                        rng, pop->gametes, pop->diploids, pop->mutations,
                        pop->mcounts, pop->N, nextN, tmu.data(), mmodels,
                        recpols, between_region_rec_rates.data(),
                        // rec b/w loci is interpreted as cM!!!!!
                        [](const gsl_rng *__r, const double __d) {
                            return gsl_ran_bernoulli(__r, __d);
                        },
                        fitness->fitness_function, pop->neutral, pop->selected,
                        f, rules_local, KTfwd::remove_neutral());
                    fwdpy::update_mutations_n(pop->mutations, pop->fixations,
                                              pop->fixation_times,
                                              pop->mut_lookup, pop->mcounts,
                                              pop->generation, 2 * nextN);
                    pop->N = nextN;
                    //fitness->update(pop);
                }
            if (interval && pop->generation
                && pop->generation % interval == 0.)
                {
                    s(pop, pop->generation);
                }
                                */
            gsl_rng_free(rng);
            // Allow a sampler to clean up after itself
            s.cleanup();
        }

		//! \deprecated
        void evolve_qtrait_mloc_cpp(
            GSLrng_t *rng, std::vector<std::shared_ptr<multilocus_t>> *pops,
            std::vector<std::unique_ptr<sampler_base>> &samplers,
            const unsigned *Nvector, const size_t Nvector_length,
            const std::vector<double> &neutral_mutation_rates,
            const std::vector<double> &selected_mutation_rates,
            const std::vector<KTfwd::extensions::shmodel> &shmodels,
            const std::vector<double> &within_region_rec_rates,
            const std::vector<double> &between_region_rec_rates,
            const double f, const double sigmaE, const double optimum,
            const double VS, const int interval,
            const multilocus_fitness &fitness);

		//! Evolve a multi-locus quant-trait system w/"regions"
        void evolve_qtrait_mloc_regions_cpp(
            GSLrng_t *rng, std::vector<std::shared_ptr<multilocus_t>> *pops,
            std::vector<std::unique_ptr<sampler_base>> &samplers,
            const unsigned *Nvector, const size_t Nvector_length,
            const internal::region_manager *rm,
            const std::vector<double> &between_region_rec_rates,
            const double f, const double sigmaE, const double optimum,
            const double VS, const int interval,
            const multilocus_fitness &fitness);
    }
}
#endif
