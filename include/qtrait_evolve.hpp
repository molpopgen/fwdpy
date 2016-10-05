
#ifndef FWDP_QTRAIT_EVOLVE_QTRAIT_SAMPLER_HPP
#define FWDP_QTRAIT_EVOLVE_QTRAIT_SAMPLER_HPP

#include "fwdpp_features.hpp"
#include "fwdpy_fitness.hpp"
#include "internal_region_manager.hpp"
#include "reserve.hpp"
#include "sampler_base.hpp"
#include "types.hpp"
#include <algorithm>
#include <future>
#include <fwdpp/diploid.hh>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/experimental/sample_diploid.hpp>
#include <type_traits>
#include <vector>

namespace fwdpy
{
    namespace qtrait
    {
        template <typename rules_t>
        void
        evolve_regions_qtrait_sampler_cpp_details(
            singlepop_t *pop, const unsigned long seed,
            const unsigned *Nvector, const size_t Nvector_len,
            const double neutral, const double selected, const double recrate,
            const double f, const double sigmaE, const double optimum,
            const double VS, std::unique_ptr<singlepop_fitness> &fitness,
            const int interval, KTfwd::extensions::discrete_mut_model &&__m,
            KTfwd::extensions::discrete_rec_model &&__recmap, sampler_base &s,
            rules_t &&rules)
        /*
          \note the gist of this implementation is from
          fwdpy/fwdpy/evolve_regions_sampler.cc
        */
        {
            gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(rng, seed);
            const unsigned simlen = unsigned(Nvector_len);
            const double mu_tot = neutral + selected;
            auto x = std::max_element(Nvector, Nvector + Nvector_len);
            assert(x != Nvector + Nvector_len);
            reserve_space(pop->gametes, pop->mutations, *x, mu_tot);
            KTfwd::extensions::discrete_mut_model m(std::move(__m));
            KTfwd::extensions::discrete_rec_model recmap(std::move(__recmap));
            rules_t model_rules(std::forward<rules_t>(rules));
            const auto recpos = KTfwd::extensions::bind_drm(
                recmap, pop->gametes, pop->mutations, rng, recrate);
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
                        pop->mcounts, pop->N, nextN, mu_tot,
                        KTfwd::extensions::bind_dmm(
                            m, pop->mutations, pop->mut_lookup, rng, neutral,
                            selected, pop->generation),
                        recpos, fitness->fitness_function, pop->neutral,
                        pop->selected, f, model_rules,
                        KTfwd::remove_neutral());
                    fwdpy::update_mutations_n(pop->mutations, pop->fixations,
                                              pop->fixation_times,
                                              pop->mut_lookup, pop->mcounts,
                                              pop->generation, 2 * nextN);
                    assert(KTfwd::check_sum(pop->gametes, 2 * nextN));
                    pop->N = nextN;
                    // fitness->update(pop);
                }
            if (interval && pop->generation
                && pop->generation % interval == 0.)
                {
                    s(pop, pop->generation);
                }
            gsl_rng_free(rng);
            // Allow a sampler to clean up after itself
            s.cleanup();
        }

        void evolve_regions_qtrait_cpp(
            GSLrng_t *rng, std::vector<std::shared_ptr<singlepop_t>> &pops,
            std::vector<std::unique_ptr<sampler_base>> &samplers,
            const unsigned *Nvector, const size_t Nvector_length,
            const double neutral, const double selected, const double recrate,
            const double f, const double sigmaE, const double optimum,
            const double VS, const int interval,
            const internal::region_manager *rm,
            const singlepop_fitness &fitness);
    }
}

#endif
