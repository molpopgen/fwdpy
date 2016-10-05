#include <algorithm>
#include <functional>
#include <future>
#include <fwdpp/diploid.hh>
#include <fwdpp/experimental/sample_diploid.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <iterator>
#include <memory>
#include <thread>
#include <type_traits>
#include <vector>

#include "evolve_regions_sampler.hpp"
#include "fwdpy_fitness.hpp"
#include "reserve.hpp"
#include "sampler_base.hpp"
#include "types.hpp"
#include "wf_rules.hpp"

using namespace std;

namespace fwdpy
{
    void
    evolve_regions_sampler_cpp_details(
        singlepop_t *pop, const unsigned long seed, const unsigned *Nvector,
        const size_t Nvector_len, const double neutral, const double selected,
        const double recrate, const double f,
        std::unique_ptr<singlepop_fitness> &fitness, const int interval,
        KTfwd::extensions::discrete_mut_model &&__m,
        KTfwd::extensions::discrete_rec_model &&__recmap, sampler_base &s,
        wf_rules rules)
    {
        const size_t simlen = Nvector_len;
        auto x = std::max_element(Nvector, Nvector + Nvector_len);
        assert(x != Nvector + Nvector_len);
        reserve_space(pop->gametes, pop->mutations, *x, neutral + selected);
        const double mu_tot = neutral + selected;
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        KTfwd::extensions::discrete_mut_model m(std::move(__m));
        KTfwd::extensions::discrete_rec_model recmap(std::move(__recmap));
        // Recombination policy: more complex than the standard case...
        const auto recpos = KTfwd::extensions::bind_drm(
            recmap, pop->gametes, pop->mutations, rng, recrate);

        wf_rules local_rules(std::move(rules));
        /*
          Update fitness model data.
          Needed for stateful fitness models and
          situtations where we evolve, stop, then
          evolve again.
        */
        // fitness->update(pop);
        for (size_t g = 0; g < simlen; ++g, ++pop->generation)
            {
                const unsigned nextN = *(Nvector + g);
                KTfwd::experimental::sample_diploid(
                    rng, pop->gametes, pop->diploids, pop->mutations,
                    pop->mcounts, pop->N, nextN, mu_tot,
                    KTfwd::extensions::bind_dmm(m, pop->mutations,
                                                pop->mut_lookup, rng, neutral,
                                                selected, pop->generation),
                    recpos, fitness->fitness_function, pop->neutral,
                    pop->selected, f, local_rules);
                pop->N = nextN;
                if (interval && pop->generation + 1
                    && (pop->generation + 1) % interval == 0.)
                    {
                        s(pop, pop->generation + 1);
                    }
                KTfwd::update_mutations(
                    pop->mutations, pop->fixations, pop->fixation_times,
                    pop->mut_lookup, pop->mcounts, pop->generation, 2 * nextN);
                // Allow fitness model to update any data that it may need
                // fitness->update(pop);
                assert(KTfwd::check_sum(pop->gametes, 2 * nextN));
            }
        // if (interval && pop->generation && (pop->generation) % interval ==
        // 0.)
        //    {
        //        s(pop, pop->generation);
        //    }
        // Update population's size variable to be the current pop size
        pop->N = unsigned(pop->diploids.size());
        // cleanup
        gsl_rng_free(rng);
        // Let the sampler clean up after itself
        s.cleanup();
    }

    void
    evolve_regions_sampler_cpp(
        GSLrng_t *rng, std::vector<std::shared_ptr<singlepop_t>> &pops,
        std::vector<std::unique_ptr<sampler_base>> &samplers,
        const unsigned *Nvector, const size_t Nvector_length,
        const double mu_neutral, const double mu_selected,
        const double littler, const double f, const int sample,
        const internal::region_manager *rm, const singlepop_fitness &fitness)
    {
        // check inputs--this is point of failure.  Throw excceptions here b4
        // getting into any threaded nonsense.
        if (mu_neutral < 0. || mu_selected < 0. || littler < 0.)
            {
                throw std::runtime_error("mutation and recombination rates "
                                         "must all be non-negative.");
            }
        if (f < 0. || f > 1.)
            throw std::runtime_error("selfing probabilty must be 0<=f<=1.");
        if (sample < 0)
            throw std::runtime_error("sampling interval must be non-negative");
        std::vector<std::thread> threads;
        wf_rules rules;
        std::vector<std::unique_ptr<singlepop_fitness>> fitnesses;
        for (std::size_t i = 0; i < pops.size(); ++i)
            {
                fitnesses.emplace_back(
                    std::unique_ptr<singlepop_fitness>(fitness.clone()));
            }
        for (std::size_t i = 0; i < pops.size(); ++i)
            {
                threads.emplace_back(std::thread(
                    evolve_regions_sampler_cpp_details,
                    pops[i].get(), gsl_rng_get(rng->get()),
                    Nvector, Nvector_length, mu_neutral, mu_selected, littler,
                    f, std::ref(fitnesses[i]), sample,
                    KTfwd::extensions::discrete_mut_model(
                        rm->nb, rm->ne, rm->nw, rm->sb, rm->se, rm->sw,
                        rm->callbacks),
                    KTfwd::extensions::discrete_rec_model(rm->rb, rm->rw,
                                                          rm->rw),
                    std::ref(*samplers[i]), rules));
            }
        for (auto &t : threads)
            t.join();
    }
}
