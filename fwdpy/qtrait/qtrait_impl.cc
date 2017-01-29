/*
  Models of quantitative traits
  Trait values are additive over 1, 1+hs, 1+2s, where s is a Gaussian deviate

  The infinitely-many sites stuff is an Cython/fwdpp-based re-implementation of
  the
  code used to generate preliminary data for R01GM115564.
*/

#include <algorithm>
#include <future>
#include <fwdpp/diploid.hh>
#include <fwdpp/experimental/sample_diploid.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <gsl/gsl_statistics_double.h>
#include <limits>
#include <memory>
#include <thread>
#include <utility>
#include <vector>

#include "internal_region_manager.hpp"
#include "qtrait_details.hpp"
#include "qtrait_evolve.hpp"
#include "qtrait_evolve_rules.hpp"
#include "sampler_additive_variance.hpp"
#include "sampler_no_sampling.hpp"
#include "types.hpp"

using namespace std;

using poptype = fwdpy::singlepop_t;

namespace fwdpy
{
    namespace qtrait
    {
        void
        evolve_regions_qtrait_cpp(
            GSLrng_t *rng, std::vector<std::shared_ptr<singlepop_t>> &pops,
            std::vector<std::unique_ptr<sampler_base>> &samplers,
            const unsigned *Nvector, const size_t Nvector_length,
            const double neutral, const double selected, const double recrate,
            const double f, const double sigmaE, const double optimum,
            const double VS, const int interval,
            const internal::region_manager *rm,
            const singlepop_fitness &fitness)
        {
            if (neutral < 0. || selected < 0. || recrate < 0.)
                {
                    throw std::runtime_error("mutation and recombination "
                                             "rates must all be "
                                             "non-negative.");
                }
            if (samplers.size() != pops.size())
                {
                    throw std::runtime_error("length of samplers != length of "
                                             "population container");
                }
            if (f < 0. || f > 1.)
                throw std::runtime_error(
                    "selfing probabilty must be 0<=f<=1.");
            if (interval < 0)
                throw std::runtime_error(
                    "sampling interval must be non-negative");
            qtrait_model_rules rules(
                sigmaE, optimum, VS,
                *std::max_element(Nvector, Nvector + Nvector_length));
            if (pops.size() > 1)
                {
                    std::vector<std::thread> threads;
                    std::vector<std::unique_ptr<singlepop_fitness>> fitnesses;
                    for (std::size_t i = 0; i < pops.size(); ++i)
                        {
                            fitnesses.emplace_back(
                                std::unique_ptr<singlepop_fitness>(
                                    fitness.clone()));
                        }
                    for (std::size_t i = 0; i < pops.size(); ++i)
                        {
                            threads.emplace_back(std::thread(
                                evolve_regions_qtrait_sampler_cpp_details<qtrait_model_rules>,
                                pops[i].get(), gsl_rng_get(rng->get()),
                                Nvector, Nvector_length, neutral, selected,
                                recrate, f, sigmaE, optimum, VS,
                                std::ref(fitnesses[i]), interval,
                                KTfwd::extensions::discrete_mut_model(
                                    rm->nb, rm->ne, rm->nw, rm->sb, rm->se,
                                    rm->sw, rm->callbacks),
                                KTfwd::extensions::discrete_rec_model(
                                    rm->rb, rm->rw, rm->rw),
                                std::ref(*samplers[i]), rules));
                        }
                    for (auto &t : threads)
                        t.join();
                }
            else
                {
                    auto fitness_ptr
                        = std::unique_ptr<singlepop_fitness>(fitness.clone());
                    evolve_regions_qtrait_sampler_cpp_details<qtrait_model_rules>(
                        pops[0].get(), gsl_rng_get(rng->get()), Nvector,
                        Nvector_length, neutral, selected, recrate, f, sigmaE,
                        optimum, VS, std::ref(fitness_ptr), interval,
                        KTfwd::extensions::discrete_mut_model(
                            rm->nb, rm->ne, rm->nw, rm->sb, rm->se, rm->sw,
                            rm->callbacks),
                        KTfwd::extensions::discrete_rec_model(rm->rb, rm->rw,
                                                              rm->rw),
                        std::ref(*samplers[0]), std::move(rules));
                }
        }
    } // ns qtrait
} // ns fwdpy
