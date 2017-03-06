#include <algorithm>
#include <future>
#include <fwdpp/diploid.hh>
#include <fwdpp/experimental/sample_diploid.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <gsl/gsl_statistics_double.h>
#include <limits>
#include <memory>
#include <set>
#include <thread>
#include <utility>
#include <vector>

#include "qtrait_evolve_mlocus.hpp"
#include "qtrait_mloc_rules.hpp"
#include "types.hpp"

using namespace std;

namespace fwdpy
{
    namespace qtrait
    {
        void
        evolve_qtrait_mloc_cpp(
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
            const multilocus_fitness &fitness)
        {
            std::set<std::size_t> vec_sizes{ neutral_mutation_rates.size(),
                                             selected_mutation_rates.size(),
                                             shmodels.size(),
                                             within_region_rec_rates.size() };
            if (samplers.size() != pops->size())
                {
                    throw std::runtime_error("length of samplers != length of "
                                             "population container");
                }
            if (vec_sizes.size() > 1)
                {
                    throw std::runtime_error("vectors of properties for each "
                                             "region must be same length");
                }
            if (std::any_of(neutral_mutation_rates.begin(),
                            neutral_mutation_rates.end(),
                            [](double d) { return d < 0.; }))
                {
                    throw std::runtime_error(
                        "neutral mutation rates must be >= 0 for all loci");
                }
            if (f < 0. || f > 1.)
                {
                    throw std::runtime_error(
                        "selfing probabilty must be 0<=f<=1.");
                }
            if (interval < 0)
                {
                    throw std::runtime_error(
                        "sampling interval must be non-negative");
                }
            qtrait_mloc_rules rules(
                sigmaE, optimum, VS,
                *std::max_element(Nvector, Nvector + Nvector_length));
            if (pops->size() > 1)
                {
                    std::vector<std::thread> threads;
                    std::vector<std::unique_ptr<multilocus_fitness>> fitnesses;
                    for (std::size_t i = 0; i < pops->size(); ++i)
                        {
                            fitnesses.emplace_back(
                                std::unique_ptr<multilocus_fitness>(
                                    fitness.clone()));
                        }
                    for (std::size_t i = 0; i < pops->size(); ++i)
                        {
                            threads.emplace_back(std::thread(
                                evolve_qtrait_mloc_cpp_details<qtrait_mloc_rules>,
                                pops->operator[](i).get(),
                                std::ref(fitnesses[i]), std::ref(*samplers[i]),
                                gsl_rng_get(rng->get()), Nvector,
                                Nvector_length, neutral_mutation_rates,
                                selected_mutation_rates, shmodels,
                                within_region_rec_rates,
                                between_region_rec_rates, f, interval, rules));
                        }
                    for (auto &t : threads)
                        {
                            t.join();
                        }
                }
            else
                {
                    auto fitness_clone
                        = std::unique_ptr<multilocus_fitness>(fitness.clone());
                    evolve_qtrait_mloc_cpp_details<qtrait_mloc_rules>(
                        pops->operator[](0).get(), fitness_clone, *samplers[0],
                        gsl_rng_get(rng->get()), Nvector, Nvector_length,
                        neutral_mutation_rates, selected_mutation_rates,
                        shmodels, within_region_rec_rates,
                        between_region_rec_rates, f, interval,
                        std::move(rules));
                }
        }

        void
        evolve_qtrait_mloc_regions_cpp(
            GSLrng_t *rng, std::vector<std::shared_ptr<multilocus_t>> *pops,
            std::vector<std::unique_ptr<sampler_base>> &samplers,
            const unsigned *Nvector, const size_t Nvector_length,
            const internal::region_manager *rm,
            const std::vector<double> &between_region_rec_rates,
            const double f, const double sigmaE, const double optimum,
            const double VS, const int interval,
            const multilocus_fitness &fitness)
        {
            if (samplers.size() != pops->size())
                {
                    throw std::runtime_error("length of samplers != length of "
                                             "population container");
                }
            if (f < 0. || f > 1.)
                {
                    throw std::runtime_error(
                        "selfing probabilty must be 0<=f<=1.");
                }
            if (interval < 0)
                {
                    throw std::runtime_error(
                        "sampling interval must be non-negative");
                }
            qtrait_mloc_rules rules(
                sigmaE, optimum, VS,
                *std::max_element(Nvector, Nvector + Nvector_length));
            if (pops->size() > 1)
                {
                    std::vector<std::thread> threads;
                    std::vector<std::unique_ptr<multilocus_fitness>> fitnesses;
                    for (std::size_t i = 0; i < pops->size(); ++i)
                        {
                            fitnesses.emplace_back(
                                std::unique_ptr<multilocus_fitness>(
                                    fitness.clone()));
                        }
                    for (std::size_t i = 0; i < pops->size(); ++i)
                        {
                            threads.emplace_back(std::thread(
                                evolve_qtrait_mloc_regions_cpp_details<qtrait_mloc_rules>,
                                pops->operator[](i).get(),
                                std::ref(fitnesses[i]), std::ref(*samplers[i]),
                                gsl_rng_get(rng->get()), Nvector,
                                Nvector_length, rm, between_region_rec_rates,
                                f, interval, rules));
                        }
                    for (auto &t : threads)
                        {
                            t.join();
                        }
                }
            else
                {
                    auto fitness_clone
                        = std::unique_ptr<multilocus_fitness>(fitness.clone());
                    evolve_qtrait_mloc_regions_cpp_details<qtrait_mloc_rules>(
                        pops->operator[](0).get(), fitness_clone, *samplers[0],
                        gsl_rng_get(rng->get()), Nvector, Nvector_length, rm,
                        between_region_rec_rates, f, interval,
                        std::move(rules));
                }
        }
    }
}
