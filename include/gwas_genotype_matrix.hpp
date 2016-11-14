#ifndef FWDPY_GWAS_GENOTYPE_MATRIX_HPP
#define FWDPY_GWAS_GENOTYPE_MATRIX_HPP

//#include "sampler_additive_variance.hpp"
#include "types.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_set>
#include <stdexcept>
namespace fwdpy
{
    namespace gwas
    {
        struct genotype_matrix
        /*!
          Intended to be coerced to 2D numpy matrix
          */
        {
            std::vector<double> G, E, neutral, causative, npos, cpos,
                q_neutral, q_causative;
            std::size_t N, n_neutral, n_causative;
            genotype_matrix()
                : G{}, E{}, neutral{}, causative{}, npos{}, cpos{},
                  q_neutral{}, q_causative{}, N{ 0 }, n_neutral{ 0 },
                  n_causative{ 0 }
            {
            }
            genotype_matrix(const genotype_matrix &) = default;
        };

        struct mut_info
        {
            std::vector<std::size_t> neut_indexes, causative_indexes;
            std::vector<double> npos, cpos, qn, qc;
            mut_info()
                : neut_indexes{}, causative_indexes{}, npos{}, cpos{}, qn{},
                  qc{}
            {
            }
        };

        inline void
        update_mut_key_sets(std::unordered_set<std::size_t> &m,
                            const KTfwd::uint_t twoN,
                            const std::vector<KTfwd::uint_t> &mcounts,
                            const std::vector<KTfwd::uint_t> &keys)
        {
            for (auto &k : keys)
                {
                    if (mcounts[k] && mcounts[k] < twoN)
                        {
                            m.emplace(k);
                        }
                }
        }

        inline mut_info
        get_mutation_indexes(const singlepop_t *pop,
                             const std::vector<std::size_t> &individuals)
        {
            mut_info m;
            std::unordered_set<std::size_t> neutral, causative;
            const double twoN = 2. * static_cast<double>(pop->diploids.size());
            for (auto &ind : individuals)
                {
                    auto &dip = pop->diploids[ind];
                    update_mut_key_sets(
                        neutral, static_cast<KTfwd::uint_t>(twoN),
                        pop->mcounts, pop->gametes[dip.first].mutations);
                    update_mut_key_sets(
                        neutral, static_cast<KTfwd::uint_t>(twoN),
                        pop->mcounts, pop->gametes[dip.second].mutations);
                    update_mut_key_sets(
                        causative, static_cast<KTfwd::uint_t>(twoN),
                        pop->mcounts, pop->gametes[dip.first].smutations);
                    update_mut_key_sets(
                        causative, static_cast<KTfwd::uint_t>(twoN),
                        pop->mcounts, pop->gametes[dip.second].smutations);
                }
            m.neut_indexes.insert(m.neut_indexes.end(), neutral.begin(),
                                  neutral.end());
            m.causative_indexes.insert(m.causative_indexes.end(),
                                       causative.begin(), causative.end());
            std::sort(m.neut_indexes.begin(), m.neut_indexes.end(),
                      [&pop](std::size_t a, std::size_t b) {
                          return pop->mutations[a].pos < pop->mutations[b].pos;
                      });
            std::sort(m.causative_indexes.begin(), m.causative_indexes.end(),
                      [&pop](std::size_t a, std::size_t b) {
                          return pop->mutations[a].pos < pop->mutations[b].pos;
                      });
            for (auto &mk : m.neut_indexes)
                {
                    m.npos.push_back(pop->mutations[mk].pos);
                    m.qn.push_back(static_cast<double>(pop->mcounts[mk])
                                   / twoN);
                }
            for (auto &mk : m.causative_indexes)
                {
                    m.cpos.push_back(pop->mutations[mk].pos);
                    m.qc.push_back(static_cast<double>(pop->mcounts[mk])
                                   / twoN);
                }
            return m;
        }

        inline void
        update_row(std::vector<double> &v,
                   const std::vector<KTfwd::uint_t> &mut_keys,
                   const std::vector<std::size_t> &indexes)
        {
            if (v.size() != indexes.size())
                {
                    throw std::runtime_error("vector sizes do not match");
                }
            for (auto &&mk : mut_keys)
                {
                    auto i = std::find(indexes.begin(), indexes.end(),
                                       static_cast<std::size_t>(mk));
                    if (i != indexes.end()) // i may equal indexes.end iff mk
                                            // refers to a fixation
                        {
                            auto idx = std::distance(indexes.begin(), i);
                            if (static_cast<std::size_t>(idx) >= v.size())
                                {
                                    throw std::runtime_error(
                                        "idx >= v.size()");
                                }
                            v[static_cast<std::size_t>(idx)] += 1.0;
                        }
                }
        }

        genotype_matrix
        make_geno_matrix(const singlepop_t *pop,
                         std::vector<std::size_t> &individuals,
                         bool maf = true)
        //! Quick and dirty for ADL & Ted.
        {
            auto mi = get_mutation_indexes(pop, individuals);
            genotype_matrix rv;
            rv.G.reserve(individuals.size());
            rv.E.reserve(individuals.size());
            std::vector<double> neutral_row(mi.neut_indexes.size(), 0.);
            std::vector<double> causative_row(mi.causative_indexes.size(), 0.);
            for (auto ind : individuals)
                {
                    auto &dip = pop->diploids[ind];
                    rv.G.push_back(dip.g);
                    rv.E.push_back(dip.e);
                    update_row(neutral_row, pop->gametes[dip.first].mutations,
                               mi.neut_indexes);
                    update_row(neutral_row, pop->gametes[dip.second].mutations,
                               mi.neut_indexes);
                    update_row(causative_row,
                               pop->gametes[dip.first].smutations,
                               mi.causative_indexes);
                    update_row(causative_row,
                               pop->gametes[dip.second].smutations,
                               mi.causative_indexes);
                    rv.neutral.insert(rv.neutral.end(), neutral_row.begin(),
                                      neutral_row.end());
                    rv.causative.insert(rv.causative.end(),
                                        causative_row.begin(),
                                        causative_row.end());
                    std::fill(neutral_row.begin(), neutral_row.end(), 0.);
                    std::fill(causative_row.begin(), causative_row.end(), 0.);
                }
            rv.n_neutral = mi.neut_indexes.size();
            rv.n_causative = mi.causative_indexes.size();
            rv.npos = std::move(mi.npos);
            rv.cpos = std::move(mi.cpos);
            rv.q_neutral = std::move(mi.qn);
            rv.q_causative = std::move(mi.qc);
            rv.N = individuals.size();
            double check = 0.;
            for (std::size_t i = 1; i < rv.neutral.size() / 2;
                 i += rv.n_neutral)
                {
                    check += rv.neutral[i];
                }
            // Final checks on rv
            if (rv.neutral.size() != (individuals.size() * rv.n_neutral))
                {
                    throw std::runtime_error(
                        "incorrect matrix size for neutral variants");
                }
            if (rv.causative.size() != (individuals.size() * rv.n_causative))
                {
                    throw std::runtime_error(
                        "incorrect matrix size for causative variants");
                }
            if (rv.npos.size() != rv.q_neutral.size()
                || rv.npos.size() != rv.n_neutral
                || rv.q_neutral.size() != rv.n_neutral)
                {
                    throw std::runtime_error(
                        "neutral: " + std::to_string(rv.n_neutral) + ' '
                        + std::to_string(rv.npos.size()) + ' '
                        + std::to_string(rv.q_neutral.size()));
                }
            return rv;
        }

        genotype_matrix
        make_geno_matrix(const singlepop_t *pop, bool maf = true)
        {
            std::vector<std::size_t> ind(pop->diploids.size());
            std::size_t i = 0;
            std::transform(ind.begin(), ind.end(), ind.begin(),
                           [&i](std::size_t indi) { return i++; });
            return make_geno_matrix(pop, ind, maf);
        }
    }
}

#endif
