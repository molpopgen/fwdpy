#ifndef FWDPY_GWAS_GENOTYPE_MATRIX_HPP
#define FWDPY_GWAS_GENOTYPE_MATRIX_HPP

#include "sampler_additive_variance.hpp"
#include "gsl.hpp"
#include "types.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_set>
#include <stdexcept>
#include <iostream>
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
            genotype_matrix() {}
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
            for (auto &ind : individuals)
                {
                    auto &dip = pop->diploids[ind];
                    update_mut_key_sets(neutral, 2 * pop->diploids.size(),
                                        pop->mcounts,
                                        pop->gametes[dip.first].mutations);
                    update_mut_key_sets(neutral, 2 * pop->diploids.size(),
                                        pop->mcounts,
                                        pop->gametes[dip.second].mutations);
                    update_mut_key_sets(causative, 2 * pop->diploids.size(),
                                        pop->mcounts,
                                        pop->gametes[dip.first].smutations);
                    update_mut_key_sets(causative, 2 * pop->diploids.size(),
                                        pop->mcounts,
                                        pop->gametes[dip.second].smutations);
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
            const double twoN = 2. * static_cast<double>(pop->diploids.size());
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
        update_matrix(gsl_matrix *m, const std::size_t row,
                      const std::vector<std::size_t> &indexes,
                      const std::vector<KTfwd::uint_t> &mut_keys)
        {
            for (auto &&k : mut_keys)
                {
                    auto i = std::find(indexes.begin(), indexes.end(), k);
                    if (i == indexes.end())
                        {
                            throw std::runtime_error(
                                "mutation key not found in indexes");
                        }
                    auto col = std::distance(indexes.begin(), i);
                    auto p = gsl_matrix_ptr(m, row, col);
                    *p += 1.0;
                }
        }

        inline std::vector<double>
        row_major_copy(const gsl_matrix *m)
        {
            std::vector<double> rv;
            for (std::size_t row = 0; row < m->size1; ++row)
                {
                    auto row_view = gsl_matrix_const_row(m, row);
                    for (std::size_t j = 0; j < row_view.vector.size; ++j)
                        {
                            rv.push_back(gsl_vector_get(&row_view.vector, j));
                        }
                }
            return rv;
        }

        genotype_matrix
        make_geno_matrix(const singlepop_t *pop,
                         std::vector<std::size_t> &individuals,
                         bool maf = true)
        //! Quick and dirty for ADL & Ted.
        {
            auto mi = get_mutation_indexes(pop, individuals);
            // Use GSL matrixes
            gsl::gsl_matrix_ptr_t gn(
                gsl_matrix_alloc(individuals.size(), mi.neut_indexes.size())),
                gc(gsl_matrix_alloc(individuals.size(),
                                    mi.causative_indexes.size()));
            gsl_matrix_set_zero(gn.get());
            gsl_matrix_set_zero(gc.get());
            // Fill the matrices, etc.
            std::size_t row = 0;
            genotype_matrix rv;
            rv.G.reserve(individuals.size());
            rv.E.reserve(individuals.size());
            for (auto ind : individuals)
                {
                    auto &dip = pop->diploids[ind];
                    rv.G.push_back(dip.g);
                    rv.E.push_back(dip.e);
                    // Fill selected matrix
                    update_matrix(gc.get(), row, mi.causative_indexes,
                                  pop->gametes[dip.first].smutations);
                    update_matrix(gc.get(), row, mi.causative_indexes,
                                  pop->gametes[dip.second].smutations);
                    // Fill neutral matrix
                    update_matrix(gn.get(), row, mi.neut_indexes,
                                  pop->gametes[dip.first].mutations);
                    update_matrix(gn.get(), row, mi.neut_indexes,
                                  pop->gametes[dip.second].mutations);
                    ++row;
                }
            // sum up first 3,000 of column 0
            auto c0 = gsl_matrix_const_column(gn.get(), 0);
            double c0sum = 0.;
            for (unsigned j = 0; j < c0.vector.size; ++j)
                {
                    c0sum += gsl_vector_get(&c0.vector, j);
                }
            std::cout << "c0sum = " << c0sum << '\n';
            for (unsigned j = 0; j < gn->size2; ++j)
                {
                    // check that this column is sound
                    auto c = gsl_matrix_const_column(gn.get(), j);
                    double sum = 0.;
                    for (unsigned k = 0; k < c.vector.size; ++k)
                        {
                            sum += gsl_vector_get(&c.vector, k);
                        }
                    unsigned ncopies = static_cast<unsigned>(sum);
                    unsigned ncopies_pop = static_cast<unsigned>(std::round(
                        mi.qn[j] * 2.
                        * static_cast<double>(pop->diploids.size())));
                    unsigned ncopies_pop2 = pop->mcounts[mi.neut_indexes[j]];
                    if (ncopies > ncopies_pop || ncopies > ncopies_pop2)
                        {
                            throw std::runtime_error(
                                std::to_string(ncopies) + ' '
                                + std::to_string(ncopies_pop) + ' '
                                + std::to_string(ncopies_pop2) + ' '
                                + std::to_string(sum) + ' '
                                + std::to_string(mi.qn[j]));
                        }
                }
            rv.neutral = row_major_copy(gn.get());
            rv.causative = row_major_copy(gc.get());
            rv.n_neutral = mi.neut_indexes.size();
            rv.n_causative = mi.causative_indexes.size();
            rv.npos = std::move(mi.npos);
            rv.cpos = std::move(mi.cpos);
            rv.q_neutral = std::move(mi.qn);
            rv.q_causative = std::move(mi.qc);
            rv.N = individuals.size();
            // Final checks on rv
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
