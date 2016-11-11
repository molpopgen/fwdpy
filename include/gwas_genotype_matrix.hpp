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
namespace fwdpy
{
    namespace gwas
    {
        struct genotype_matrix
        /*!
          Intended to be coerced to 2D numpy matrix
        */
        {
            std::vector<double> G, E, neutral, causative, npos, cpos;
            std::size_t N, n_neutral, n_causative;
            genotype_matrix(std::vector<double> &&G_, std::vector<double> &&E_,
                            std::vector<double> &&n_, std::vector<double> &&c_,
                            std::vector<double> &&np_,
                            std::vector<double> &&cp_, std::size_t N_,
                            std::size_t nn_, std::size_t ns_)
                : G(std::move(G_)), E(std::move(E_)), neutral(std::move(n_)),
                  causative(std::move(c_)), npos(std::move(np_)),
                  cpos(std::move(cp_)), N(N_), n_neutral(nn_), n_causative(ns_)
            {
            }
            genotype_matrix() {}
            genotype_matrix(const genotype_matrix &) = default;
        };

        struct mut_info
        {
            std::vector<std::size_t> neut_indexes, causative_indexes;
            std::vector<double> npos, cpos;
            mut_info() : neut_indexes{}, causative_indexes{}, npos{}, cpos{} {}
        };

        inline mut_info
        get_mutation_indexes(const singlepop_t *pop,
                             const std::vector<std::size_t> &individuals)
        {
            mut_info m;
            std::unordered_set<std::size_t> neut_indexes, causative_indexes;
            std::unordered_set<double> npos, cpos;
            for (auto &ind : individuals)
                {
                    auto &dip = pop->diploids[ind];
                    for (auto &m : pop->gametes[dip.first].mutations)
                        {
                            if (pop->mcounts[m]
                                && pop->mcounts[m] < 2 * pop->diploids.size())
                                {
                                    neut_indexes.insert(m);
                                    npos.insert(pop->mutations[m].pos);
                                }
                        }
                    for (auto &m : pop->gametes[dip.second].mutations)
                        {
                            if (pop->mcounts[m]
                                && pop->mcounts[m] < 2 * pop->diploids.size())
                                {
                                    neut_indexes.insert(m);
                                    npos.insert(pop->mutations[m].pos);
                                }
                        }
                    for (auto &m : pop->gametes[dip.first].smutations)
                        {
                            if (pop->mcounts[m]
                                && pop->mcounts[m] < 2 * pop->diploids.size())
                                {
                                    causative_indexes.insert(m);
                                    cpos.insert(pop->mutations[m].pos);
                                }
                        }
                    for (auto &m : pop->gametes[dip.second].smutations)
                        {
                            if (pop->mcounts[m]
                                && pop->mcounts[m] < 2 * pop->diploids.size())
                                {
                                    causative_indexes.insert(m);
                                    cpos.insert(pop->mutations[m].pos);
                                }
                        }
                }
            m.neut_indexes.insert(m.neut_indexes.end(), neut_indexes.begin(),
                                  neut_indexes.end());
            m.causative_indexes.insert(m.causative_indexes.end(),
                                       causative_indexes.begin(),
                                       causative_indexes.end());
			m.npos.insert(m.npos.end(),npos.begin(),npos.end());
			m.cpos.insert(m.cpos.end(),cpos.begin(),cpos.end());
            /*for (std::size_t i = 0; i < pop->mcounts.size(); ++i)
                {
                    if (pop->mcounts[i]
                        && pop->mcounts[i] < 2 * pop->diploids.size())
                        {
                            if (pop->mutations[i].neutral)
                                {
                                    m.neut_indexes.push_back(i);
                                    m.npos.push_back(pop->mutations[i].pos);
                                }
                            else
                                {
                                    m.causative_indexes.push_back(i);
                                    m.cpos.push_back(pop->mutations[i].pos);
                                }
                        }
                }
                                */
            return m;
        }

        genotype_matrix
        make_geno_matrix(const singlepop_t *pop,
                         std::vector<std::size_t> &individuals,
                         bool maf = true)
        //! Quick and dirty for ADL & Ted.
        {
            auto mi = get_mutation_indexes(pop,individuals);
            // Use GSL matrixes
            gsl::gsl_matrix_ptr_t gn(gsl_matrix_alloc(pop->diploids.size(),
                                                      mi.neut_indexes.size())),
                gc(gsl_matrix_alloc(pop->diploids.size(),
                                    mi.causative_indexes.size()));
            gsl_matrix_set_zero(gn.get());
            gsl_matrix_set_zero(gc.get());
            // Fill the matrices, etc.
            std::size_t row = 0;
            std::vector<double> G, E;
            G.reserve(pop->diploids.size());
            E.reserve(pop->diploids.size());
            // for (const auto &dip : pop->diploids)
            for (auto ind : individuals)
                {
                    auto &dip = pop->diploids[ind];
                    G.push_back(dip.g);
                    E.push_back(dip.e);
                    // Fill selected matrix
                    if (!mi.causative_indexes.empty())
                        {
                            for (auto k : pop->gametes[dip.first].smutations)
                                {
                                    auto i = std::find(
                                        mi.causative_indexes.begin(),
                                        mi.causative_indexes.end(), k);
                                    if (i != mi.causative_indexes.end())
                                        {
                                            std::size_t col = std::distance(
                                                mi.causative_indexes.begin(),
                                                i);
                                            if (col >= gc->size2)
                                                throw std::out_of_range(
                                                    "column index out of "
                                                    "range: "
                                                    + std::to_string(col)
                                                    + ">=" + std::to_string(
                                                                 gc->size2));
                                            auto m = gsl_matrix_ptr(gc.get(),
                                                                    row, col);
                                            *m += 1.0;
                                        }
                                }
                            for (auto k : pop->gametes[dip.second].smutations)
                                {
                                    auto i = std::find(
                                        mi.causative_indexes.begin(),
                                        mi.causative_indexes.end(), k);
                                    if (i != mi.causative_indexes.end())
                                        {
                                            std::size_t col = std::distance(
                                                mi.causative_indexes.begin(),
                                                i);
                                            if (col >= gc->size2)
                                                throw std::out_of_range(
                                                    "column index out of "
                                                    "range: "
                                                    + std::to_string(col)
                                                    + ">=" + std::to_string(
                                                                 gc->size2));
                                            auto m = gsl_matrix_ptr(gc.get(),
                                                                    row, col);
                                            *m += 1.0;
                                        }
                                }
                        }
                    // Fill neutral matrix
                    if (!mi.neut_indexes.empty())
                        {
                            for (auto k : pop->gametes[dip.first].mutations)
                                {
                                    auto i
                                        = std::find(mi.neut_indexes.begin(),
                                                    mi.neut_indexes.end(), k);
                                    if (i != mi.neut_indexes.end())
                                        {
                                            std::size_t col = std::distance(
                                                mi.neut_indexes.begin(), i);
                                            if (col >= gn->size2)
                                                throw std::out_of_range(
                                                    "column index out of "
                                                    "range: "
                                                    + std::to_string(col)
                                                    + ">=" + std::to_string(
                                                                 gn->size2));
                                            auto m = gsl_matrix_ptr(gn.get(),
                                                                    row, col);
                                            *m += 1.0;
                                        }
                                }
                            for (auto k : pop->gametes[dip.second].mutations)
                                {
                                    auto i
                                        = std::find(mi.neut_indexes.begin(),
                                                    mi.neut_indexes.end(), k);
                                    if (i != mi.neut_indexes.end())
                                        {
                                            std::size_t col = std::distance(
                                                mi.neut_indexes.begin(), i);
                                            if (col >= gn->size2)
                                                throw std::out_of_range(
                                                    "column index out of "
                                                    "range: "
                                                    + std::to_string(col)
                                                    + ">=" + std::to_string(
                                                                 gn->size2));
                                            auto m = gsl_matrix_ptr(gn.get(),
                                                                    row, col);
                                            *m += 1.0;
                                        }
                                }
                        }
                    ++row;
                }
            return genotype_matrix(
                std::move(G), std::move(E),
                (!mi.neut_indexes.empty())
                    ? std::vector<double>(gn->data,
                                          gn->data + (gn->size1 * gn->size2))
                    : std::vector<double>(),
                (!mi.causative_indexes.empty())
                    ? std::vector<double>(gc->data,
                                          gc->data + (gc->size1 * gc->size2))
                    : std::vector<double>(),
                std::move(mi.npos), std::move(mi.cpos), pop->diploids.size(),
                mi.neut_indexes.size(), mi.causative_indexes.size());
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
