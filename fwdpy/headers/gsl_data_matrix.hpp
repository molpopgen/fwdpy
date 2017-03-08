#ifndef FWDPY_GSL_DATA_MATRIX_HPP_
#define FWDPY_GSL_DATA_MATRIX_HPP_

#include <fwdpp/forward_types.hpp>
#include <zlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <memory>
#include <cstddef>
#include <sstream>
#include "types.hpp"
#include "gsl.hpp"
namespace fwdpy
{
    namespace gsl_data_matrix
    {
        struct geno_matrix
        {
            std::vector<double> G;
            gsl::gsl_matrix_ptr_t m;
            std::size_t ncol, nrow;
            geno_matrix(std::size_t nrow_, std::size_t ncol_)
                : G{}, m{ gsl_matrix_alloc(nrow_, ncol_) }, ncol{ ncol_ },
                  nrow{ nrow_ }
            {
            }
        };

        inline void
        emplace_move(
            std::vector<std::pair<KTfwd::uint_t, std::unique_ptr<geno_matrix>>>
                &v,
            std::pair<KTfwd::uint_t, std::unique_ptr<geno_matrix>> &p)
        {
            v.emplace_back(std::move(p));
        }

        inline void
        write_geno_matrix(const geno_matrix *m, const KTfwd::uint_t generation,
                          std::string stub, const int repid,
                          const bool keep_origin)
        {
            stub += ".generation" + std::to_string(generation) + ".rep"
                    + std::to_string(repid) + ".gz";
            gzFile gzout = gzopen(stub.c_str(), "w");
            std::ostringstream buffer;
            int nwrites = 0;
            for (std::size_t row = 0; row < m->nrow; ++row)
                {
                    auto row_view = gsl_matrix_const_row(m->m.get(), row);
                    buffer << m->G[row] << '\t';
                    for (std::size_t col
                         = 0 + static_cast<size_t>(keep_origin == false);
                         col < m->ncol; ++col)
                        {
                            buffer << gsl_vector_get(&row_view.vector, col);
                            if (col < m->ncol - 1)
                                buffer << '\t';
                        }
                    buffer << '\n';
                    ++nwrites;
                    if (nwrites == 10)
                        {
                            gzwrite(gzout, buffer.str().c_str(),
                                    buffer.str().size());
                            buffer.str(std::string());
                            nwrites = 0;
                        }
                }
            if (nwrites)
                {
                    gzwrite(gzout, buffer.str().c_str(), buffer.str().size());
                }
            gzclose(gzout);
        }

        template <typename pop_t>
        std::vector<KTfwd::uint_t>
        get_mut_keys(const pop_t *pop, const bool sort_freq = false,
                     const bool sort_esize = false)
        {
            std::vector<KTfwd::uint_t> mut_keys; // array of keys for each
                                                 // segregating, non-neutral
                                                 // variant
            std::set<KTfwd::uint_t>
                ucounts; // the number of unique frequency bins...
            for (std::size_t i = 0; i < pop->mutations.size(); ++i)
                {
                    // first check is to avoid extinct variants that fwdpp will
                    // recycle later.
                    // The second avoids fixed variants
                    if (pop->mcounts[i]
                        && (pop->mcounts[i] < 2 * pop->diploids.size())
                        && !pop->mutations[i].neutral)
                        {
                            mut_keys.push_back(i);
                            ucounts.insert(pop->mcounts[i]);
                        }
                }
            if (sort_freq)
                {
                    // Now, I need to sort based on frequency, descending order
                    std::sort(mut_keys.begin(), mut_keys.end(),
                              [&pop](KTfwd::uint_t a, KTfwd::uint_t b) {
                                  return pop->mcounts[a] > pop->mcounts[b];
                              });
                    if (sort_esize)
                        {
                            // Within each frequency class, sort in descending
                            // order
                            // via
                            // |effect size|...
                            for (auto uc : ucounts)
                                {
                                    auto itr_b = std::find_if(
                                        mut_keys.begin(), mut_keys.end(),
                                        [&pop, uc](KTfwd::uint_t a) {
                                            return pop->mcounts[a] == uc;
                                        });
                                    auto itr_e = std::find_if(
                                        itr_b + 1, mut_keys.end(),
                                        [&pop, uc](KTfwd::uint_t a) {
                                            return pop->mcounts[a] != uc;
                                        });
                                    std::sort(
                                        itr_b, itr_e, [&pop](KTfwd::uint_t a,
                                                             KTfwd::uint_t b) {
                                            return std::fabs(
                                                       pop->mutations[a].s)
                                                   > std::fabs(
                                                         pop->mutations[b].s);
                                        });
                                }
                        }
                }
            else if (sort_esize)
                // Simply sort based on |esize|, largest first
                {
                    std::sort(mut_keys.begin(), mut_keys.end(),
                              [&pop](KTfwd::uint_t a, KTfwd::uint_t b) {
                                  return std::fabs(pop->mutations[a].s)
                                         > std::fabs(pop->mutations[b].s);
                              });
                }
            return mut_keys;
        }
        template <typename pop_t>
        void
        update_row_details(gsl_matrix *m, const typename pop_t::gamete_t &g,
                           const pop_t *pop,
                           const std::vector<KTfwd::uint_t> &mut_keys,
                           const size_t row)
        {
            for (auto &&k : g.smutations)
                {
                    if (pop->mcounts[k] < 2 * pop->N) // skip fixations!!!
                        {
                            if (!pop->mcounts[k])
                                throw std::runtime_error(
                                    "extinct mutation encountered: "
                                    + std::string(__FILE__) + ", "
                                    + std::to_string(__LINE__));
                            auto i = std::find(mut_keys.begin(),
                                               mut_keys.end(), k);
                            if (i == mut_keys.end())
                                throw std::runtime_error(
                                    "mutation key not found: "
                                    + std::string(__FILE__) + ", "
                                    + std::to_string(__LINE__) + ", "
                                    + "mcount = "
                                    + std::to_string(pop->mcounts[k]));
                            std::size_t col
                                = std::distance(mut_keys.begin(), i);
                            if (col + 1 >= m->size2)
                                throw std::runtime_error(
                                    "second dimension out of range: "
                                    + std::string(__FILE__) + ", "
                                    + std::to_string(__LINE__));
                            auto mp = gsl_matrix_ptr(m, row, col + 1);
                            *mp += 1.0; // update counts
                        }
                }
        }
        template <typename pop_t, typename diploid_t>
        void
        update_matrix_counts_details(
            gsl_matrix *m, const pop_t *pop,
            const std::vector<KTfwd::uint_t> &mut_keys, const diploid_t &dip,
            const size_t row)
        {
            update_row_details(m, pop->gametes[dip.first], pop, mut_keys, row);
            update_row_details(m, pop->gametes[dip.second], pop, mut_keys,
                               row);
        }
        template <typename pop_t>
        void
        update_matrix_counts(const pop_t *pop,
                             const std::vector<KTfwd::uint_t> &mut_keys,
                             gsl_matrix *rv)
        /*!
         * Fills rv with an 0,1,2 matrix of derived mutation counts.
         * Order of mutations is based on values in mut_keys, which are indexes
         * to pop->mutations/pop->mcounts.
         *
         * \note rv Should be zeroed out and have mut_keys.size()+1 columns.
         * Column 0 is set to 1.0
         */
        {
            // Fill the matrix
            std::size_t row = 0;
            for (const auto &dip : pop->diploids)
                {
                    gsl_matrix_set(rv, row, 0,
                                   1.0); // set column 0 to a value of 1.0
                    update_matrix_counts_details(rv, pop, mut_keys, dip, row);
                    row++;
                }
        }
        template <>
        inline void
        update_matrix_counts<multilocus_t>(
            const multilocus_t *pop,
            const std::vector<KTfwd::uint_t> &mut_keys, gsl_matrix *rv)
        /*!
          Return a 0,1,2 matrix of counts of causative alleles in each diploid.

          Specialization for fwdpy::multilocus_t
        */
        {
            // Fill the matrix
            std::size_t row = 0;
            for (const auto &dip : pop->diploids)
                {
                    gsl_matrix_set(rv, row, 0,
                                   1.0); // set column 0 to a value of 1.0
                    for (const auto &locus : dip)
                        {
                            update_matrix_counts_details(rv, pop, mut_keys,
                                                         locus, row);
                        }
                    row++;
                }
        }
    }
}

#endif
