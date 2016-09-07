#ifndef FWDPY_HAPLOTYPE_MATRIX_HPP
#define FWDPY_HAPLOTYPE_MATRIX_HPP

#include "types.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fwdpy
{
    struct haplotype_matrix
    {
        //! Indexes of neutral markers in matrix
        std::vector<std::size_t> n;
        //! Indexes of selected markers in matrix
        std::vector<std::size_t> s;
        //! Positions of neutral markers
        std::vector<double> np;
        //! Frequencies of neutral markers
        std::vector<double> nf;
        //! Positions of selected markers
        std::vector<double> sp;
        //! Frequencies of selected markers
        std::vector<double> sf;
        //! Genetic value
        std::vector<double> G;
        //! Random value
        std::vector<double> E;
        //! fitness
        std::vector<double> w;
        //! Effect sizes of mutations
        std::vector<double> esizes;
        //! Dominances of mutations
        std::vector<double> h;
        std::size_t nrow, ncol_n, ncol_s;
    };

    template <typename mcont_t>
    std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
    get_mut_keys_common(const mcont_t &mutations,
                        const std::vector<KTfwd::uint_t> &mcounts,
                        const KTfwd::uint_t twoN,
                        const std::unordered_set<std::size_t> &n,
                        const std::unordered_set<std::size_t> &s,
                        haplotype_matrix &hm)
    {
        std::pair<std::vector<std::size_t>, std::vector<std::size_t>> rv;
        std::copy(n.begin(), n.end(), std::back_inserter(rv.first));
        std::copy(s.begin(), s.end(), std::back_inserter(rv.second));

        std::sort(rv.first.begin(), rv.first.end(),
                  [&mutations](const std::size_t i, const std::size_t j) {
                      return mutations[i].pos < mutations[j].pos;
                  });
        std::sort(rv.second.begin(), rv.second.end(),
                  [&mutations](const std::size_t i, const std::size_t j) {
                      return mutations[i].pos < mutations[j].pos;
                  });
        for (auto &ni : rv.first)
            {
                hm.np.push_back(mutations[ni].pos);
                hm.nf.push_back(double(mcounts[ni]) / double(twoN));
            }

        for (auto &si : rv.second)
            {
                hm.sp.push_back(mutations[si].pos);
                hm.esizes.push_back(mutations[si].s);
                hm.h.push_back(mutations[si].h);
                hm.sf.push_back(double(mcounts[si]) / double(twoN));
            }
        return rv;
    }

    template <typename dipvec_t, typename gcont_t, typename mcont_t>
    std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
    get_mut_keys(const dipvec_t &diploids, const gcont_t &gametes,
                 const mcont_t &mutations,
                 const std::vector<unsigned> &mcounts,
                 const std::vector<std::size_t> diploids_sample,
                 haplotype_matrix &hm)
    /*!
      Workhorse function for a single deme
    */
    {
        std::unordered_set<std::size_t> n, s;
        for (auto &dip : diploids_sample)
            {
                if (dip >= diploids.size())
                    throw std::out_of_range("diploid index out of range");
                for (auto &m : gametes[diploids[dip].first].mutations)
                    {
                        if (mcounts[m]
                            < 2 * diploids.size()) // avoid fixed variants
                            n.insert(m);
                    }
                for (auto &m : gametes[diploids[dip].second].mutations)
                    {
                        if (mcounts[m] < 2 * diploids.size())
                            n.insert(m);
                    }
                for (auto &m : gametes[diploids[dip].first].smutations)
                    {
                        if (mcounts[m] < 2 * diploids.size())
                            s.insert(m);
                    }
                for (auto &m : gametes[diploids[dip].second].smutations)
                    {
                        if (mcounts[m] < 2 * diploids.size())
                            s.insert(m);
                    }
            }
        auto rv = get_mut_keys_common(mutations, mcounts, 2 * diploids.size(),
                                      n, s, hm);
        hm.nrow = 2 * diploids_sample.size();
        hm.ncol_n = rv.first.size();
        hm.ncol_s = rv.second.size();
        return rv;
    }

    template <typename dipvec_t, typename gcont_t, typename mcont_t>
    std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
    get_mut_keys_mloc(const dipvec_t &diploids, const gcont_t &gametes,
                      const mcont_t &mutations,
                      const std::vector<unsigned> &mcounts,
                      const std::vector<std::size_t> diploids_sample,
                      haplotype_matrix &hm)
    {
        std::unordered_set<std::size_t> n, s;
        for (auto &dip : diploids_sample)
            {
                if (dip >= diploids.size())
                    throw std::out_of_range("diploid index out of range");
                for (auto &locus : diploids[dip])
                    {
                        for (auto &m : gametes[locus.first].mutations)
                            {
                                if (mcounts[m]
                                    < 2 * diploids.size()) // avoid fixed
                                                           // variants
                                    n.insert(m);
                            }
                        for (auto &m : gametes[locus.second].mutations)
                            {
                                if (mcounts[m] < 2 * diploids.size())
                                    n.insert(m);
                            }
                        for (auto &m : gametes[locus.first].smutations)
                            {
                                if (mcounts[m] < 2 * diploids.size())
                                    s.insert(m);
                            }
                        for (auto &m : gametes[locus.second].smutations)
                            {
                                if (mcounts[m] < 2 * diploids.size())
                                    s.insert(m);
                            }
                    }
            }
        auto rv = get_mut_keys_common(mutations, mcounts, 2 * diploids.size(),
                                      n, s, hm);
        hm.nrow = 2 * diploids_sample.size();
        hm.ncol_n = rv.first.size();
        hm.ncol_s = rv.second.size();
        return rv;
    }

    inline void
    update_matrix(const std::size_t key, const std::size_t row,
                  const std::vector<std::size_t> &keys,
                  std::vector<std::size_t> &indexes)
    {
        auto i = std::find(keys.begin(), keys.end(), key);
        if (i == keys.end())
            throw std::runtime_error("fatal error: mutation key not found: "
                                     + std::to_string(key));
        indexes.emplace_back(
            std::size_t(row * keys.size() + std::distance(keys.begin(), i)));
    }

    template <typename dipvec_t, typename gcont_t, typename mcont_t>
    haplotype_matrix
    make_haplotype_matrix_single_deme_details(
        const dipvec_t &diploids, const gcont_t &gametes,
        const mcont_t &mutations, const std::vector<unsigned> &mcounts,
        const std::vector<std::size_t> diploids_sample)
    {
        haplotype_matrix rv;
        // Step 1: get mutation keys/positions, no. rows/columns
        auto keys = get_mut_keys(diploids, gametes, mutations, mcounts,
                                 diploids_sample, rv);

        // Step 2: fill matrices
        std::size_t row = 0;
        for (auto &dip : diploids_sample)
            {
                rv.G.push_back(diploids[dip].g);
                rv.E.push_back(diploids[dip].e);
                rv.w.push_back(diploids[dip].w);
                for (auto &m : gametes[diploids[dip].first].mutations)
                    {
                        if (mcounts[m] < 2 * diploids.size())
                            update_matrix(m, row, keys.first, rv.n);
                    }
                for (auto &m : gametes[diploids[dip].first].smutations)
                    {
                        if (mcounts[m] < 2 * diploids.size())
                            update_matrix(m, row, keys.second, rv.s);
                    }
                ++row;
                for (auto &m : gametes[diploids[dip].second].mutations)
                    {
                        if (mcounts[m] < 2 * diploids.size())
                            update_matrix(m, row, keys.first, rv.n);
                    }
                for (auto &m : gametes[diploids[dip].second].smutations)
                    {
                        if (mcounts[m] < 2 * diploids.size())
                            update_matrix(m, row, keys.second, rv.s);
                    }
                ++row;
            }
        if (row != rv.nrow)
            throw std::runtime_error("fatal error: row != rv.nrow");
        return rv;
    }
    haplotype_matrix
    make_haplotype_matrix(const singlepop_t *pop,
                          const std::vector<std::size_t> &diploids)
    /*!
      Create a haplotype matrix from a single deme from a specified set of
      individuals
    */
    {
        return make_haplotype_matrix_single_deme_details(
            pop->diploids, pop->gametes, pop->mutations, pop->mcounts,
            diploids);
    }

    haplotype_matrix
    make_haplotype_matrix(const metapop_t *pop,
                          const std::vector<std::size_t> &diploids,
                          const std::size_t deme)
    {
        if (deme > pop->diploids.size())
            throw std::out_of_range("deme index out of range");
        return make_haplotype_matrix_single_deme_details(
            pop->diploids[deme], pop->gametes, pop->mutations, pop->mcounts,
            diploids);
    }
    haplotype_matrix
    make_haplotype_matrix(const multilocus_t *pop,
                          const std::vector<std::size_t> &diploids)
    /*!
      Create a haplotype matrix from a single deme (multi-locus) from a
      specified set of individuals
    */
    {
        haplotype_matrix rv;
        auto keys
            = get_mut_keys_mloc(pop->diploids, pop->gametes, pop->mutations,
                                pop->mcounts, diploids, rv);
        std::size_t row = 0;
        for (auto &dip : diploids)
            {
                rv.G.push_back(pop->diploids[dip][0].g);
                rv.E.push_back(pop->diploids[dip][0].e);
                rv.w.push_back(pop->diploids[dip][0].w);
                for (auto &locus : pop->diploids[dip])
                    {
                        for (auto &m : pop->gametes[locus.first].mutations)
                            {
                                if (pop->mcounts[m] < 2 * pop->diploids.size())
                                    update_matrix(m, row, keys.first, rv.n);
                            }
                        for (auto &m : pop->gametes[locus.first].smutations)
                            {
                                if (pop->mcounts[m] < 2 * pop->diploids.size())
                                    update_matrix(m, row, keys.second, rv.s);
                            }
                        ++row;
                        for (auto &m : pop->gametes[locus.second].mutations)
                            {
                                if (pop->mcounts[m] < 2 * pop->diploids.size())
                                    update_matrix(m, row, keys.first, rv.n);
                            }
                        for (auto &m : pop->gametes[locus.second].smutations)
                            {
                                if (pop->mcounts[m] < 2 * pop->diploids.size())
                                    update_matrix(m, row, keys.second, rv.s);
                            }
                        ++row;
                    }
            }
        return rv;
    }

    inline void
    update_genotype_matrix(const std::multimap<std::size_t, std::size_t> &m,
                           const std::size_t hap_row,
                           const std::size_t geno_row, const std::size_t ncol,
                           std::vector<std::size_t> &Aa,
                           std::vector<std::size_t> &aa)
    {
        // find mutationgs in this rows
        auto r1 = m.equal_range(hap_row);
        auto r2 = m.equal_range(hap_row + 1);
        std::unordered_map<std::size_t, std::size_t> counts;
        for (auto j = r1.first; j != r1.second; ++j)
            {
                counts[j->second]++;
            }
        for (auto j = r2.first; j != r2.second; ++j)
            {
                counts[j->second]++;
            }
        for (const auto &j : counts)
            {
                if (j.second == 1)
                    {
                        Aa.push_back(geno_row * ncol + j.first);
                    }
                else if (j.second == 2)
                    {
                        aa.push_back(geno_row * ncol + j.first);
                    }
            }
    }

    std::map<std::string, std::vector<std::size_t>>
    make_genotype_matrix(const haplotype_matrix &hm)
    //! Take haplotype matrix and return info needed to make a genotype matrix
    {
        std::vector<std::size_t> nAa, naa, sAa, saa;

        // record row/column positions of each element in hm
        std::multimap<std::size_t, std::size_t> n, s;
        auto position = [&hm](const std::size_t i, const double ncol) {
            double t = double(i + 1) / ncol;
            double row = std::floor(t);
            double col = (t - row) * ncol - 1.0;
            return std::make_pair(std::size_t(row), std::size_t(col));
        };
        for (const auto &i : hm.n)
            n.insert(position(i, double(hm.ncol_n)));
        for (const auto &i : hm.s)
            s.insert(position(i, double(hm.ncol_s)));

        unsigned nrow = 0;
        for (std::size_t i = 0; i < hm.nrow; i += 2)
            {
                update_genotype_matrix(n, i, nrow, hm.ncol_n, nAa, naa);
                update_genotype_matrix(s, i, nrow, hm.ncol_s, sAa, saa);
                ++nrow;
            }

        return std::map<std::string, std::vector<std::size_t>>{
            { "nAa", std::move(nAa) },
            { "naa", std::move(naa) },
            { "sAa", std::move(sAa) },
            { "saa", std::move(saa) }
        };
    }
}

#endif
