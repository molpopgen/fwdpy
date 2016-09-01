#ifndef FWDPY_SAMPLE_N_HPP
#define FWDPY_SAMPLE_N_HPP

#include "sample.hpp"
#include "sampler_base.hpp"
#include "types.hpp"
#include <Sequence/SimData.hpp>
#include <algorithm>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/poptypes/tags.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
namespace fwdpy
{
    struct detailed_deme_sample
    {
        KTfwd::sep_sample_t genotypes;
        std::vector<std::pair<double, double>> sh;
        template <typename T1, typename T2>
        detailed_deme_sample(T1 &&t1, T2 &&t2)
            : genotypes(std::forward<T1>(t1)), sh(std::forward<T2>(t2))
        {
        }
    };

    class sample_n
        : public sampler_base // take a sample of size n from a population
                              /*
                                \brief A "sampler" that takes a sample of n gametes from a population
                                \ingroup samplers
                              */
    {
      public:
        using final_t
            = std::vector<std::pair<KTfwd::sep_sample_t, popsample_details>>;
        virtual void
        operator()(const singlepop_t *pop, const unsigned generation)
        {
            auto s = KTfwd::sample_separate(r.get(), *pop, nsam, removeFixed);
            if (nfile != NULL)
                {
                    gzFile gz = gzopen(nfile, "ab");
                    std::ostringstream o;
                    o << Sequence::SimData(s.first.begin(), s.first.end())
                      << '\n';
                    gzwrite(gz, o.str().c_str(), o.str().size());
                    gzclose(gz);
                }
            if (sfile != NULL)
                {
                    gzFile gz = gzopen(sfile, "ab");
                    std::ostringstream o;
                    o << Sequence::SimData(s.second.begin(), s.second.end())
                      << '\n';
                    gzwrite(gz, o.str().c_str(), o.str().size());
                    gzclose(gz);
                }
            std::vector<std::pair<double, double>> sh;
            for (const auto &i : s.second)
                {
                    auto itr = std::find_if(
                        pop->mutations.begin(), pop->mutations.end(),
                        [&i](const singlepop_t::mutation_t &m) noexcept {
                            return m.pos == i.first;
                        });
                    if (itr != pop->mutations.end())
                        {
                            sh.emplace_back(itr->s, itr->h);
                        }
                    else
                        {
                            itr = std::find_if(
                                pop->fixations.begin(), pop->fixations.end(),
                                [&i](const singlepop_t::mutation_t
                                         &m) noexcept {
                                    return m.pos == i.first;
                                });
                            if (itr != pop->fixations.end())
                                {
                                    sh.emplace_back(itr->s, itr->h);
                                }
                        }
                }
            auto details = get_sh_details(
                s.second, pop->mutations, pop->fixations, pop->mcounts,
                2 * pop->diploids.size(), generation);
            rv.emplace_back(std::move(s), std::move(details));
        }

        virtual void
        operator()(const multilocus_t *pop, const unsigned generation)
        {
            auto s = KTfwd::sample_separate(r.get(), *pop, nsam, removeFixed,
                                            locus_boundaries);
            if (nfile != NULL)
                {
                    gzFile gz = gzopen(nfile, "ab");
                    std::ostringstream o;
                    for (auto &&i : s)
                        {
                            o.str(std::string());
                            o << Sequence::SimData(i.first.begin(),
                                                   i.first.end())
                              << '\n';
                            gzwrite(gz, o.str().c_str(), o.str().size());
                        }
                    gzclose(gz);
                }
            if (sfile != NULL)
                {
                    gzFile gz = gzopen(sfile, "ab");
                    std::ostringstream o;
                    for (auto &&i : s)
                        {
                            o.str(std::string());
                            o << Sequence::SimData(i.second.begin(),
                                                   i.second.end())
                              << '\n';
                            gzwrite(gz, o.str().c_str(), o.str().size());
                        }
                    gzclose(gz);
                }
            std::vector<detailed_deme_sample> vds;
            for (unsigned i = 0; i < s.size(); ++i)
                {
                    std::vector<std::pair<double, double>> sh;
                    for (const auto &si : s[i].second)
                        {
                            auto itr = std::find_if(
                                pop->mutations.begin(), pop->mutations.end(),
                                [&si](const singlepop_t::mutation_t
                                          &m) noexcept {
                                    return m.pos == si.first;
                                });
                            sh.emplace_back(itr->s, itr->h);
                        }
                    auto details = get_sh_details(
                        s[i].second, pop->mutations, pop->fixations,
                        pop->mcounts, 2 * pop->diploids.size(), generation);
                    rv.emplace_back(std::move(s[i]), std::move(details));
                }
        }

        final_t
        final() const
        {
            return rv;
        }
        explicit sample_n(
            unsigned nsam_, const gsl_rng *r_, const bool rfixed = true,
            const char *neutral_file = NULL, const char *selected_file = NULL,
            const std::vector<std::pair<double, double>> &boundaries
            = std::vector<std::pair<double, double>>(),
            const bool append = true)
            : rv(final_t()), nsam(nsam_), r(GSLrng_t(gsl_rng_get(r_))),
              nfile(neutral_file), sfile(selected_file),
              locus_boundaries(boundaries), removeFixed(rfixed)
        /*!
          Note the implementation of this constructor!!

          By taking a gsl_rng * from outside, we are able to guarantee
          that this object is reproducibly seeded to the extent that
          this constructor is called in a reproducible order.
        */
        {
            if (!append)
                {
                    gzFile gz = gzopen(neutral_file, "wb");
                    if (gz == NULL)
                        {
                            throw std::runtime_error("could not open "
                                                     + std::string(nfile)
                                                     + " in 'wb' mode");
                        }
                    gzclose(gz);
                    gz = gzopen(selected_file, "wb");
                    if (gz == NULL)
                        {
                            throw std::runtime_error("could not open "
                                                     + std::string(sfile)
                                                     + " in 'wb' mode");
                        }
                    gzclose(gz);
                }
        }

      private:
        final_t rv;
        const unsigned nsam;
        GSLrng_t r;
        const char *nfile, *sfile;
        const std::vector<std::pair<double, double>> locus_boundaries;
        const bool removeFixed;
    };
}

#endif
