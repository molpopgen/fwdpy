#ifndef FWDPY_SAMPLE_N_HPP
#define FWDPY_SAMPLE_N_HPP

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

      private:
        final_t rv;
        const unsigned nsam;
        GSLrng_t r;
        const std::string nfile, sfile;
        const std::vector<std::pair<double, double>> locus_boundaries;
        const bool removeFixed, recordSamples, recordDetails;

        void
        remove_redundant_selected_fixations(KTfwd::sep_sample_t &sample)
        /* In fwdpy, we have a slight issue:
         * 1. Quant-trait sims do copy selected fixations
         *    into fixations/fixation times, while also
         *    leaving them in the pop
         * 2. Thus, when removeFixed = false, we can end up w/two
         *    copies of such mutations.  Essentially, we are using
         *    fwdpp as it was not intended to be used :).
         * 3. This function cleans up the mess we've made
         */
        {
            sample.second.erase(
                std::unique(sample.second.begin(), sample.second.end()),
                sample.second.end());
        }

        void
        write_sample_details(gzFile gz, const Sequence::SimData &d)
        {
            std::stringstream o;
            o << d << '\n';
            gzwrite(gz, o.str().c_str(), o.str().size());
        }

        void
        write_sample(const std::string &fn, const KTfwd::sample_t &s)
        {
            gzFile gz = gzopen(fn.c_str(), "ab");
            write_sample_details(gz, Sequence::SimData(s.begin(), s.end()));
            gzclose(gz);
        }

        void
        write_mloc_sample(const std::string &fn,
                          const std::vector<KTfwd::sep_sample_t> &s,
                          const bool neutral)
        {
            gzFile gz = gzopen(fn.c_str(), "ab");
            for (auto &&si : s)
                {
                    if (neutral)
                        {
                            write_sample_details(
                                gz, Sequence::SimData(si.first.begin(),
                                                      si.first.end()));
                        }
                    else
                        {
                            write_sample_details(
                                gz, Sequence::SimData(si.second.begin(),
                                                      si.second.end()));
                        }
                }
            gzclose(gz);
        }

      public:
        virtual void
        operator()(const singlepop_t *pop, const unsigned generation)
        {
            auto s = KTfwd::sample_separate(r.get(), *pop, nsam, removeFixed);
            remove_redundant_selected_fixations(s);
            if (!nfile.empty())
                {
                    write_sample(nfile, s.first);
                }
            if (!sfile.empty())
                {
                    write_sample(sfile, s.second);
                }
            if (recordDetails)
                {
                    auto details = get_sh_details(
                        s.second, pop->mutations, pop->fixations,
                        pop->fixation_times, pop->mcounts,
                        pop->diploids.size(), generation, 0);
                    if (recordSamples)
                        {
                            rv.emplace_back(std::move(s), std::move(details));
                        }
                    else
                        {
                            rv.emplace_back(final_t::value_type::first_type(),
                                            std::move(details));
                        }
                }
            else if (recordSamples)
                {
                    rv.emplace_back(std::move(s),
                                    final_t::value_type::second_type());
                }
        }

        virtual void
        operator()(const multilocus_t *pop, const unsigned generation)
        {
            auto s = KTfwd::sample_separate(r.get(), *pop, nsam, removeFixed,
                                            locus_boundaries);
            for (auto &si : s)
                {
                    remove_redundant_selected_fixations(si);
                }
            if (!nfile.empty())
                {
                    write_mloc_sample(nfile, s, true);
                }
            if (!sfile.empty())
                {
                    write_mloc_sample(sfile, s, false);
                }
            for (unsigned i = 0; i < s.size(); ++i)
                {
                    if (recordDetails)
                        {
                            auto details = get_sh_details(
                                s[i].second, pop->mutations, pop->fixations,
                                pop->fixation_times, pop->mcounts,
                                pop->diploids.size(), generation, i);
                            if (recordSamples)
                                {
                                    rv.emplace_back(std::move(s[i]),
                                                    std::move(details));
                                }
                            else
                                {
                                    rv.emplace_back(
                                        final_t::value_type::
                                            first_type(),
                                        std::move(details));
                                }
                        }
                    else if (recordSamples)
                        {
                            rv.emplace_back(std::move(s[i]),
                                            final_t::value_type::
                                                second_type());
                        }
                }
        }
        final_t
        final() const
        {
            return rv;
        }
        explicit sample_n(
            unsigned nsam_, const gsl_rng *r_, const std::string &neutral_file,
            const std::string &selected_file, const bool rfixed = true,
            const bool rec_samples = true, const bool rec_sh = true,
            const std::vector<std::pair<double, double>> &boundaries
            = std::vector<std::pair<double, double>>(),
            const bool append = true)
            : rv(final_t()), nsam(nsam_), r(GSLrng_t(gsl_rng_get(r_))),
              nfile(neutral_file), sfile(selected_file),
              locus_boundaries(boundaries), removeFixed(rfixed),
              recordSamples(rec_samples), recordDetails(rec_sh)
        /*!
          Note the implementation of this constructor!!

          By taking a gsl_rng * from outside, we are able to guarantee
          that this object is reproducibly seeded to the extent that
          this constructor is called in a reproducible order.
        */
        {
            if (!append)
                {
                    if (!neutral_file.empty())
                        {
                            gzFile gz = gzopen(neutral_file.c_str(), "wb");
                            if (gz == NULL)
                                {
                                    throw std::runtime_error(
                                        "could not open " + std::string(nfile)
                                        + " in 'wb' mode");
                                }
                            gzclose(gz);
                        }
                    if (!selected_file.empty())
                        {
                            gzFile gz = gzopen(selected_file.c_str(), "wb");
                            if (gz == NULL)
                                {
                                    throw std::runtime_error(
                                        "could not open " + std::string(sfile)
                                        + " in 'wb' mode");
                                }
                            gzclose(gz);
                        }
                }
        }
    };
}

#endif
