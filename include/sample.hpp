/*!
  \file sample.hpp
  \brief Functions related to taking samples from populations.
*/
#ifndef __FWDPY_SAMPLE_HPP__
#define __FWDPY_SAMPLE_HPP__

#include "types.hpp"

namespace fwdpy
{
    struct popsample_details
    {
        std::vector<double> s, h, p;
        std::vector<unsigned> dcount, origin, generation, ftime, locus;
        std::vector<decltype(KTfwd::mutation_base::xtra)> label;
        popsample_details(
            std::vector<double> &&s_, std::vector<double> &&h_,
            std::vector<double> &&p_, std::vector<unsigned> &&dcount_,
            std::vector<unsigned> &&origin_,
            std::vector<unsigned> &&generation_,
            std::vector<unsigned> &&ftime_, std::vector<unsigned> &&locus_,
            std::vector<decltype(KTfwd::mutation_base::xtra)> &&label_)
            : s(std::move(s_)), h(std::move(h_)), p(std::move(p_)),
              dcount(std::move(dcount_)), origin(std::move(origin_)),
              generation(std::move(generation_)), ftime(std::move(ftime_)),
              locus(std::move(locus_)), label(std::move(label_))
        {
        }
        popsample_details()
            : s(std::vector<double>()), h(std::vector<double>()),
              p(std::vector<double>()), dcount(std::vector<unsigned>()),
              origin(std::vector<unsigned>()),
              generation(std::vector<unsigned>()),
              ftime(std::vector<unsigned>()), locus(std::vector<unsigned>()),
              label(decltype(label)())
        {
        }
    };

    inline std::pair<singlepop_t::mcont_t::const_iterator, bool>
    find_variant(const singlepop_t::mcont_t &mutations,
                 const singlepop_t::mcont_t &fixations,
                 const std::pair<double, std::string> &site)
    {
        auto mitr = std::find_if(fixations.begin(), fixations.end(),
                                 [&site](const singlepop_t::mutation_t &m) {
                                     return site.first == m.pos;
                                 });
        if (mitr != fixations.end())
            {
                return std::make_pair(mitr, true);
            }

        mitr = std::find_if(mutations.begin(), mutations.end(),
                            [&site](const singlepop_t::mutation_t &m) {
                                return site.first == m.pos;
                            });
        if (mitr == mutations.end()) // BAD
            {
                throw std::runtime_error("Variant at position "
                                         + std::to_string(site.first)
                                         + " could not be found");
            }
        return std::make_pair(mitr, false);
    }

    inline popsample_details
    get_sh_details(const std::vector<std::pair<double, std::string>> &sample,
                   const singlepop_t::mcont_t &mutations,
                   const std::vector<KTfwd::popgenmut> &fixations,
                   const std::vector<KTfwd::uint_t> &fixation_times,
                   const singlepop_t::mcount_t &mcounts, const size_t &twoN,
                   const unsigned &gen, const unsigned &locus_num)
    {
        std::vector<double> s, h, p;
        std::vector<unsigned> dcount, origin, generation, ftime, locus;
        std::vector<std::uint16_t> label;
        for (const auto &site : sample)
            {
                auto x = find_variant(mutations, fixations, site);
                s.push_back(x.first->s);
                h.push_back(x.first->h);
                origin.push_back(x.first->g);
                if (x.second)
                    {
                        p.push_back(1.0);
                        std::size_t idx = static_cast<std::size_t>(
                            std::distance(fixations.begin(), x.first));
                        ftime.push_back(fixation_times[idx] - fixations[idx].g
                                        + 1);
                    }
                else
                    {
                        p.push_back(double(mcounts[std::distance(
                                        mutations.begin(), x.first)])
                                    / (2.0 * double(twoN)));
                        ftime.push_back(std::numeric_limits<unsigned>::max());
                    }

                dcount.push_back( // count of derived allele in sample
                    std::count(site.second.begin(), site.second.end(), '1'));
                label.push_back(
                    x.first->xtra); // This is the 'label' assigned to a
                                    // mutation -- See Regions.pyx for
                                    // details.
            }
        if (s.empty())
            {
                // We add in meaningless data
                // to indicate no selected variants for this time period
                s.push_back(std::numeric_limits<double>::quiet_NaN());
                h = s;
                p = s;
                ftime.push_back(std::numeric_limits<unsigned>::max());
                dcount = ftime;
                origin = ftime;
                label.push_back(std::numeric_limits<std::uint16_t>::max());
                generation.push_back(gen);
                locus.push_back(locus_num);
            }
        else
            {
                generation.resize(s.size(), gen);
                locus.resize(s.size(), locus_num);
            }
        return popsample_details(std::move(s), std::move(h), std::move(p),
                                 std::move(dcount), std::move(origin),
                                 std::move(generation), std::move(ftime),
                                 std::move(locus), std::move(label));
    }
    /*!
      \brief Get detailed info about mutations in a sample.
      \note Definition in fwdpy/fwdpy/sample.cc
    */
    popsample_details
    get_sh(const std::vector<std::pair<double, std::string>> &ms_sample,
           const std::vector<KTfwd::popgenmut> &mutations,
           const std::vector<KTfwd::popgenmut> &fixations,
           const std::vector<KTfwd::uint_t> &fixation_times,
           const std::vector<KTfwd::uint_t> &mcounts, const unsigned &ttlN,
           const unsigned &generation, const unsigned &locus);
}

#endif
