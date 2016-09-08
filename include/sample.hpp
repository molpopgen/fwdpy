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
        std::vector<unsigned> age, dcount, generation, locus;
        std::vector<decltype(KTfwd::mutation_base::xtra)> label;
        popsample_details(
            std::vector<double> &&s_, std::vector<double> &&h_,
            std::vector<double> &&p_, std::vector<unsigned> &&age_,
            std::vector<unsigned> &&dcount_,
            std::vector<unsigned> &&generation_,
            std::vector<unsigned> &&locus_,
            std::vector<decltype(KTfwd::mutation_base::xtra)> &&label_)
            : s(std::move(s_)), h(std::move(h_)), p(std::move(p_)),
              age(std::move(age_)), dcount(std::move(dcount_)),
              generation(std::move(generation_)), locus(std::move(locus_)),
              label(std::move(label_))
        {
        }
    };

    inline popsample_details
    get_sh_details(const std::vector<std::pair<double, std::string>> &sample,
                   const singlepop_t::mcont_t &mutations,
                   const std::vector<KTfwd::popgenmut> &fixations,
                   const singlepop_t::mcount_t &mcounts, const size_t &twoN,
                   const unsigned &gen, const unsigned &locus_num)
    {
        std::vector<double> s, h, p;
        std::vector<unsigned> age, dcount, generation, locus;
        std::vector<std::uint16_t> label;
        for (const auto &site : sample)
            {
                // try to find mutation amongst
                // segregating variants
                auto mitr
                    = std::find_if(mutations.begin(), mutations.end(),
                                   [&site](const singlepop_t::mutation_t &m) {
                                       return site.first == m.pos;
                                   });
                // if that failed, try fixations
                auto mitr2
                    = std::find_if(fixations.begin(), fixations.end(),
                                   [&site](const singlepop_t::mutation_t &m) {
                                       return site.first == m.pos;
                                   });
                if (mitr == mutations.end())
                    {
                        if (mitr2 == fixations.end())
                            {
                                throw std::runtime_error(
                                    "mutation in sample not found in either "
                                    "mutations nor fixations");
                            }
                        else
                            {
                                mitr = mitr2;
                            }
                    }
                s.push_back(mitr->s);
                h.push_back(mitr->h);
                p.push_back(
                    double(mcounts[std::distance(mutations.begin(), mitr)])
                    / (2.0 * double(twoN)));
                age.push_back((gen - mitr->g)); // mutation age--this is
                                                // correct b/c of def'n of
                                                // 'gen' in the pop
                                                // objects!
                dcount.push_back( // count of derived allele in sample
                    std::count(site.second.begin(), site.second.end(), '1'));
                label.push_back(
                    mitr->xtra); // This is the 'label' assigned to a
                                 // mutation -- See Regions.pyx for
                                 // details.
            }
        if (s.empty())
            {
				//We add in meaningless data
				//to indicate no selected variants for this time period
                s.push_back(std::numeric_limits<double>::quiet_NaN());
                h = s;
                p = s;
                age.push_back(std::numeric_limits<unsigned>::max());
                dcount = age;
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
                                 std::move(age), std::move(dcount),
                                 std::move(generation), std::move(locus),
                                 std::move(label));
    }
    /*!
      \brief Get detailed info about mutations in a sample.
      \note Definition in fwdpy/fwdpy/sample.cc
    */
    popsample_details
    get_sh(const std::vector<std::pair<double, std::string>> &ms_sample,
           const std::vector<KTfwd::popgenmut> &mutations,
           const std::vector<KTfwd::popgenmut> &fixations,
           const std::vector<KTfwd::uint_t> &mcounts, const unsigned &ttlN,
           const unsigned &generation, const unsigned &locus);
}

#endif
