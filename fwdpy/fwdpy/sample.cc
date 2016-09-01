#include "sample.hpp"
#include <algorithm>
#include <fwdpp/diploid.hh>
#include <limits>
#include <stdexcept>

namespace fwdpy
{
    popsample_details
    get_sh_details(const std::vector<std::pair<double, std::string>> &sample,
                   const singlepop_t::mcont_t &mutations,
                   const std::vector<KTfwd::popgenmut> &fixations,
                   const singlepop_t::mcount_t &mcounts, const size_t &twoN,
                   const unsigned &gen)
    {
        std::vector<double> s, h, p;
        std::vector<unsigned> age, dcount, generation;
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
                s.push_back(std::numeric_limits<double>::quiet_NaN());
                h = s;
                p = s;
                age.push_back(std::numeric_limits<unsigned>::max());
                dcount = age;
                label.push_back(std::numeric_limits<std::uint16_t>::max());
                generation.push_back(gen);
            }
        else
            {
                generation.resize(s.size(), gen);
            }
        return popsample_details(std::move(s), std::move(h), std::move(p),
                                 std::move(age), std::move(dcount),
                                 std::move(generation), std::move(label));
    }

    popsample_details
    get_sh(const std::vector<std::pair<double, std::string>> &samples,
           const std::vector<KTfwd::popgenmut> &mutations,
           const std::vector<KTfwd::popgenmut> &fixations,
           const singlepop_t::mcount_t &mcounts, const KTfwd::uint_t &ttlN,
           const unsigned generation)
    {
        return get_sh_details(samples, mutations, fixations, mcounts, ttlN,
                              generation);
    }
}
