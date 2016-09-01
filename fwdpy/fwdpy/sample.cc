#include "sample.hpp"
#include <algorithm>
#include <fwdpp/diploid.hh>
#include <limits>
#include <stdexcept>

namespace fwdpy
{
    popsample_details
    get_sh(const std::vector<std::pair<double, std::string>> &samples,
           const std::vector<KTfwd::popgenmut> &mutations,
           const std::vector<KTfwd::popgenmut> &fixations,
           const singlepop_t::mcount_t &mcounts, const KTfwd::uint_t &ttlN,
           const unsigned &generation, const unsigned &locusID)
    {
        return get_sh_details(samples, mutations, fixations, mcounts, ttlN,
                              generation, locusID);
    }
}
