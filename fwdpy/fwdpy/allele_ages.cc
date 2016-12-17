#include "allele_ages.hpp"
#include <algorithm>
#include <algorithm>
#include <iostream>
#include <stdexcept>

using namespace std;

namespace fwdpy
{
    std::vector<std::pair<unsigned, double>>::const_iterator
    find_max_element(const std::vector<std::pair<unsigned, double>> &traj)
    {
        // only do this check if possible
        // to fail
        using element_t = std::pair<unsigned, double>;
        return max_element(traj.begin(), traj.end(),
                           [](const element_t &a, const element_t &b) {
                               return a.second <= b.second;
                           });
    }

    vector<allele_age_data_t>
    allele_ages_details(const selected_mut_tracker::final_t &trajectories,
                        const double minfreq, const unsigned minsojourn)
    {
        if (minfreq < 0.0)
            throw runtime_error("minfreq must be >= 0.0");
        vector<allele_age_data_t> rv;
        for (const auto &t : trajectories)
            {
                if (t.second.empty())
                    {
                        throw runtime_error("frequency vector empty");
                    }
                if (t.second.size() >= minsojourn)
                    {
                        auto mfi = find_max_element(t.second);
                        if (mfi->second
                            >= minfreq) // it hit the right minimum frequency
                            {
                                rv.emplace_back(t.first.esize, mfi->second,
                                                t.second.back().second,
                                                t.first.origin,
                                                t.second.size());
                            }
                    }
            }
        return rv;
    }

    selected_mut_tracker::final_t
    merge_trajectories_details(const selected_mut_tracker::final_t &traj1,
                               const selected_mut_tracker::final_t &traj2)
    {
        selected_mut_tracker::final_t rv(traj1.begin(), traj1.end());
        for (auto &&t : traj2)
            {
                auto x = std::find_if(
                    rv.begin(), rv.end(),
                    [&t](const selected_mut_tracker::final_t::value_type &xi) {
                        return xi.first == t.first;
                    });
                if (x == rv.end())
                    {
                        rv.push_back(t);
                    }
                else
                    {
                        x->second.insert(x->second.end(), t.second.begin(),
                                         t.second.end());
                    }
            }
        return rv;
    }
}
