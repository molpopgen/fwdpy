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
        std::vector<unsigned> age, dcount, generation;
        std::vector<decltype(KTfwd::mutation_base::xtra)> labels;
        popsample_details(
            std::vector<double> &&s_, std::vector<double> &&h_,
            std::vector<double> &&p_, std::vector<unsigned> &&age_,
            std::vector<unsigned> &&dcount_,
            std::vector<unsigned> &&generation_,
            std::vector<decltype(KTfwd::mutation_base::xtra)> &&labels_)
            : s(std::move(s_)), h(std::move(h_)), p(std::move(p_)),
              age(std::move(age_)), dcount(std::move(dcount_)),
              generation(std::move(generation_)), labels(std::move(labels_))
        {
        }
    };
    /*!
      \brief Get detailed info about mutations in a sample.
      \note Definition in fwdpy/fwdpy/sample.cc
    */
    popsample_details
    get_sh(const std::vector<std::pair<double, std::string>> &ms_sample,
           const std::vector<KTfwd::popgenmut> &mutations,
           const std::vector<KTfwd::popgenmut> &fixations,
           const std::vector<KTfwd::uint_t> &mcounts, const unsigned &ttlN,
           const unsigned &generation);
}

#endif
