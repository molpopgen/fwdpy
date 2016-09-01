/*!
  \file sample.hpp
  \brief Functions related to taking samples from populations.
*/
#ifndef __FWDPY_SAMPLE_HPP__
#define __FWDPY_SAMPLE_HPP__

#include "types.hpp"

namespace fwdpy
{
    /*!
      \brief Get detailed info about mutations in a sample.
      \note Definition in fwdpy/fwdpy/sample.cc
    */
    void get_sh(const std::vector<std::pair<double, std::string>> &ms_sample,
                const std::vector<KTfwd::popgenmut> &mutations,
                const std::vector<KTfwd::popgenmut> &fixations,
                const std::vector<KTfwd::uint_t> &mcounts,
                const unsigned &ttlN, const unsigned &generation,
                std::vector<double> *s, std::vector<double> *h,
                std::vector<double> *p, std::vector<double> *a,
				std::vector<unsigned> *c,
                std::vector<decltype(KTfwd::mutation_base::xtra)> *l);
}

#endif
