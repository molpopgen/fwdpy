#ifndef FWDPY_NO_SAMPLING_HPP
#define FWDPY_NO_SAMPLING_HPP

#include "sampler_base.hpp"
/*!
  \defgroup samplers Function objects to record info from simulations at
  regular intervals.
*/

namespace fwdpy
{
    struct no_sampling : public sampler_base
    /*!
      \brief A "sampler" that does nothing
      \ingroup samplers
    */
    {
        using final_t = void;
        virtual void operator()(const singlepop_t *, const unsigned) final{};
        virtual void operator()(const multilocus_t *, const unsigned) final{};
        virtual void operator()(const metapop_t *, const unsigned) final{};

        final_t
        final() const
        {
            return;
        }
    };
}

#endif
