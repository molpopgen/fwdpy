#ifndef FWDPY_NO_SAMPLING_HPP
#define FWDPY_NO_SAMPLING_HPP

/*!
  \defgroup samplers Function objects to record info from simulations at regular intervals.
*/

namespace fwdpy
{
  struct no_sampling
  /*!
    \brief A "sampler" that does nothing
    \ingroup samplers
  */
  {
    using final_t = void;
    inline void operator()(const singlepop_t * ,gsl_rng * ,
			   const unsigned )
    {
      return;
    }

    final_t final() const
    {
      return;
    }
  };
}

#endif
