#ifndef FWDPY_NO_SAMPLING_HPP
#define FWDPY_NO_SAMPLING_HPP

namespace fwdpy
{
  struct no_sampling
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
