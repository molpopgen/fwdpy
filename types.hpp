#ifndef __FWDPY_TYPES__
#define __FWDPY_TYPES__

#include <fwdpp/sugar.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
namespace fwdpy {
  using GSLrng_t = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;
  using singlepop_t = KTfwd::singlepop_serialized<KTfwd::popgenmut,KTfwd::mutation_writer,KTfwd::mutation_reader<KTfwd::popgenmut>>;
}

#endif
