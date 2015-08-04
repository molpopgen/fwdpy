#ifndef __FWDPY_TYPES__
#define __FWDPY_TYPES__

#include <fwdpp/sugar.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
namespace fwdpy {
  using GSLrng_t = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;

  struct singlepop_t : public KTfwd::singlepop_serialized<KTfwd::popgenmut,KTfwd::mutation_writer,KTfwd::mutation_reader<KTfwd::popgenmut>>
  {
    using base = KTfwd::singlepop_serialized<KTfwd::popgenmut,KTfwd::mutation_writer,KTfwd::mutation_reader<KTfwd::popgenmut>>;
    unsigned generation;
    singlepop_t(const unsigned & N) : base(N),generation(0)
      {
      }
  };

  struct metapop_t : public KTfwd::metapop_serialized<KTfwd::popgenmut,KTfwd::mutation_writer,KTfwd::mutation_reader<KTfwd::popgenmut> >
  {
    using base = KTfwd::metapop_serialized<KTfwd::popgenmut,KTfwd::mutation_writer,KTfwd::mutation_reader<KTfwd::popgenmut> >;
    unsigned generation;
    metapop_t( std::initializer_list<unsigned> Ns ) : base(Ns), generation(0)
    {
    }
  };
  
}

#endif
