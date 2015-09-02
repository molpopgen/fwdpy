#ifndef __FWDPY_TYPES__
#define __FWDPY_TYPES__

#include <fwdpp/sugar.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <memory>
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

  struct popvector
  {
    using poptype = singlepop_t;
    using ptr_t = std::unique_ptr<poptype>;
    std::vector<ptr_t> pops;
    popvector(const unsigned npops,
	      const unsigned N) : pops(std::vector<ptr_t>())
    {
      for(unsigned i=0;i<npops;++i)
	{
	  pops.emplace_back(ptr_t(new poptype(N)));
	}
    }
    const singlepop_t * get( size_t i ) const
    {
      return pops[i].get();
    }
    int size() const
    {
      return int(pops.size());
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
