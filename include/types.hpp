#ifndef __FWDPY_TYPES__
#define __FWDPY_TYPES__

#include <fwdpp/tags/diploid_tags.hpp>
#include <fwdpp/sugar.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <memory>
namespace fwdpy {
  using GSLrng_t = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;

  using gamete_t = KTfwd::gamete_base<KTfwd::popgenmut>;
  using glist_t = std::list<gamete_t,std::allocator<gamete_t>>;
  
  struct diploid_t : public KTfwd::tags::custom_diploid_t
  {
    using first_type = glist_t::iterator;
    using second_type = glist_t::iterator;
    first_type first;
    second_type second;
    double g,e,w;
    diploid_t() : first(first_type()),second(second_type()),g(0.),e(0.),w(0.) {}
    diploid_t(first_type g1, first_type g2) : first(g1),second(g2),g(0.),e(0.),w(0.) {}
  };

  struct diploid_writer
  {
    using result_type = void;
    template<typename iterator>
    inline result_type operator()( iterator dip, std::ostream & o ) const
    {
      o.write( reinterpret_cast<const char *>(&dip->g),sizeof(double));
      o.write( reinterpret_cast<const char *>(&dip->e),sizeof(double));
      o.write( reinterpret_cast<const char *>(&dip->w),sizeof(double));
    }
  };

  struct diploid_reader
  {
    using result_type = void;
    template<typename iterator>
    inline result_type operator()( iterator dip, std::istream & i ) const
    {
      i.read( reinterpret_cast<char *>(&dip->g),sizeof(double));
      i.read( reinterpret_cast<char *>(&dip->e),sizeof(double));
      i.read( reinterpret_cast<char *>(&dip->w),sizeof(double));
    }
  };
  
  struct singlepop_t :  public KTfwd::singlepop_serialized<KTfwd::popgenmut,
							   KTfwd::mutation_writer,
							   KTfwd::mutation_reader<KTfwd::popgenmut>,
							   diploid_t,
							   diploid_writer,
							   diploid_reader>
  {
    using base = KTfwd::singlepop_serialized<KTfwd::popgenmut,
					     KTfwd::mutation_writer,
					     KTfwd::mutation_reader<KTfwd::popgenmut>,
					     diploid_t,
					     diploid_writer,
					     diploid_reader>;
    using trajtype = std::map< std::pair<unsigned,std::pair<double,double> >, std::vector<double> >;
    unsigned generation;
    trajtype trajectories;
    singlepop_t(const unsigned & N) : base(N),generation(0),
				      trajectories(trajtype())
    {
    }
    unsigned gen() const
    {
      return generation;
    }
    unsigned popsize() const
    {
      return N;
    }
    int sane() const
    {
      return int(N == diploids.size());
    }
    //defind in fwdpy/src/poptypes.cpp
    void updateTraj();
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
