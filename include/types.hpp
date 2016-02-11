#ifndef __FWDPY_TYPES__
#define __FWDPY_TYPES__

#include <map>
#include <tuple>
#include <memory>
#include <vector>
#include <array>
#include <fwdpp/tags/diploid_tags.hpp>
#include <fwdpp/sugar.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <gsl/gsl_statistics_double.h>

/*
  NAMESPACE POLLUTION!!!!!

  Types defined here are in the global namespace.

  This is bad, BUT, it allows for auto-conversion of
  struct to dict via Cython.

  Currently, Cython will fail to compile auto-conversion
  code for structs declared inside a C++ namespace.
*/

struct qtrait_stats_cython
{
  std::string stat;
  double value;
  unsigned generation;
  qtrait_stats_cython(std::string _stat,
		      double _v, unsigned _g) : stat(std::move(_stat)),value(_v),generation(_g)
  {
  }
};

namespace fwdpy {
  using GSLrng_t = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;

  using mcont_t = std::vector<KTfwd::popgenmut>;
  using gamete_t = KTfwd::gamete;
  using gcont_t = std::vector<gamete_t>;

  struct diploid_t : public KTfwd::tags::custom_diploid_t
  {
    using first_type = std::size_t;
    using second_type = std::size_t;
    first_type first;
    second_type second;
    double g,e,w;
    diploid_t() noexcept : first(first_type()),second(second_type()),g(0.),e(0.),w(0.) {}
    diploid_t(first_type g1, first_type g2) noexcept : first(g1),second(g2),g(0.),e(0.),w(0.) {}
  };

  using dipvector_t = std::vector<diploid_t>;

  struct diploid_writer
  {
    using result_type = void;
    template<typename diploid_t>
    inline result_type operator()( const diploid_t & dip, std::ostream & o ) const
    {
      o.write( reinterpret_cast<const char *>(&dip.g),sizeof(double));
      o.write( reinterpret_cast<const char *>(&dip.e),sizeof(double));
      o.write( reinterpret_cast<const char *>(&dip.w),sizeof(double));
    }
  };

  struct diploid_reader
  {
    using result_type = void;
    template<typename diploid_t>
    inline result_type operator()( diploid_t & dip, std::istream & i ) const
    {
      i.read( reinterpret_cast<char *>(&dip.g),sizeof(double));
      i.read( reinterpret_cast<char *>(&dip.e),sizeof(double));
      i.read( reinterpret_cast<char *>(&dip.w),sizeof(double));
    }
  };

  enum class traj_key_values : std::size_t
  {
    deme,origin,pos,esize
      };

  using trajectories_key_t = std::tuple<unsigned,unsigned,double,double>;
  using trajectories_t = std::map< trajectories_key_t , std::vector<double> >;

  enum class qtrait_stat_names : std::size_t { GEN,VG,VE,PLF,LE,MAXEXP,EBAR,WBAR };

  using qtrait_stats_t = std::vector<std::array<double,8>>;

  struct singlepop_t :  public KTfwd::singlepop<KTfwd::popgenmut,diploid_t>
  {
    using base = KTfwd::singlepop<KTfwd::popgenmut,diploid_t>;
    unsigned generation;
    singlepop_t(const unsigned & N) : base(N),generation(0)
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
  };

  std::vector<qtrait_stats_cython> convert_qtrait_stats( const singlepop_t * pop );

  struct metapop_t : public KTfwd::metapop<KTfwd::popgenmut,diploid_t>
  {
    using base = KTfwd::metapop<KTfwd::popgenmut,diploid_t>;
    unsigned generation;
    metapop_t( const std::vector<unsigned> &Ns ) : base(&Ns[0],Ns.size()), generation(0)
    {
    }
    unsigned gen() const
    {
      return generation;
    }
    std::vector<unsigned> popsizes() const
    {
      return Ns;
    }
    int sane() const
    {
      for(unsigned i=0;i<diploids.size();++i)
	{
	  if(diploids[i].size()!=Ns[i]) return 0;
	}
      return 1;
    }
    int size() const
    {
      return int(diploids.size());
    }
  };

  //Types based on KTfwd::generalmut_vec
  using gamete_gm_vec_t = KTfwd::gamete;
  using glist_gm_vec_t = std::vector<gamete_gm_vec_t>;
  using mlist_gm_vec_t = std::vector<KTfwd::generalmut_vec>;

  struct diploid_gm_vec_t : public KTfwd::tags::custom_diploid_t
  {
    using first_type = std::size_t;
    using second_type = std::size_t;
    first_type first;
    second_type second;
    double g,e,w;
    diploid_gm_vec_t() : first(first_type()),second(second_type()),g(0.),e(0.),w(0.) {}
    diploid_gm_vec_t(first_type g1, first_type g2) : first(g1),second(g2),g(0.),e(0.),w(0.) {}
  };

  using dipvector_gm_vec_t = std::vector<diploid_gm_vec_t>;

  struct singlepop_gm_vec_t :  public KTfwd::singlepop<KTfwd::generalmut_vec,diploid_t>
  {
    using base = KTfwd::singlepop<KTfwd::generalmut_vec,diploid_t>;
    unsigned generation;
    singlepop_gm_vec_t(const unsigned & N) : base(N),generation(0)
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
  };
}

#endif
