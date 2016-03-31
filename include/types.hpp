/*! 
  \file types.hpp

  \brief Wrappers around fwdpp types for the fwdpy Python package.
*/

#ifndef __FWDPY_TYPES__
#define __FWDPY_TYPES__

#include <map>
#include <tuple>
#include <memory>
#include <vector>
#include <fwdpp/tags/diploid_tags.hpp>
#include <fwdpp/sugar.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <gsl/gsl_statistics_double.h>

namespace fwdpy {

  /*!
    Random number generator.

    This is a std::unique_ptr wrapper to a gsl_rng * initialized
    as a Mersenne twister type (gsl_rng_mt19937).
  */
  using GSLrng_t = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;

  //! Typedef for mutation container
  using mcont_t = std::vector<KTfwd::popgenmut>;
  //! Typedef for gamete type
  using gamete_t = KTfwd::gamete;
  //! Typedef for gamete container
  using gcont_t = std::vector<gamete_t>;

  struct diploid_t : public KTfwd::tags::custom_diploid_t
  /*!
    \brief Custom diploid type.
  */
  {
    using first_type = std::size_t;
    using second_type = std::size_t;
    //! First gamete.  A gamete is vector<size_t> where the elements are indexes to a population's gamete container
    first_type first;
    //! Second gamete. A gamete is vector<size_t> where the elements are indexes to a population's gamete container
    second_type second;
    //! Genetic component of trait value.  This is not necessarily written to by a simulation.
    double g;
    //! Random component of trait value.  This is not necessarily written to by a simulation.
    double e;
    //! Fitness.  This is not necessarily written to by a simulation.
    double w;
    //! Constructor
    diploid_t() noexcept : first(first_type()),second(second_type()),g(0.),e(0.),w(0.) {}
    //! Construct from two indexes to gametes
    diploid_t(first_type g1, first_type g2) noexcept : first(g1),second(g2),g(0.),e(0.),w(0.) {}
  };

  //! Typedef for container of diploids
  using dipvector_t = std::vector<diploid_t>;

  //! Allows serialization of diploids.
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

  //! Allows de-serialization of diploids.
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

  struct singlepop_t :  public KTfwd::singlepop<KTfwd::popgenmut,diploid_t>
  /*!
    \brief Single-deme object where mutations have single effect size and dominance.

    This is the C++ representation of a single-deme simulation where a KTfwd::popgenmut
    is the mutation type, and a custom diploid (fwdpy::diploid_t) is the diploid type.

    This type inherits from fwdpp's type KTfwd::singlepop, and the main documentation
    is found in the fwdpp reference manual. 

    Here is a brief overview of the most important members of the base class.

    The format is name: type, comment

    mutations: vector<popgenmut>, contains a mix of extinct, segregating, and possibly fixed variants, depending on the simulation type.
    mcounts: vector<unsigned>, is the same length as mutations, and records the count (frequency, as an integer) of each mutation.
    gametes: vector<gamete>, contains a mix of extinct and extant gametes
    diploids: vector<diploid>
    fixations: vector<popgenmut>, is filled by simulations that "prine" fixations when they occur
    fixation_times: vector<unsigned>, is filled by simulations that "prine" fixations when they occur

    \note Internally, fwdppy uses extinct elements in containers for "object recycling."

    Further:

    1. For each diploid, first and second are indexes into gametes
    2. Within each gamete, mutations and smutations contain indexes to "neutral" and "selected" mutations, resepectively, in mutations
  */
  {
    //! Typedef for base type
    using base = KTfwd::singlepop<KTfwd::popgenmut,diploid_t>;
    //! The current generation.  Start counting from zero
    unsigned generation;
    //! Constructor takes number of diploids as argument
    singlepop_t(const unsigned & N) : base(N),generation(0)
    {
    }

    unsigned gen() const
    /*!
      \return current generation.

      This is mostly useful on the Cython side
    */
    {
      return generation;
    }
    unsigned popsize() const
    /*!
      \return current population size.

      This is mostly useful on the Cython side
    */
    {
      return N;
    }
    int sane() const
    /*!
      \return int(N == diploids.size())

      Useful on Cython side to check that N is 
      getting updated during simulations...
    */
    {
      return int(N == diploids.size());
    }
  };


  struct metapop_t : public KTfwd::metapop<KTfwd::popgenmut,diploid_t>
  /*!
    \brief Multi-deme object where mutations have single effect size and dominance.

    This is the C++ representation of a single-deme simulation where a KTfwd::popgenmut
    is the mutation type, and a custom diploid (fwdpy::diploid_t) is the diploid type.

    This type inherits from fwdpp's type KTfwd::singlepop, and the main documentation
    is found in the fwdpp reference manual. 

    Here is a brief overview of the most important members of the base class.

    The format is name: type, comment

    mutations: vector<popgenmut>, contains a mix of extinct, segregating, and possibly fixed variants, depending on the simulation type.
    mcounts: vector<unsigned>, is the same length as mutations, and records the count (frequency, as an integer) of each mutation.
    gametes: vector<gamete>, contains a mix of extinct and extant gametes
    diploids: vector<vector<diploid>>, each vector is a deme
    fixations: vector<popgenmut>, is filled by simulations that "prine" fixations when they occur
    fixation_times: vector<unsigned>, is filled by simulations that "prine" fixations when they occur

    \note Internally, fwdppy uses extinct elements in containers for "object recycling."

    Further:

    1. For each diploid, first and second are indexes into gametes
    2. Within each gamete, mutations and smutations contain indexes to "neutral" and "selected" mutations, resepectively, in mutations
  */
  {
    //! Typedef for base type
    using base = KTfwd::metapop<KTfwd::popgenmut,diploid_t>;
    //! Current generation.  Start counting from 0
    unsigned generation;
    //! Constructor takes list of deme sizes are aregument
    metapop_t( const std::vector<unsigned> &Ns ) : base(&Ns[0],Ns.size()), generation(0)
    {
    }

    //! Constructor takes list of deme sizes are aregument
    metapop_t( const std::initializer_list<unsigned> & Ns ) : base(Ns), generation(0)
    {
    }

    //! Construct from a fwdpy::singlepop_t
    metapop_t(const singlepop_t & p) : base(p),generation(p.generation)
    {
    }
    
    unsigned gen() const
    /*!
      \return current generation.
      
      This is mostly useful on the Cython side
    */
    {
      return generation;
    }
    std::vector<unsigned> popsizes() const
    /*!
      \return list of current deme sizese
      
      This is mostly useful on the Cython side
    */
    {
      return Ns;
    }
    int sane() const
    /*!
      \return true of demes[i].size() == Ns[i] for all i. Returns false otherwise.
    */
    {
      for(unsigned i=0;i<diploids.size();++i)
	{
	  if(diploids[i].size()!=Ns[i]) return 0;
	}
      return 1;
    }
    int size() const
    /*!
      \return Number of demes
    */
    {
      return int(diploids.size());
    }
  };

  //Types based on KTfwd::generalmut_vec  //! Typedef for gamete type

  //! Typedef for gamete container
  using gcont_gm_vec_t = std::vector<gamete_t>;

  struct singlepop_gm_vec_t :  public KTfwd::singlepop<KTfwd::generalmut_vec,diploid_t>
  /*!
    \brief Single-deme object where mutations contain vector<double> for internal data.
    
    See fwdpy::singlepop_t documentation for details, which are the same as for this type.
  */ 
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
