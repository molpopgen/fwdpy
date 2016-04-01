#include "metapop.hpp"
#include <fwdpp/sugar/demography.hpp>
namespace fwdpy
{
  size_t migpop(const size_t & source_pop, const gsl_rng * r, const double & mig_prob)
  {
  }

  void re_init_mpop( metapop_t * mpop, const singlepop_t * pop)
  {
  }

  void copy_deme( metapop_t * mpop, const std::size_t i )
  {
    KTfwd::copy_pop(*mpop,i);
  }
  
  void remove_deme( metapop_t * mpop, const std::size_t i )
  {
    KTfwd::remove_pop(*mpop,i);
  }
  
  void merge_demes(metapop_t  * mpop, const std::size_t i, const std::size_t j)
  {
    KTfwd::merge_pops(*mpop,i,j);
  }
  
  void split_deme(const gsl_rng * r, metapop_t * mpop, const std::size_t i, const unsigned N_new, const bool replacement )
  {
    KTfwd::split_pop(r,*mpop,i,N_new,replacement);
  }
  
  void admix_demes(const gsl_rng * r, metapop_t * mpop, const std::size_t i, const std::size_t j, const double prop_i, const bool replacement)
  {
    KTfwd::admix_pops(r,*mpop,i,j,prop_i,replacement);
  }
}
