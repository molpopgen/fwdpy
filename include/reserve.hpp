/*!
  \file reserve.hpp
  \brief Reserve memory for mutations and gametes in a population object
*/
#ifndef FWDPY_RESERVE_HPP
#define FWDPY_RESERVE_HPP

#include <cmath>
#include <type_traits>
#include <fwdpp/type_traits.hpp>

namespace fwdpy
{
  template<typename gcont_t,
	   typename mcont_t>
  void reserve_space( gcont_t & gametes,
		      mcont_t & mutations,
		      const unsigned N,
		      const double ttl_mutrate )
  {
    static_assert( KTfwd::traits::is_gamete_t<typename gcont_t::value_type>::value,
		   "gcont_t must be a container of gametes" );
    static_assert( KTfwd::traits::is_mutation_t<typename mcont_t::value_type>::value,
		   "mcont_t must be a container of mutations" );
    if( gametes.capacity() < std::size_t(4*N) )
      {
	gametes.reserve(std::size_t(4*N));
      }
    double theta = 4.*double(N)*ttl_mutrate;
    //Expected number of variants in a W-F population
    //under "standard" modeling assumptions.
    //The expression for ES is from Watterson 1975.
    double ES = std::log(2*N)*theta + (2./3.)*theta;
    if( mutations.capacity() < std::size_t(ES) )
      {
	mutations.reserve(std::size_t(ES));
      }
  }

}

#endif
