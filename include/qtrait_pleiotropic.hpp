#ifndef __FWDPY_QTRAITS_HPP__
#define __FWDPY_QTRAITS_HPP__


#include <vector>
#include <map>
#include <utility>
#include <string>
#include <limits>
#include <fwdpp/extensions/callbacks.hpp>
#include "types.hpp"
#include "internal_region_manager.hpp"
#include "qtrait_details.hpp"

namespace fwdpy
{
  namespace qtrait
  {    
    struct ew_mut_details
    {
      double s,e,p;
      ew_mut_details() : s(std::numeric_limits<double>::quiet_NaN()),
			 e(std::numeric_limits<double>::quiet_NaN()),
			 p(std::numeric_limits<double>::quiet_NaN())
      {
      }
      ew_mut_details(const double __s,
		     const double __e,
		     const double __p) : s(__s),e(__e),p(__p)
      {
      }
      ew_mut_details(const ew_mut_details &) = default;
      ew_mut_details(ew_mut_details &&) = default;
      ew_mut_details & operator=(const ew_mut_details &)=default;
    };

    std::map<double,ew_mut_details> ew2010_assign_effects(GSLrng_t * rng,
							  const fwdpy::singlepop_t * pop,
							  const double tau,
							  const double sigma);
    std::vector<double> ew2010_traits_cpp(const fwdpy::singlepop_t * pop,
					  const std::map<double,ew_mut_details> & effects);
  }
}

#endif
