/*!
  \file internal_region_manager.hpp
  \brief Details of how "regions" are implemented in fwdpy.
*/

#ifndef __FWDPY_INTERNAL_HPP__
#define __FWDPY_INTERNAL_HPP__

#include <vector>
#include <fwdpp/extensions/callbacks.hpp>

namespace fwdpy
{
  namespace internal
  {
    struct region_manager
    /*!
      \brief How regions are implemented in fwdpy.

      This type is used in fwdpy/internal/internal.pyx to convert user input
      to something that fwdpp can use.

      Within fwdpp, the members of this struct are bound to objects defined in
      fwdpp/extensions/regions.hpp

      \note These callbacks are for the case where a mutation 
      has a single 's' and 'h' value.
    */
    {
      //! Vector of callbacks for assiging s/h to mutations
      std::vector<KTfwd::extensions::shmodel> callbacks;
      /*
	Key to deciphering what is below:
	n* = neutral
	s* = selected
	r* = recombination
	*b = beginning of region
	*e = end of region
	*w = weight assigned to region
       */
      std::vector<double> nb,ne,nw,sb,se,sw,rb,re,rw;
      region_manager() : callbacks(std::vector<KTfwd::extensions::shmodel>()),
			 nb(std::vector<double>()),
			 ne(std::vector<double>()),
			 nw(std::vector<double>()),
			 sb(std::vector<double>()),
			 se(std::vector<double>()),
			 sw(std::vector<double>()),
			 rb(std::vector<double>()),
			 re(std::vector<double>()),
			 rw(std::vector<double>())
      {
      }
    };
  }
}

#endif
