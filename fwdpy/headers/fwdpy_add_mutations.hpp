/*!
  \file fwdpy_add_mutations.hpp

  Wrappers for the variadic template functions
  in fwdpp/sugar/add_mutation.hpp
*/

#ifndef FWDPY_ADD_MUTATIONS_HPP
#define FWDPY_ADD_MUTATIONS_HPP

#include "types.hpp"

namespace fwdpy
{
    std::size_t add_mutation_cpp(singlepop_t *pop,
                                 const std::vector<std::size_t> &indlist,
                                 const std::vector<short> &clist,
                                 const double pos, const double s,
                                 const double h);
}

#endif
