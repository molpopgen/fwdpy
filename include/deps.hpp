/*!
  \file deps.hpp
  \brief Declare functions to return version numbers and citation info for Python package.
*/

#include <vector>
#include <string>

namespace fwdpy
{
  std::vector<std::string> fwdpy_version();
  void fwdpy_citation();
}
