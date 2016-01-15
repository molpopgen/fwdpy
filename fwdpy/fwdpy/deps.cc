#include <deps.hpp>
#include <iostream>
using namespace std;

namespace fwdpy {
  vector<string> fwdpy_dependencies()
  {
    vector<string> rv{string(FWDPP_VERSION),string(GSL_VERSION)};
    return rv;
  }

  vector<string> fwdpy_version()
  {
    return vector<string>(1,string(PACKAGE_VERSION));
  }
}

