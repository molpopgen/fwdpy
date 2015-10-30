#include <config.h>
#include <deps.hpp>

using namespace std;

namespace fwdpy {
  vector<string> fwdpy_dependencies()
  {
    vector<string> rv;
    rv.push_back(string(FWDPP_VERSION));
    rv.push_back(string(GSL_VERSION));
    return rv;
  }

  vector<string> fwdpy_version()
  {
    return vector<string>(1,string(PACKAGE_VERSION));
  }
}

