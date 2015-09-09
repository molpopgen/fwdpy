#include <config.h>
#include <deps.hpp>

using namespace std;

namespace fwdpy {
  vector<string> dependencies()
  {
    vector<string> rv;
    rv.push_back(string(FWDPP_VERSION));
    rv.push_back(string(LIBSEQ_VERSION));
    rv.push_back(string(GSL_VERSION));
    return rv;
  }

  vector<string> version()
  {
    return vector<string>(1,string(PACKAGE_VERSION));
  }
}

