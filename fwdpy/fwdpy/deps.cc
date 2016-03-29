#include <iostream>
#include "deps.hpp"

using namespace std;

namespace fwdpy {
  static const std::string fwdpp_citation = R"(@Article{BIBTEXKEY,
  author =   {K. R. Thornton},
  title =    {A C++ Template Library for Efficient Forward-Time Population Genetic Simulation of Large Populations},
  journal =      {Genetics},
  year =     {2014},
  OPTkey =   {},
  volume =   {198},
  OPTnumber =    {},
  pages =    {157-166},
  OPTmonth =     {},
  OPTnote =      {},
  annote =   {doi:/10.1534/genetics.114.165019}
})";

  vector<string> fwdpy_version()
  {
    return vector<string>(1,string(PACKAGE_VERSION));
  }

  void fwdpy_citation()
  {
    std::cout << "If you use fwdpy for your research, please cite the following:\n"
	      << "You are using fwdpy version " << PACKAGE_VERSION << '\n'
	      << "You should also cite the fwdpp paper, which is the brains behind this package\n"
	      << "You can get the fwdpp version used with the terminal command fwdpp-config --version\n"
	      << "The citation for fwdpp is:\n"
	      << fwdpp_citation
	      << '\n';
  }
}

