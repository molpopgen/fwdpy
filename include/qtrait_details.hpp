#ifndef __FWDPY_QTRAIT_DETAILS_HPP__
#define __FWDPY_QTRAIT_DETAILS_HPP__


#include <fwdpp/experimental/sample_diploid.hpp>
#include <fwdpp/extensions/regions.hpp>
#include "types.hpp"

//Namespace pollution due to Cython/struct/namespace bug:
struct qtrait_sample_info_t
{
  KTfwd::sep_sample_t genotypes;
  std::vector<std::pair<double, double> > sh;
  template<typename T1,typename T2>
  qtrait_sample_info_t( T1 && a, T2 && b ) noexcept :
					    genotypes(std::forward<T1>(a)),
						sh(std::forward<T2>(b))
  {
  }
  qtrait_sample_info_t() noexcept : genotypes(KTfwd::sep_sample_t()),
				    sh(std::vector<std::pair<double,double>>())
  {
  }
};

namespace fwdpy {
  namespace qtrait {
  }
} //namespace
#endif
