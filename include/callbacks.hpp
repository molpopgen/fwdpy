#include <fwdpp/extensions/callbacks.hpp>

namespace fwdpy {
  void make_constant_s(KTfwd::extensions::shmodel * s, const double scoeff);
  void make_uniform_s(KTfwd::extensions::shmodel * s, const double lo, const double hi);
  void make_exp_s(KTfwd::extensions::shmodel * s, const double mean);
  void make_gamma_s(KTfwd::extensions::shmodel * s, const double mean, const double shape);
  void make_gaussian_s(KTfwd::extensions::shmodel * s, const double sd);
  void make_constant_h(KTfwd::extensions::shmodel * s, const double h);
}
