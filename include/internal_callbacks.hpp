/*! 
  \file internal_callbacks.hpp

  \brief Wrappers to namespace KTfwd::extension from fwdpp.  Used in fwdpy.internal module
*/
#include <fwdpp/extensions/callbacks.hpp>

namespace fwdpy {
  namespace internal {
    //! Callback wrapper
    void make_constant_s(KTfwd::extensions::shmodel * s, const double scoeff);
    //! Callback wrapper
    void make_uniform_s(KTfwd::extensions::shmodel * s, const double lo, const double hi);
    //! Callback wrapper
    void make_exp_s(KTfwd::extensions::shmodel * s, const double mean);
    //! Callback wrapper
    void make_gamma_s(KTfwd::extensions::shmodel * s, const double mean, const double shape);
    //! Callback wrapper
    void make_gaussian_s(KTfwd::extensions::shmodel * s, const double sd);
    //! Callback wrapper
    void make_constant_h(KTfwd::extensions::shmodel * s, const double h);
  }
}
