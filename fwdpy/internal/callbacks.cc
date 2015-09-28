#include <fwdpp/extensions/callbacks.hpp>

#include <functional>

namespace fwdpy
{
  namespace internal
  {
    void make_gamma_s(KTfwd::extensions::shmodel * s, const double mean, const double shape)
    {
      s->s = std::bind(KTfwd::extensions::gamma(mean,shape),std::placeholders::_1);
    }

    void make_uniform_s(KTfwd::extensions::shmodel * s, const double lo, const double hi)
    {
      s->s = std::bind(KTfwd::extensions::uniform(lo,hi),std::placeholders::_1);
    }
  
    void make_exp_s(KTfwd::extensions::shmodel * s, const double mean)
    {
      s->s = std::bind(KTfwd::extensions::exponential(mean),std::placeholders::_1);
    }
  
    void make_constant_s(KTfwd::extensions::shmodel * s, const double scoeff)
    {
      s->s = std::bind(KTfwd::extensions::constant(scoeff),std::placeholders::_1);
    }
  
    void make_gaussian_s(KTfwd::extensions::shmodel * s, const double sd)
    {
      s->s = std::bind(KTfwd::extensions::gaussian(sd),std::placeholders::_1);
    }
  
    void make_constant_h(KTfwd::extensions::shmodel * s, const double h)
    {
      s->h = std::bind(KTfwd::extensions::constant(h),std::placeholders::_1);
    }
  }
}
