/*
  This file provides lightweight wrappers
  around calls to GSL.

  The intendend use is to implement 
  standard models of distributions on selection/
  dominance effects of mutations.

  These structures are currently NOT
  exported via inst/include/foRward.

  Rather, the intended use is that they provide
  the implementation details of foRward::shmodel,
  which is found in inst/includefoRward/shmodel.hpp.

  Within the Rcpp world, these structures should be 
  transparent to the user, returned via:
  as<wrap><foRward::shmode>( SEXP )
 */
#ifndef __FORWARD_SHMODELWRAPPERS_HPP__
#define __FORWARD_SHMODELWRAPPERS_HPP__

#include <gsl/gsl_randist.h>
#include <cmath>

namespace fwdpy {
  struct constsh
  {
    double x;
    constsh(const double & __x) : x(__x)
    {
      if(!std::isfinite(x)) {
	Rcpp::stop("value must be finite");
	}
    }
    inline double operator()(gsl_rng *) const
    {
      return x;
    }
  };

  struct exps
  {
    double mean;
    exps(const double & m) : mean(m)
    {
      if(!std::isfinite(mean))
	{
	  Rcpp::stop("mean must be finite");
	}
      if(mean==0.)
	{
	  Rcpp::stop("mean must not equal 0");
	}
    }
    inline double operator()(gsl_rng * r) const
    {
      return gsl_ran_exponential(r,mean);
    }
  };

  struct uniformsh
  {
    double mn,mx;
    uniformsh(const double & __mn,
	     const double & __mx) : mn(__mn),mx(__mx)
    {
      if(!std::isfinite(mn) || !std::isfinite(mx))
	{
	  Rcpp::stop("min and max of range must both be finite");
	}
      if(mn>mx)
	{
	  Rcpp::stop("min must be <= max");
	}
    }
    inline double operator()(gsl_rng * r) const
    {
      return gsl_ran_flat(r,mn,mx);
    }
  };

  struct betash
  {
    double a,b,factor;
    betash(const double & __a,
	   const double & __b,
	   const double & __f) : a(__a),b(__b),factor(__f)
    {
      if(!std::isfinite(a) || a <= 0.)
	{
	  Rcpp::stop("a must be > 0.");
	}
      if(!std::isfinite(b) || b <= 0.)
	{
	  Rcpp::stop("b must be > 0.");
	}
      if(!std::isfinite(factor) || !(factor>0.))
	{
	  Rcpp::stop("scaling factor must be finite and > 0");
	}
    }
    inline double operator()(gsl_rng * r) const
    {
      return factor*gsl_ran_beta(r,a,b);
    }
  };

  struct gaussiansh
  {
    double sd;
    gaussiansh(const double & __sd) : sd(__sd)
    {
      if(sd == 0.) Rcpp::stop("sd must not equal 0");
      if(!std::isfinite(sd)) Rcpp::stop("sd must be finite");
    }
    inline double operator()(const gsl_rng * r) const
    {
      return gsl_ran_gaussian(r,sd);
    }
  };

  struct gammas
  {
    double mean,shape;
    gammas(const double & __m,
	   const double & __s ) : mean(__m),shape(__s)
    {
      if(!std::isfinite(mean))
	{
	  Rcpp::stop("mean and sign must both be finite");
	}
      if(shape <= 0)
	{
	  Rcpp::stop("shape must be >= 0");
	}
    }
    inline double operator()(const gsl_rng * r) const
    {
      return gsl_ran_gamma(r,shape,mean/shape);
    }
  };
}
#endif
