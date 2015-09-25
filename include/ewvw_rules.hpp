#ifndef __FWDPY_EWVW_RULES_HPP__
#define __FWDPY_EWVW_RULES_HPP__

#include <gsl/gsl_statistics_double.h>
#include <cmath>
#include <types.hpp>
#include <fwdpp/diploid.hh>
#include <iostream>
/*
  Custom "rules" policy for the "Eyre-Walker 2010" scheme where the trait is not 100% of variation in fitness

  v(w) at this trait is N(0,1)
  h2w is the proportion of variance in fitness due to this trait.  This is the "heritability in fitness due to variation in this trait".

  This allows for optimum shifts.
*/

namespace fwdpy
{
  namespace qtrait
  {
    struct ewvw_rules
    {
      gsl_rng * rng;
      mutable double wbar,h2w;
      mutable std::vector<double> fitnesses,gterms;
      mutable KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
      const double optimum;
      ewvw_rules(gsl_rng * r,
		 const double __h2w,
		 const double __optimum = 0.,
		 const unsigned __maxN = 100000) :rng(r),
						  wbar(0.),
						  h2w(__h2w),
						  fitnesses(std::vector<double>(__maxN)),
						  gterms(std::vector<double>(__maxN)),
						  lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr)),
						  optimum(__optimum)
      {
      }

      template<typename T,typename ff>
      void w( const T * diploids, const ff & fitness_func ) const
      {
	unsigned N_curr = diploids->size();
	if(fitnesses.size() < N_curr)
	  {
	    fitnesses.resize(N_curr);
	    gterms.resize(N_curr);
	  }

	auto itr = diploids->cbegin();
	for( unsigned i = 0 ; i < N_curr ; ++i,++itr )
	  {
	    itr->first->n=0;
	    itr->second->n=0;
	    //The fitness value at this trait is a unit Gaussian w.r.t. the optimum.
	    gterms[i] = itr->g;
	    //exp( -pow(itr->g - optimum,2.)/2. );
	  }
	assert(itr == diploids->cend());
	/*
	  Variance in fitness due to this trait.  It is critical that "G" values for 
	  this trait be normalized
	*/
	double mglocus = gsl_stats_mean(&gterms[0],1,N_curr);
	double sdglocus = gsl_stats_sd(&gterms[0],1,N_curr);
	std::transform(gterms.cbegin(),gterms.cbegin()+N_curr,
		       fitnesses.begin(),
		       [mglocus,sdglocus](const double g) {
			 return exp( -( pow((g-mglocus)/sdglocus,2.)/2.) );
		       });
	double vwlocus = gsl_stats_variance(&fitnesses[0],1,N_curr);
	/*
	  Rest of variance in fitness (effects of other loci in linkage equilibrium w/this trait
	  and/or V(E)
	*/
	double vwrest = vwlocus*(1.-h2w)/h2w;
	double sigma_rest = sqrt(vwrest);
	std::cerr << h2w << ' ' << mglocus << ' ' << sdglocus << ' ' << vwlocus << ' ' << ' ' << vwrest << ' ' << sigma_rest << ' ';
	double vwtot = vwlocus+vwrest;
	wbar = 0.;
	//Now, calculate the new, modified trait values
	std::transform(gterms.cbegin(),gterms.cbegin()+N_curr,
		       gterms.begin(),
		       [this,sigma_rest,vwtot](const double gi) {
			 double gttl = gi + gsl_ran_gaussian_ziggurat(this->rng,sigma_rest);
			 return gttl;
			 //double w = exp( -pow(gttl-optimum,2.)/(2.*vwtot) );
			 //assert( std::isfinite(w) );
			 //this->wbar += w;
			 //return w;
		       });
	//re-use variables for global fitness calcs
	mglocus = gsl_stats_mean(&gterms[0],1,N_curr);
	sdglocus = gsl_stats_sd(&gterms[0],1,N_curr);
	std::transform(gterms.cbegin(),gterms.cbegin()+N_curr,fitnesses.begin(),
		       [this,mglocus,sdglocus,vwlocus,vwtot](const double d) {
			 double w = exp( -pow((d-mglocus)/sdglocus,2.)/(2.*vwtot) );
			 this->wbar+=w;
			 return w;
		       });
	std::cerr <<"| "<< mglocus << ' ' << sdglocus << ' ' << gsl_stats_variance(&fitnesses[0],1,N_curr) << '\n';
	wbar /= double(N_curr);
	lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(N_curr,&fitnesses[0]));
      }

      //! \brief Pick parent one
      inline size_t pick1(gsl_rng * r) const
      {
	return gsl_ran_discrete(r,lookup.get());
      }
  
      //! \brief Pick parent 2.  Parent 1's data are passed along for models where that is relevant
      template<typename diploid_itr_t>
      inline size_t pick2(gsl_rng * r, const size_t & p1, diploid_itr_t, const double & f ) const
      {
	return (f==1. ||(f>0.&&gsl_rng_uniform(r) < f)) ? p1 : gsl_ran_discrete(r,lookup.get());
      }
  
      //! \brief Update some property of the offspring based on properties of the parents
      template<typename offspring_itr_t, typename parent_itr_t>
      void update(gsl_rng * r, offspring_itr_t offspring, parent_itr_t, parent_itr_t) const
      {
	//'g' is the part of fitness due to current trait.
	//it is additive with dominance.
	offspring->g = KTfwd::site_dependent_fitness()(offspring->first,offspring->second,
						       [=](double & fitness,const fwdpy::singlepop_t::mlist_t::const_iterator & mut)
						       {
							 fitness += (2.*mut->s);
						       },
						       [](double & fitness,const fwdpy::singlepop_t::mlist_t::const_iterator & mut)
						       {
							 fitness += (mut->h*mut->s);
						       },
						       0.);
	return;
      }
    };
  } //ns qtrait
} //ns fwdpy


#endif

