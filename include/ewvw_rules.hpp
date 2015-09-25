#ifndef __FWDPY_EWVW_RULES_HPP__
#define __FWDPY_EWVW_RULES_HPP__

#include <gsl/gsl_statistics_double.h>
#include <cmath>
#include <types.hpp>
#include <fwdpp/diploid.hh>

/*
  Custom "rules" policy for the "Eyre-Walker 2010" scheme where the trait is not 100% of variation in fitness
*/

namespace fwdpy
{
  namespace qtrait
  {
    struct ewvw_rules
    {
      gsl_rng * rng;
      mutable double wbar,tau,h2w;
      mutable std::vector<double> fitnesses,gterms,zi;
      mutable KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
      const double optimum;
      ewvw_rules(gsl_rng * r,
		 const double __tau,
		 const double __h2w,
		 const double __optimum = 0.,
		 const unsigned __maxN = 100000) :rng(r),
						  wbar(0.),
						  tau(__tau),
						  h2w(__h2w),
						  fitnesses(std::vector<double>(__maxN)),
						  gterms(std::vector<double>(__maxN)),
						  zi(std::vector<double>(__maxN)),
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
	    zi.resize(N_curr);
	  }

	auto itr = diploids->cbegin();
	for( unsigned i = 0 ; i < N_curr ; ++i,++itr )
	  {
	    itr->first->n=0;
	    itr->second->n=0;
	    gterms[i] = itr->g;
	  }
	//Get z-scores for genetic values in fitness
	double meang = gsl_stats_mean(&gterms[0],1,N_curr);
	double vwlocus = gsl_stats_variance(&gterms[0],1,N_curr);
	std::transform(gterms.begin(),gterms.begin()+N_curr,
		       zi.begin(),[meang,vwlocus](const double gi) {
			 return (gi-meang)/vwlocus;
		       });
	//variance in Z
	double vz = gsl_stats_variance(&zi[0],1,N_curr);
	double sigma = vz*(1.-h2w)/h2w;
	double sigmaw = sqrt(vz+pow(sigma,2.));
	//update the z-scores to reflect individual deviation w.r.to total variance in fitness
	//update fitnesses, too.
	wbar = 0.;
	std::transform(zi.begin(),zi.begin()+N_curr,
		       fitnesses.begin(),[this,sigma,sigmaw](const double zi) {
			 double zprime = zi + gsl_ran_gaussian_ziggurat(rng,sigma);
			 double w = gsl_ran_gaussian_pdf(zprime,sigmaw);
			 wbar += w;
			 return w;
		       });
	assert(itr == diploids->cend());
	wbar/=double(N_curr);
	//Update variance in fitness due to this locus

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

