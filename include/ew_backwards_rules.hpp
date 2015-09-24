#ifndef __FWDPY_EW_BACKWARDS_RULES_HPP__
#define __FWDPY_EW_BACKWARDS_RULES_HPP__

#include <gsl/gsl_statistics_double.h>
#include <cmath>
#include <types.hpp>
#include <fwdpp/diploid.hh>

/*
  Custom "rules" policy for the "backwards Eyre-Walker 2010" scheme
*/

namespace fwdpy
{
  namespace qtrait
  {
    struct qtrait_model_rules
    {
      mutable double wbar,tau,h2w,vwlocus;
      mutable std::vector<double> fitnesses,gterms;
      mutable KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
      const double optimum;
      qtrait_model_rules(const double __tau,
			 const double __h2w,
			 const double __optimum = 0.,
			 const unsigned __maxN = 100000) :wbar(0.),
							  tau(__tau),
							  h2w(__h2w),
							  vwlocus(0.),
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
	if(fitnesses.size() < N_curr) fitnesses.resize(N_curr);
	wbar = 0.;
	auto itr = diploids->cbegin();
	for( unsigned i = 0 ; i < N_curr ; ++i,++itr )
	  {
	    itr->first->n=0;
	    itr->second->n=0;
	    fitnesses[i] = itr->w;
	    gterms[i] = itr->g;
	    wbar += itr->w;
	  }
	assert(itr == diploids->cend());
	wbar/=double(N_curr);
	//Update variance in fitness due to this locus
	vwlocus = gsl_stats_variance(&gterms[0],1,N_curr);
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
	//Here, 'e' is used as the additional component of variance in fitness
	offspring->e = gsl_ran_gaussian_ziggurat(r,std::sqrt(vwlocus*(1.-h2w)/h2w));
	offspring->w = std::max(0.,offspring->g + offspring->e);
	return;
      }
    };
  } //ns qtrait
} //ns fwdpy


#endif

