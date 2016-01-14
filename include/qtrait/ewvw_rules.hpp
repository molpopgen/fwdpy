#ifndef __FWDPY_EWVW_RULES_HPP__
#define __FWDPY_EWVW_RULES_HPP__

#include <gsl/gsl_statistics_double.h>
#include <cmath>
#include <types.hpp>
#include <fwdpp/diploid.hh>
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
      mutable double wbar,vs_ttl,sigmae, optimum;;
      mutable std::vector<double> fitnesses;//,gterms;
      mutable KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
      ewvw_rules(const double __vs_ttl,
		 const double __sigmae,
		 const double __optimum = 0.,
		 const unsigned __maxN = 100000) :wbar(0.),
						  vs_ttl(__vs_ttl),
						  sigmae(__sigmae),
						  optimum(__optimum),
						  fitnesses(std::vector<double>(__maxN)),
						  //gterms(std::vector<double>(__maxN)),
						  lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr))
      {
      }

      template<typename dipcont_t,
	       typename gcont_t,
	       typename mcont_t,
	       typename ff>
      void w( const dipcont_t & diploids,
	      gcont_t & gametes,
	      const mcont_t & mutations,
	      const ff & fitness_func ) const
      {
	unsigned N_curr = diploids.size();
	if(fitnesses.size() < N_curr) fitnesses.resize(N_curr);
	wbar = 0.;
	for( unsigned i = 0 ; i < N_curr ; ++i )
	  {
	    gametes[diploids[i].first].n = gametes[diploids[i].second].n = 0;
	    fitnesses[i] = diploids[i].w;
	    wbar += diploids[i].w;
	  }
	wbar/=double(N_curr);
	lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(N_curr,&fitnesses[0]));
      }

      //! \brief Pick parent one
      inline size_t pick1(gsl_rng * r) const
      {
	return gsl_ran_discrete(r,lookup.get());
      }
  
      //! \brief Pick parent 2.  Parent 1's data are passed along for models where that is relevant
      template<typename diploid_t,typename gcont_t,typename mcont_t>
      inline size_t pick2(gsl_rng * r, const size_t & p1, const double & f,
			  const diploid_t &, const gcont_t &, const mcont_t &) const
      {
	return (f==1. ||(f>0.&&gsl_rng_uniform(r) < f)) ? p1 : gsl_ran_discrete(r,lookup.get());
      }
  
      //! \brief Update some property of the offspring based on properties of the parents
      template<typename diploid_t,typename gcont_t,typename mcont_t>
      void update(gsl_rng * r,diploid_t & offspring, const diploid_t &, const diploid_t &,
		  const gcont_t & gametes, const mcont_t & mutations) const noexcept
      {
	//'g' is the part of fitness due to current trait.
	//it is additive with dominance.
	offspring.g = KTfwd::site_dependent_fitness()(gametes[offspring.first],gametes[offspring.second],mutations,
						       [=](double & fitness,const typename mcont_t::value_type & mut) noexcept
						       {
							 fitness += (2.*mut.s);
						       },
						       [](double & fitness,const typename mcont_t::value_type & mut) noexcept
						       {
							 fitness += (mut.h*mut.s);
						       },
						       0.);
	offspring.e = gsl_ran_gaussian_ziggurat(r,sigmae);
	double dev = (offspring.g+offspring.e-optimum);
	offspring.w = std::exp(-(dev*dev)/(2.*vs_ttl));
	return;
      }
    };
  } //ns qtrait
} //ns fwdpy


#endif

