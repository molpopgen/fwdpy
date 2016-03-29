#ifndef __FWDPY_HOCRULES_HPP__
#define __FWDPY_HOCRULES_HPP__

#include <cmath>
#include <fwdpp/diploid.hh>
#include "types.hpp"
/*
  Custom "rules" policy for single-region House-of-Cards simulations.

  This file is partly KT having fun, but also an effort to stress-test fwdpp's
  "experimental" API.

  The advantage of this struct is:
  1. An offspring has its G and E values automatically assigned.  This allows us to record the
  exact fitnesses/heritabilities used in the simulation over time.
*/

namespace fwdpy
{
  namespace qtrait
  {
    struct qtrait_model_rules
    {
      mutable double wbar;
      const double sigE,optimum,VS;
      mutable std::vector<double> fitnesses;
      mutable KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
      qtrait_model_rules(const double & __sigE,
			 const double & __optimum,
			 const double & __VS,
			 const unsigned __maxN = 100000) :wbar(0.),
							  sigE(__sigE),
							  optimum(__optimum),
							  VS(__VS),
							  fitnesses(std::vector<double>(__maxN)),
							  lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr))
      {
      }

      qtrait_model_rules(qtrait_model_rules &&) = default;
      
      qtrait_model_rules(const qtrait_model_rules & rhs) : wbar(rhs.wbar),sigE(rhs.sigE),optimum(rhs.optimum),
							   VS(rhs.VS),fitnesses(rhs.fitnesses),
							   lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr))
      {
	if(!fitnesses.empty())
	  lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(fitnesses.size(),&fitnesses[0]));
      }
      template<typename T,
	       typename gcont_t,
	       typename mcont_t,
	       typename ff>
      void w( const T & diploids,
	      gcont_t & gametes,
	      const mcont_t & mutations,
	      const ff & fitness_func ) const
      {
	auto N_curr = diploids.size();
	if(fitnesses.size() < N_curr) fitnesses.resize(N_curr);
	wbar = 0.;
	for( size_t i = 0 ; i < N_curr ; ++i )
	  {
	    gametes[diploids[i].first].n = gametes[diploids[i].second].n = 0;
	    fitnesses[i] = diploids[i].w;
	    wbar += diploids[i].w;
	  }
	wbar/=double(N_curr);
	lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(N_curr,&fitnesses[0]));
      }

      //! \brief Pick parent one
      inline size_t pick1(const gsl_rng * r) const
      {
	return gsl_ran_discrete(r,lookup.get());
      }
  
      //! \brief Pick parent 2.  Parent 1's data are passed along for models where that is relevant
      template<typename diploid_t,typename gcont_t,typename mcont_t>
      inline size_t pick2(const gsl_rng * r, const size_t & p1, const double & f,
			  const diploid_t &, const gcont_t &, const mcont_t &) const
      {
	return (f==1. ||(f>0.&&gsl_rng_uniform(r) < f)) ? p1 : gsl_ran_discrete(r,lookup.get());
      }
  
      //! \brief Update some property of the offspring based on properties of the parents
      template<typename diploid_t,typename gcont_t,typename mcont_t>
      void update(const gsl_rng * r,diploid_t & offspring, const diploid_t &, const diploid_t &,
		  const gcont_t & gametes, const mcont_t & mutations) const noexcept
      {
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
	offspring.e = gsl_ran_gaussian_ziggurat(r,sigE);
	double dev = (offspring.g+offspring.e-optimum);
	offspring.w = std::exp( -(dev*dev)/(2.*VS) );
	assert( std::isfinite(offspring.w) );
	return;
      }
    };

    struct gbr_model_rules
    /*
      "Gene-based" recessive model of Thornton, Foran, and Long (2013) PLoS Genetics
    */
    {
      mutable double wbar;
      const double sigE,optimum,VS;
      mutable std::vector<double> fitnesses;
      mutable KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
      gbr_model_rules(const double & __sigE,
			 const double & __optimum,
			 const double & __VS,
			 const unsigned __maxN = 100000) :wbar(0.),
							  sigE(__sigE),
							  optimum(__optimum),
							  VS(__VS),
							  fitnesses(std::vector<double>(__maxN)),
							  lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr))
      {
      }

      gbr_model_rules(gbr_model_rules &&) = default;
      
      gbr_model_rules(const gbr_model_rules & rhs) : wbar(rhs.wbar),sigE(rhs.sigE),optimum(rhs.optimum),
						     VS(rhs.VS),fitnesses(rhs.fitnesses),
						     lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr))
      {
	if(!fitnesses.empty())
	  lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(fitnesses.size(),&fitnesses[0]));
      }
      template<typename T,
	       typename gcont_t,
	       typename mcont_t,
	       typename ff>
      void w( const T & diploids,
	      gcont_t & gametes,
	      const mcont_t & mutations,
	      const ff & fitness_func ) const
      {
	auto N_curr = diploids.size();
	if(fitnesses.size() < N_curr) fitnesses.resize(N_curr);
	wbar = 0.;
	for( size_t i = 0 ; i < N_curr ; ++i )
	  {
	    gametes[diploids[i].first].n = gametes[diploids[i].second].n = 0;
	    fitnesses[i] = diploids[i].w;
	    wbar += diploids[i].w;
	  }
	wbar/=double(N_curr);
	lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(N_curr,&fitnesses[0]));
      }

      //! \brief Pick parent one
      inline size_t pick1(const gsl_rng * r) const
      {
	return gsl_ran_discrete(r,lookup.get());
      }
  
      //! \brief Pick parent 2.  Parent 1's data are passed along for models where that is relevant
      template<typename diploid_t,typename gcont_t,typename mcont_t>
      inline size_t pick2(const gsl_rng * r, const size_t & p1, const double & f,
			  const diploid_t &, const gcont_t &, const mcont_t &) const
      {
	return (f==1. ||(f>0.&&gsl_rng_uniform(r) < f)) ? p1 : gsl_ran_discrete(r,lookup.get());
      }
  
      //! \brief Update some property of the offspring based on properties of the parents
      template<typename diploid_t,typename gcont_t,typename mcont_t>
      void update(const gsl_rng * r,diploid_t & offspring, const diploid_t &, const diploid_t &,
		  const gcont_t & gametes, const mcont_t & mutations) const noexcept
      {
	offspring.g = std::accumulate( gametes[offspring.first].smutations.cbegin(),
				       gametes[offspring.first].smutations.cend(),0.0,
				       [&mutations](const double & d,const std::size_t & i) noexcept
				       {
					 return d + mutations[i].s;
				       } ) +
	  std::accumulate( gametes[offspring.second].smutations.cbegin(),
			   gametes[offspring.second].smutations.cend(),0.0,
			   [&mutations](const double & d,const std::size_t & i) noexcept
			   {
			     return d + mutations[i].s;
			   } );
	offspring.e = gsl_ran_gaussian_ziggurat(r,sigE);
	double dev = (offspring.g+offspring.e-optimum);
	offspring.w = std::exp( -(dev*dev)/(2.*VS) );
	assert( std::isfinite(offspring.w) );
	return;
      }
    };
  } //namespace qtrait
} //namespace fwdpy
#endif
