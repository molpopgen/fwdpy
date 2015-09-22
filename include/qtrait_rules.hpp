#ifndef __FWDPY_HOCRULES_HPP__
#define __FWDPY_HOCRULES_HPP__

#include <types.hpp>
#include <fwdpp/diploid.hh>
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
      mutable std::vector<double> fitnesses;
      mutable KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
      const double sigE,optimum;
      qtrait_model_rules(const double & __sigE,
			 const double & __optimum,
			 const unsigned __maxN = 100000) :wbar(0.),
							  fitnesses(std::vector<double>(__maxN)),
							  lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr)),
							  sigE(__sigE),
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
	    wbar += itr->w;
	  }
	assert(itr == diploids->cend());
	wbar/=double(N_curr);
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
	offspring->e = gsl_ran_gaussian_ziggurat(r,sigE);
	double dev = (offspring->g+offspring->e-optimum);
	offspring->w = std::exp( -(dev*dev)/2. );
	assert( std::isfinite(offspring->w) );
	return;
      }
    };

  } //namespace qtrait
} //namespace fwdpy
#endif
