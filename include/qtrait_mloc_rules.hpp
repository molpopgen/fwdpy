#ifndef FWDPY_QTRAIT_MLOC_RULES_HPP
#define FWDPY_QTRAIT_MLOC_RULES_HPP

#include <cmath>
#include <fwdpp/diploid.hh>

namespace fwdpy
{
  namespace qtrait
  {
    struct qtrait_mloc_rules
    {
      mutable double wbar;
      const double sigE,optimum,VS;
      mutable std::vector<double> fitnesses;

      mutable KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
      //! \brief Constructor
      qtrait_mloc_rules(const double & __sigE,
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

      qtrait_mloc_rules(qtrait_mloc_rules &&) = default;
      
      qtrait_mloc_rules(const qtrait_mloc_rules & rhs) : wbar(rhs.wbar),sigE(rhs.sigE),optimum(rhs.optimum),
							 VS(rhs.VS),fitnesses(rhs.fitnesses),
							 lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr))
      {
	if(!fitnesses.empty())
	  lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(fitnesses.size(),&fitnesses[0]));
      }
    
      //! \brief The "fitness manager"
      template<typename dipcont_t,
	       typename gcont_t,
	       typename mcont_t,
	       typename fitness_func>
      void w(const dipcont_t & diploids,
	     gcont_t & gametes,
	     const mcont_t & mutations,
	     const fitness_func & ) const
      {
	unsigned N_curr = diploids.size();
	if(fitnesses.size() < N_curr) fitnesses.resize(N_curr);
	wbar = 0.;

	for(unsigned i=0;i<N_curr;++i)
	  {
	    for ( auto region : diploids[i] )
	      {
		gametes[region.first].n=gametes[region.second].n=0;
	      }

	    //the g/e/w fields will be populated via update()
	    fitnesses[i] = diploids[i][0].w;
	    wbar += fitnesses[i];
	  }

	wbar /= double(diploids.size());

	/*!
	  Black magic alert:
	  fwdpp_internal::gsl_ran_discrete_t_ptr contains a std::unique_ptr wrapping the GSL pointer.
	  This type has its own deleter, which is convenient, because
	  operator= for unique_ptrs automagically calls the deleter before assignment!
	  Details: http://www.cplusplus.com/reference/memory/unique_ptr/operator=

	  This only works b/c the rhs of the expression below may be treated as an rvalue reference.
	*/
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
	return ((f==1.)||(f>0.&&gsl_rng_uniform(r) < f)) ? p1 : gsl_ran_discrete(r,lookup.get());
      }

      //! \brief Update some property of the offspring based on properties of the parents
      template<typename diploid_t,typename gcont_t,typename mcont_t>
      void update(const gsl_rng * r, diploid_t & offspring, const diploid_t & ,
		  const diploid_t &,
		  const gcont_t & gametes,
		  const mcont_t & mutations) const
      {
	offspring[0].g=0.;
	//The genetic model is additive w/in and b/w loci, so we do repeated calls to
	//fwdpp's additive fxn:
	for( const auto & locus : offspring )
	  {
	    offspring[0].g += KTfwd::site_dependent_fitness()(gametes[locus.first],gametes[locus.second],mutations,
							      [=](double & fitness,const typename mcont_t::value_type & mut) noexcept
							      {
								fitness += (2.*mut.s);
							      },
							      [](double & fitness,const typename mcont_t::value_type & mut) noexcept
							      {
								fitness += (mut.h*mut.s);
							      },
							      0.);
	  }
	offspring[0].e = gsl_ran_gaussian_ziggurat(r,sigE);
	double dev = (offspring[0].g+offspring[0].e-optimum);
	offspring[0].w = std::exp( -(dev*dev)/(2.*VS) );
	assert( std::isfinite(offspring[0].w) );
      }
    };

    struct qtrait_mloc_pm_rules
    {
      mutable double wbar;
      const double sigE,optimum,VS;
      const std::vector<double> SLd;
      const double SLp;
      const std::vector<double> MLd;
      const double MLp;
      mutable std::vector<double> fitnesses;

      mutable KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
      //! \brief Constructor
      qtrait_mloc_pm_rules(const double & __sigE,
			const double & __optimum,
			const double & __VS,
			const std::vector<double> & __SLd,
			const double & __SLp,
			const std::vector<double> & __MLd,
			const double & __MLp,
			const unsigned __maxN = 100000) :wbar(0.),
							 sigE(__sigE),
							 optimum(__optimum),
							 VS(__VS),
							 SLd(__SLd),
							 SLp(__SLp),
							 MLd(__MLd),
							 MLp(__MLp),
							 fitnesses(std::vector<double>(__maxN)),
							 lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr))
      {
      }

      qtrait_mloc_pm_rules(qtrait_mloc_pm_rules &&) = default;
      
      qtrait_mloc_pm_rules(const qtrait_mloc_pm_rules & rhs) : wbar(rhs.wbar),sigE(rhs.sigE),optimum(rhs.optimum),
							 VS(rhs.VS),SLd(rhs.SLd), SLp(rhs.SLp), MLd(rhs.MLd),MLp(rhs.MLp),
							 fitnesses(rhs.fitnesses),
							 lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr))
      {
	if(!fitnesses.empty())
	  lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(fitnesses.size(),&fitnesses[0]));
      }
    

      //! \brief The "fitness manager"
      template<typename dipcont_t,
	       typename gcont_t,
	       typename mcont_t,
	       typename fitness_func>
      void w(const dipcont_t & diploids,
	     gcont_t & gametes,
	     const mcont_t & mutations,
	     const fitness_func & ) const
      {
	unsigned N_curr = diploids.size();
	if(fitnesses.size() < N_curr) fitnesses.resize(N_curr);
	wbar = 0.;

	for(unsigned i=0;i<N_curr;++i)
	  {
	    for ( auto region : diploids[i] )
	      {
		gametes[region.first].n=gametes[region.second].n=0;
	      }

	    //the g/e/w fields will be populated via update()
	    fitnesses[i] = diploids[i][0].w;
	    wbar += fitnesses[i];
	  }

	wbar /= double(diploids.size());

	/*!
	  Black magic alert:
	  fwdpp_internal::gsl_ran_discrete_t_ptr contains a std::unique_ptr wrapping the GSL pointer.
	  This type has its own deleter, which is convenient, because
	  operator= for unique_ptrs automagically calls the deleter before assignment!
	  Details: http://www.cplusplus.com/reference/memory/unique_ptr/operator=

	  This only works b/c the rhs of the expression below may be treated as an rvalue reference.
	*/
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
	return ((f==1.)||(f>0.&&gsl_rng_uniform(r) < f)) ? p1 : gsl_ran_discrete(r,lookup.get());
      }

      //! \brief Update some property of the offspring based on properties of the parents
      template<typename diploid_t,typename gcont_t,typename mcont_t>
      void update(const gsl_rng * r, diploid_t & offspring, const diploid_t & ,
		  const diploid_t &,
		  const gcont_t & gametes,
		  const mcont_t & mutations) const
      {
	offspring[0].g=0.;
	//The genetic model is additive w/in and b/w loci, so we do repeated calls to
	//fwdpp's additive f fxn:
	double SL_G=0;
	size_t j=0;
	for( const auto & locus : offspring )
	  {

	  	double a=0;
	  	double b=0;
        for ( const auto & i : gametes[locus.first].smutations ){
          a += mutations[i].s;
        }
        for ( const auto & i : gametes[locus.second].smutations ){
          b += mutations[i].s;
        }
        //nested power mean = (sum ( weight e^p))^(1/p)
	  	SL_G+= MLd[j]*( pow( pow( (SLd[0]*pow(a,SLp) + SLd[1]*pow(b,SLp)),1./SLp), MLp ) );
	  	j+=1;
	  }
	offspring[0].g = pow(SL_G,1./MLp);
	offspring[0].e = gsl_ran_gaussian_ziggurat(r,sigE);
	double dev = (offspring[0].g+offspring[0].e-optimum);
	offspring[0].w = std::exp( -(dev*dev)/(2.*VS) );
	assert( std::isfinite(offspring[0].w) );
      }
    };
  }
}

#endif
