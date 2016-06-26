#ifndef FWDPY_FITNESS_MODELS_HPP
#define FWDPY_FITNESS_MODELS_HPP

#include "types.hpp"
#include <fwdpp/fitness_models.hpp>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace fwdpy
{
  //Single-region fitness function signatures
  /*! 
    Function pointer representing policy signature for site-dependent fitness models.
    The first argument is a non-const reference to the current "fitness".  The function
    must update this value appropriately given the data in the second argument. This is 
    for a single-region simulation
  */
  using genotype_fitness_updater = void(*)(double &, const mcont_t::value_type &);
  //! "Finalizer" for site-based fitness schemes for single region
  using fitness_function_finalizer = double(*)(double);
  //! "Finalizer" for haplotype-based fitness schemes for single region
  using haplotype_fitness_fxn_finalizer = double(*)(double,double);
  //! Policy signature for haplotype-dependent models for single region
  using haplotype_fitness_fxn = double(*)(const gamete_t &, const mcont_t &);

  //! C++11 signature for a fitness function for a single regions. Not exposed to Python (yet).
  using single_region_fitness_fxn = std::function<double(const diploid_t &,
							 const gcont_t &,
							 const mcont_t &)>;

  //Multi-locus fitness functions signatures
  
  //! C++11 signature for a multi-locus fitness function. Not exposed to Python (yet).
  using multi_locus_fitness_fxn = std::function<double(const std::vector<diploid_t> &,const gcont_t &, const mcont_t &)>;
  //! Function pointer representing policy for fitness functions for multi-locus simulations
  using mlocus_fitness_fxn = double(*)(const std::vector<diploid_t> &,const gcont_t &, const mcont_t &);
  
  struct site_dependent_fitness_wrapper
  /*!
    fwdpp's site_dependent_fitness::operator() is allowed to return values < 0, 
    making it usable for calculating "genetic values" in addition to "fitness". 
    However, that is bad if we're going to allow users to define custom fitness functions,
    hence this light wrapper.
  */
  {
    using result_type = double;
    template< typename diploid2dispatch,
	      typename gcont_t,
	      typename mcont_t,
	      typename fitness_updating_policy_hom,
	      typename fitness_updating_policy_het>
    inline result_type operator()( const diploid2dispatch & dip,
				   const gcont_t & gametes,
				   const mcont_t & mutations,
				   const fitness_updating_policy_hom & fpol_hom,
				   const fitness_updating_policy_het & fpol_het,
				   fitness_function_finalizer f,
				   const double & starting_fitness = 1. ) const noexcept
    {
      auto x = KTfwd::site_dependent_fitness()(gametes[dip.first].smutations.cbegin(),
					       gametes[dip.first].smutations.cend(),
					       gametes[dip.second].smutations.cbegin(),
					       gametes[dip.second].smutations.cend(),
					       mutations,fpol_hom,fpol_het,starting_fitness);
      return f(x);
    }
  };
  
  struct singlepop_fitness
  /*!
    Base class for fitness schemes for single-deme simulations
  */
  {
    using fitness_fxn_t = single_region_fitness_fxn;

    //! The fitness function itself
    fitness_fxn_t fitness_function;

    /*!
      Placeholder for future functionality
    */
    virtual void update(const singlepop_t *) {}

    //! Allows us to allocate on stack in Cython
    singlepop_fitness() : fitness_function(fitness_fxn_t()) {}
    //! Constructor is a sink for a fitness_fxn_t 
    singlepop_fitness(fitness_fxn_t ff) : fitness_function(std::move(ff)) {}
  };

  struct multilocus_fitness
  /*!
    Base class for fitness schemes for single-deme simulations 
    of multiple, partially-linked regions
  */
  {
    //! This is the form of a fitness function for a single-deme simulation
    using fitness_fxn_t = multi_locus_fitness_fxn;

    //! The fitness function itself
    fitness_fxn_t fitness_function;

    /*!
      Placeholder for future functionality
    */
    virtual void update(const multilocus_t *) {}

    //! Allows us to allocate on stack in Cython
    multilocus_fitness() : fitness_function(fitness_fxn_t()) {}
    //! Constructor is a sink for a fitness_fxn_t 
    multilocus_fitness(fitness_fxn_t ff) : fitness_function(std::move(ff)) {}
  };
    
  inline singlepop_fitness make_additive_fitness(double scaling = 2.0)
  /*!
    Standard additive model w/dominance for a single region
    
    Fitnesses are  1, 1+h*s, 1+scaling*s for AA,Aa, and aa
  */
  {
    return singlepop_fitness(std::bind(KTfwd::additive_diploid(),
				       std::placeholders::_1,
				       std::placeholders::_2,
				       std::placeholders::_3,
				       scaling));
  }
  
  inline singlepop_fitness make_multiplicative_fitness(double scaling = 2.0)
  /*!
    Standard multiplicative model w/dominance for a single region

    Fitnesses are  1, 1+h*s, 1+scaling*s for AA,Aa, and aa
  */
  {
    return singlepop_fitness(std::bind(KTfwd::multiplicative_diploid(),
				       std::placeholders::_1,
				       std::placeholders::_2,
				       std::placeholders::_3,
				       scaling));
  }
  
  inline singlepop_fitness make_custom_fitness(genotype_fitness_updater Aa,
					       genotype_fitness_updater aa,
					       fitness_function_finalizer wfinal,
					       double starting_fitness)
  {
    return singlepop_fitness( std::bind(site_dependent_fitness_wrapper(),
					std::placeholders::_1,
					std::placeholders::_2,
					std::placeholders::_3,
					aa,Aa,wfinal,starting_fitness));
  }
  
  inline singlepop_fitness make_custom_haplotype_fitness(haplotype_fitness_fxn h,
							 haplotype_fitness_fxn_finalizer f)
  {
    return singlepop_fitness(std::bind(KTfwd::haplotype_dependent_fitness(),
				       std::placeholders::_1,
				       std::placeholders::_2,
				       std::placeholders::_3,
				       h,f));
  }
  
  inline multilocus_fitness make_mloc_additive_fitness(double scaling = 2.0)
  /*!
    Additive within loci w/dominance, and then additive across loci
  */
  {
    return multilocus_fitness([scaling](const std::vector<diploid_t> & diploid,
					const gcont_t & gametes,
					const mcont_t & mutations)
			      {
				double w = 0.0;
				for(const auto & locus : diploid)
				  {
				    w+= KTfwd::additive_diploid()(gametes[locus.first],
								  gametes[locus.second],
								  mutations);
				  }
				return w;
			      });
  }

  inline multilocus_fitness make_mloc_additive_trait(double scaling = 2.0)
  /*!
    Additive within loci w/dominance, and then additive across loci
  */
  {
    return multilocus_fitness([scaling](const std::vector<diploid_t> & diploid,
					const gcont_t & gametes,
					const mcont_t & mutations)
			      {
				double w = 0.0;
				for(const auto & locus : diploid)
				  {
				    w+= KTfwd::site_dependent_fitness()(gametes[locus.first],
									gametes[locus.second],
									mutations,
									[&](double & fitness,const KTfwd::popgenmut & mut) noexcept
									{
									  fitness += (1. + scaling*mut.s);
									},
									[](double & fitness,const KTfwd::popgenmut & mut) noexcept
									{
									  fitness += (1. + mut.h*mut.s);
									},
									0.);
				  }
				return w;
			      });
  }

  inline multilocus_fitness make_mloc_multiplicative_fitness(double scaling = 2.0)
  /*!
    Additive within loci w/dominance, and then additive across loci
  */
  {
    return multilocus_fitness([scaling](const std::vector<diploid_t> & diploid,
					const gcont_t & gametes,
					const mcont_t & mutations)
			      {
				double w = 0.0;
				for(const auto & locus : diploid)
				  {
				    w+= KTfwd::multiplicative_diploid()(gametes[locus.first],
									gametes[locus.second],
									mutations);
				  }
				return w;
			      });
  }

  inline multilocus_fitness make_mloc_multiplicative_trait(double scaling = 2.0)
  /*!
    Multiplicative within loci w/dominance, and then additive across loci
  */
  {
    return multilocus_fitness([scaling](const std::vector<diploid_t> & diploid,
					const gcont_t & gametes,
					const mcont_t & mutations)
			      {
				double w = 0.0;
				for(const auto & locus : diploid)
				  {
				    w+= KTfwd::site_dependent_fitness()(gametes[locus.first],
									gametes[locus.second],
									mutations,
									[&](double & fitness,const KTfwd::popgenmut & mut) noexcept
									{
									  fitness *= (1. + scaling*mut.s);
									},
									[&mutations](double & fitness,const KTfwd::popgenmut & mut) noexcept
									{
									  fitness *= (1. + mut.h*mut.s);
									},
									1.);
				  }
				return w-1.0;
			      });
  }

  inline multilocus_fitness make_mloc_gbr_trait()
  /*!
    "GBR" model within loci, additive across loci
  */
  {
    return multilocus_fitness([](const std::vector<diploid_t> & diploid,
				 const gcont_t & gametes,
				 const mcont_t & mutations)
			      {
				double w = 0.0;
				for(const auto & locus : diploid)
				  {
				    auto h1 = std::accumulate( gametes[locus.first].smutations.cbegin(),
							       gametes[locus.first].smutations.cend(),0.0,
							       [&mutations](const double & d,const std::size_t & i) noexcept
							       {
								 return d + mutations[i].s;
							       } );
				    auto h2 = std::accumulate( gametes[locus.second].smutations.cbegin(),
							       gametes[locus.second].smutations.cend(),0.0,
							       [&mutations](const double & d,const std::size_t & i) noexcept
							       {
								 return d + mutations[i].s;
							       } );
				    w += std::sqrt(h1*h2);
				  }
				return w;
			      });
  }

  inline multilocus_fitness make_mloc_power_mean_trait(const double SLp,const double MLp,
						       const std::vector<double> & SLd,
						       const std::vector<double> & MLd)
  {
    return multilocus_fitness([&](const std::vector<diploid_t> & diploid,
				  const gcont_t & gametes,
				  const mcont_t & mutations)
			      {
				double w = 0.0;
				std::size_t j = 0;
				for(const auto & locus : diploid)
				  {
				    auto h1 = std::accumulate(gametes[locus.first].smutations.cbegin(),
							      gametes[locus.first].smutations.cend(),0.0,
							      [&mutations](const double & d,const std::size_t & i) noexcept
							      {
								return d + mutations[i].s;
							      } );
				    auto h2 = std::accumulate(gametes[locus.second].smutations.cbegin(),
							      gametes[locus.second].smutations.cend(),0.0,
							      [&mutations](const double & d,const std::size_t & i) noexcept
							      {
								return d + mutations[i].s;
							      } );
				    w += MLd[j]*( std::pow( std::pow( (SLd[0]*std::pow(h1,SLp) + SLd[1]*std::pow(h2,SLp)),1./SLp), MLp ) );
				    j++;
				  }
				return std::pow(w,1./MLp);
			      });
  }

  inline multilocus_fitness make_mloc_custom_fitness(mlocus_fitness_fxn f)
  {
    return multilocus_fitness(f);
  }
}


#endif
