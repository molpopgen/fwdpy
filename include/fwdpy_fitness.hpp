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
  /*! 
    Function pointer representing policy signature for site-dependent fitness models.
    The first argument is a non-const reference to the current "fitness".  The function
    must update this value appropriately given the data in the second argument.
  */
  using genotype_fitness_updater = void(*)(double &, const mcont_t::value_type &);
  //! "Finalizer" for site-based fitness schemes
  using fitness_function_finalizer = double(*)(double);
  //! "Finalizer" for haplotype-based fitness schemes
  using hap_fitness_function_finalizer = double(*)(double,double);
  //! Polycy signature for haplotype-dependent models
  using haplotype_fitness_fxn = double(const gamete_t &, const mcont_t &);
  
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

  struct haplotype_dependent_fitness_wrapper
  {
    using result_type = double;
    template< typename diploid_t,
	      typename gcont_t,
	      typename mcont_t,
	      typename haplotype_policy,
	      typename diploid_policy >
    inline result_type operator()(const diploid_t & dip,
				  const gcont_t & gametes,
				  const mcont_t & mutations,
				  const haplotype_policy & hpol,
				  const diploid_policy & dpol) const noexcept
    {
      return KTfwd::haplotype_dependent_fitness()(gametes[dip.first],gametes[dip.second],
						  mutations,hpol,dpol);
    }
  };
  
  struct singlepop_fitness
  /*!
    Base class for fitness schemes for single-deme simulations
  */
  {
    //! This is the form of a fitness function for a single-deme simulation
    using fitness_fxn_t = std::function<double(const diploid_t &,
					       const gcont_t &,
					       const mcont_t &)>;

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
    using fitness_fxn_t = std::function<double(const std::vector<diploid_t> &,const gcont_t &, const mcont_t &)>;

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

  inline singlepop_fitness make_gbr_fitness()
  /*!
    "GBR" model of Thornton et al., 2013, for a single region
  */
  {
    return singlepop_fitness([](const diploid_t & diploid,
				const gcont_t & gametes,
				const mcont_t & mutations)
			     {
			       auto h1 = std::accumulate( gametes[diploid.first].smutations.cbegin(),
							  gametes[diploid.first].smutations.cend(),0.0,
							  [&mutations](const double & d,const std::size_t & i) noexcept
							  {
							    return d + mutations[i].s;
							  } );
			       auto h2 = std::accumulate( gametes[diploid.second].smutations.cbegin(),
							  gametes[diploid.second].smutations.cend(),0.0,
							  [&mutations](const double & d,const std::size_t & i) noexcept
							  {
							    return d + mutations[i].s;
							  } );
			       return std::sqrt(h1*h2);
			     });
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
							 hap_fitness_function_finalizer f)
  {
    return singlepop_fitness(std::bind(haplotype_dependent_fitness_wrapper(),
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

  inline multilocus_fitness make_mloc_multiplicative_fitness(double scaling = 2.0)
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
				    w+= KTfwd::multiplicative_diploid()(gametes[locus.first],
									gametes[locus.second],
									mutations);
				  }
				return w;
			      });
  }

  inline multilocus_fitness make_mloc_gbr_fitness()
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
}


#endif
