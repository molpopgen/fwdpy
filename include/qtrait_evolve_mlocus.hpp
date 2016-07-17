#ifndef FWDPY_QTRAIT_EVOLVE_MLOCUS_HPP
#define FWDPY_QTRAIT_EVOLVE_MLOCUS_HPP

#include <future>
#include <vector>
#include <algorithm>
#include <functional>
#include <exception>
#include <set>
#include <type_traits>
#include <fwdpp/diploid.hh>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <fwdpp/experimental/sample_diploid_mloc.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include "types.hpp"
#include "reserve.hpp"
#include "sampler_base.hpp"
#include "fwdpy_fitness.hpp"
#include "sampler_base.hpp"
#include "fwdpp_features.hpp"

namespace fwdpy
{
  namespace qtrait
  {
    template<typename rules_type>
    inline void
    evolve_qtrait_mloc_cpp_details(fwdpy::multilocus_t * pop,
				   std::unique_ptr<multilocus_fitness> & fitness,
				   sampler_base & s,
				   const unsigned long seed,
				   const unsigned * Nvector,
				   const size_t Nvector_len,
				   const std::vector<double> & neutral_mutation_rates,
				   const std::vector<double> & selected_mutation_rates,
				   const std::vector<KTfwd::extensions::shmodel> & effects_dominance,
				   const std::vector<double> & within_region_rec_rates,
				   const std::vector<double> & between_region_rec_rates,
				   const double f,
				   const int interval,
				   rules_type && rules)
    {
      //TODO: check that vector input sizes ok.  Ideally, this needs to be done b4 this point, else we
      //have threads throwing exceptions...
    
      //Reserve space
      auto x = std::max_element(Nvector,Nvector+Nvector_len);
      assert(x!=Nvector+Nvector_len);
      reserve_space(pop->gametes,pop->mutations,*x,
		    std::accumulate(neutral_mutation_rates.begin(),neutral_mutation_rates.end(),0.) +
		    std::accumulate(selected_mutation_rates.begin(),selected_mutation_rates.end(),0.));

      //Get local rng 4 this thread
      gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
      gsl_rng_set(rng,seed);
    
      //Establish recombination maps--uniform w/in each locus
      std::vector<std::function<std::vector<double>(const multilocus_t::gamete_t &,
						    const multilocus_t::gamete_t &,
						    const multilocus_t::mcont_t &)> > recpols;
      unsigned i = 0;
      for( auto ri : within_region_rec_rates )
	{
	  recpols.emplace_back( std::bind(KTfwd::poisson_xover(),rng,ri,double(i),double(i+1),
					  std::placeholders::_1,std::placeholders::_2,std::placeholders::_3) );
	  ++i;
	}

      //Establish mutation models--uniform process w/in each locus,
      //w/mutations affecting trait value having DFE N(0,sigma_mus[i]) at the i-th locus
      i=0;
      std::vector<std::function<std::size_t(std::queue<std::size_t> &,multilocus_t::mcont_t &)> > mmodels;
      for( auto mi : neutral_mutation_rates )
	{
	  mmodels.emplace_back(std::bind(KTfwd::infsites(),
					 std::placeholders::_1,std::placeholders::_2,
					 rng,std::ref(pop->mut_lookup),&pop->generation,
					 mi, //mutation rate
					 selected_mutation_rates[i], //mutation rate
					 [&rng,i](){return gsl_ran_flat(rng,double(i),double(i+1));},   //mutation pos'n
					 [&effects_dominance,i,rng](){return effects_dominance[i].s(rng);},
					 [&effects_dominance,i,rng](){return effects_dominance[i].h(rng);}));
	  ++i;
	}

      //Total mutation rates
      auto tmu(neutral_mutation_rates);
      std::transform(tmu.begin(),tmu.end(),selected_mutation_rates.begin(),
		     tmu.begin(),
		     [](double a,double b){return a+b;});

      auto rules_local(std::forward<rules_type>(rules));
      //evolve...
      const unsigned simlen = unsigned(Nvector_len);
      fitness->update(pop);
      for( unsigned g = 0 ; g < simlen ; ++g, ++pop->generation )
	{
	  const unsigned nextN = *(Nvector+g);
	  if (interval && pop->generation &&pop->generation%interval==0.)
	    {
	      s(pop,pop->generation);
	    }
	  KTfwd::experimental::sample_diploid(rng,
					      pop->gametes,
					      pop->diploids,
					      pop->mutations,
					      pop->mcounts,
					      pop->N,nextN,
					      tmu.data(),
					      mmodels,
					      recpols,
					      between_region_rec_rates.data(),
					      //rec b/w loci is interpreted as cM!!!!!
					      [](const gsl_rng * __r, const double __d){return gsl_ran_bernoulli(__r,__d);},
					      fitness->fitness_function,
					      pop->neutral,
					      pop->selected,
					      f,
					      rules_local,
					      KTfwd::remove_neutral());
	  fwdpy::update_mutations_n(pop->mutations,pop->fixations,pop->fixation_times,pop->mut_lookup,pop->mcounts,pop->generation,2*nextN);
	  pop->N=nextN;
	  fitness->update(pop);
	}
      if (interval && pop->generation &&pop->generation%interval==0.)
	{
	  s(pop,pop->generation);
	}
      gsl_rng_free(rng);
      //Allow a sampler to clean up after itself
      s.cleanup();
    }
    
    void evolve_qtrait_mloc_cpp( GSLrng_t * rng,
				 std::vector<std::shared_ptr<multilocus_t> > * pops,
				 std::vector<std::unique_ptr<sampler_base> > & samplers,
				 const unsigned * Nvector,
				 const size_t Nvector_length,
				 const std::vector<double> & neutral_mutation_rates,
				 const std::vector<double> & selected_mutation_rates,
				 const std::vector<KTfwd::extensions::shmodel> & shmodels,
				 const std::vector<double> & within_region_rec_rates,
				 const std::vector<double> & between_region_rec_rates,
				 const double f,
				 const double sigmaE,
				 const double optimum,
				 const double VS,
				 const int interval,
				 const multilocus_fitness & fitness
				 );
  }
}
#endif
