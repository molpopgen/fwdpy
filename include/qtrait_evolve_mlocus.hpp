#ifndef FWDPY_QTRAIT_EVOLVE_MLOCUS_HPP
#define FWDPY_QTRAIT_EVOLVE_MLOCUS_HPP

#include <future>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <fwdpp/diploid.hh>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <fwdpp/experimental/sample_diploid_mloc.hpp>
#include "types.hpp"
#include "reserve.hpp"
#include "sampler_pop_properties.hpp"
#include "sampler_sample_n.hpp"
#include "sampler_selected_mut_tracker.hpp"

namespace fwdpy
{
  namespace qtrait
  {
    //Fitness function
    struct no_selection_multi
    {
      typedef double result_type;
      inline double operator()( const multilocus_t::dipvector_t::value_type &,
				const multilocus_t::gcont_t &,
				const multilocus_t::mcont_t & ) const
      {
	return 1.;
      }
    };

    //13 args...
    template<typename sampler,typename rules_type>
    inline typename sampler::final_t
    evolve_qtrait_mloc_sampler_details(fwdpy::multilocus_t * pop,
				       const unsigned long seed,
				       const unsigned * Nvector,
				       const size_t Nvector_len,
				       const std::vector<double> & neutral_mutation_rates,
				       const std::vector<double> & selected_mutation_rates,
				       const std::vector<double> & sigma_mus,
				       const std::vector<double> & within_region_rec_rates,
				       const std::vector<double> & between_region_rec_rates,
				       const double f,
				       const int interval,
				       sampler && isampler,
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
					 [&sigma_mus,&rng,i](){return gsl_ran_gaussian_ziggurat(rng,sigma_mus[i]);}, //effect size of non-neutral mutants
					 [](){return 1.;})); //dominance of non-neutral mutants
	  ++i;
	}

      //Total mutation rates
      auto tmu(neutral_mutation_rates);
      std::transform(tmu.begin(),tmu.end(),selected_mutation_rates.begin(),
		     tmu.begin(),
		     [](double a,double b){return a+b;});
    
      sampler s(std::move(isampler));
    
      //evolve...
      const unsigned simlen = unsigned(Nvector_len);
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
					      [](const gsl_rng * __r, const double __d){ return gsl_ran_binomial(__r,__d,1); },
					      std::bind(no_selection_multi(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
					      pop->neutral,
					      pop->selected,
					      f,
					      rules);
	  pop->N=nextN;
	}
      if (interval && pop->generation &&pop->generation%interval==0.)
	{
	  s(pop,pop->generation);
	}
      gsl_rng_free(rng);
      return s.final();
    }

    template<typename sampler,typename rules_t,class... Args>
    inline
    typename std::enable_if<std::is_void<typename sampler::final_t>::value,
			    void>::type
    evolve_qtrait_mloc_async_wrapper( GSLrng_t * rng,
				      std::vector<std::shared_ptr<multilocus_t> > * pops,
				      const unsigned * Nvector,
				      const size_t Nvector_len,
				      const std::vector<double> & neutral_mutation_rates,
				      const std::vector<double> & selected_mutation_rates,
				      const std::vector<double> & sigma_mus,
				      const std::vector<double> & within_region_rec_rates,
				      const std::vector<double> & between_region_rec_rates,
				      const double f,
				      const int sample,
				      const rules_t & rules,
				      Args&&... args)
    {
      using future_t = std::future<typename sampler::final_t>;
      std::vector<future_t> futures;
      for(std::size_t i=0;i<pops->size();++i)
	{
	  rules_t rules_thread(rules);
	  futures.emplace_back( std::async(std::launch::async,
					   evolve_qtrait_mloc_sampler_details<sampler,rules_t>,
					   pops->operator[](i).get(),
					   gsl_rng_get(rng->get()),
					   Nvector,
					   Nvector_len,
					   neutral_mutation_rates,
					   selected_mutation_rates,
					   sigma_mus,
					   within_region_rec_rates,
					   between_region_rec_rates,
					   f,
					   sample,
					   sampler(std::forward<Args>(args)...),
					   std::move(rules_thread)
				      )
				);	
	}
      for(std::size_t i=0;i<futures.size();++i ) futures[i].get();
    }

    template<typename sampler,typename rules_t,class... Args>
    inline
    typename std::enable_if<!std::is_void<typename sampler::final_t>::value,
			    std::vector<typename sampler::final_t> >::type
    evolve_qtrait_mloc_async_wrapper( GSLrng_t * rng,
				      std::vector<std::shared_ptr<multilocus_t> > * pops,
				      const unsigned * Nvector,
				      const size_t Nvector_len,
				      const std::vector<double> & neutral_mutation_rates,
				      const std::vector<double> & selected_mutation_rates,
				      const std::vector<double> & sigma_mus,
				      const std::vector<double> & within_region_rec_rates,
				      const std::vector<double> & between_region_rec_rates,
				      const double f,
				      const int sample,
				      const rules_t & rules,
				      Args&&... args)
    {
      using future_t = std::future<typename sampler::final_t>;
      std::vector<future_t> futures;
      for(std::size_t i=0;i<pops->size();++i)
	{
	  rules_t rules_thread(rules);
	  futures.emplace_back( std::async(std::launch::async,
					   evolve_qtrait_mloc_sampler_details<sampler,rules_t>,
					   pops->operator[](i).get(),
					   gsl_rng_get(rng->get()),
					   Nvector,
					   Nvector_len,
					   neutral_mutation_rates,
					   selected_mutation_rates,
					   sigma_mus,
					   within_region_rec_rates,
					   between_region_rec_rates,
					   f,
					   sample,
					   sampler(std::forward<Args>(args)...),
					   std::move(rules_thread)
					   )
				);	
	}
      std::vector<typename sampler::final_t> rv(futures.size());
      for(std::size_t i=0;i<futures.size();++i ) rv[i]=futures[i].get();
      return rv;
    }

    //CONCRETE FXNS BELOW THAT CYTHON CAN HANDLE
  
    //Take samples over time
    std::vector<sample_n<multilocus_t>::final_t>
    evolve_qtrait_mloc_sample_async( GSLrng_t * rng,
				     GSLrng_t * rng_sample,
				     std::vector<std::shared_ptr<multilocus_t> > * pops,
				     const unsigned * Nvector,
				     const size_t Nvector_length,
				     const std::vector<double> & neutral_mutation_rates,
				     const std::vector<double> & selected_mutation_rates,
				     const std::vector<double> & sigma_mus,
				     const std::vector<double> & within_region_rec_rates,
				     const std::vector<double> & between_region_rec_rates,
				     const double f,
				     const double sigmaE,
				     const double optimum,
				     const double VS,
				     const int sample,
				     const unsigned nsam);
  }
}
#endif
