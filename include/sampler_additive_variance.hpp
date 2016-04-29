#ifndef FWDPY_SAMPLER_ADDITIVE_VARIANCE_HPP
#define FWDPY_SAMPLER_ADDITIVE_VARIANCE_HPP

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_statistics_double.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <memory>
#include "types.hpp"

namespace fwdpy
{
  struct VGdata
  {
    double pos,esize,VG,SS,pSS,freq;
    unsigned origin,generation;
    VGdata(double p,double e,double vg,
	   double ss, double pss,
	   double q,
	   unsigned o,unsigned g) : pos(p),
				    esize(e),
				    VG(vg),
				    SS(ss),
				    pSS(pss),
				    freq(q),
				    origin(o),
				    generation(g)
    {
    }
  };

  //Enable use of smart pointers for GSL types
  struct gsl_matrix_deleter
  {
    void operator()( gsl_matrix * l ) noexcept 
    {
      gsl_matrix_free(l);
    }
  };

  struct gsl_vector_deleter
  {
    void operator()( gsl_vector * l ) noexcept 
    {
      gsl_vector_free(l);
    }
  };

  using gsl_vector_ptr_t = std::unique_ptr< gsl_vector, gsl_vector_deleter >;
  using gsl_matrix_ptr_t = std::unique_ptr< gsl_matrix, gsl_matrix_deleter >;
  
  struct additive_variance
  {
    using final_t = std::vector<VGdata>;
    
    template<typename pop_t>
    inline void operator()(const pop_t * pop,
			   unsigned generation)
    {
      if(generation%100==0.) std::cerr << generation << '\n';
      //Genetic values for each diploid
      double VG;
      auto G = fillG(pop,&VG);
      auto mut_keys = get_mut_keys(pop);
      //0,1,2 count of mutations affecting fitness
      auto genotypes = make_variant_matrix(pop,mut_keys);

      std::size_t NROW = pop->diploids.size();
      std::size_t NCOL = genotypes->size2; //Neat... 
      //Allocate placeholder variables
      gsl_vector_ptr_t tau(gsl_vector_alloc(NCOL));
      gsl_vector_ptr_t sums(gsl_vector_alloc(NROW));
      gsl_vector_ptr_t SS(gsl_vector_alloc(NROW));
      gsl_matrix_ptr_t Q(gsl_matrix_alloc(NROW,NROW));
      gsl_matrix_ptr_t R(gsl_matrix_alloc(NROW,NCOL));
      
      //QR decomposition...
      gsl_linalg_QR_decomp(genotypes.get(), tau.get());
      //Get the Q and R matrix separately
      gsl_linalg_QR_unpack(genotypes.get(), tau.get(), Q.get(), R.get());
      //Multiply t(Q) %*% b and store in sums
      gsl_blas_dgemv(CblasTrans, 1.0, Q.get(), G.get(), 0.0, sums.get());

      //square the sums
      double SumOfSquares=0.0;
      std::vector<double> vSumOfSquares;
      for(std::size_t i = 0; i < NROW; i++)
	{
	  auto s = gsl_vector_get(sums.get(),i);
	  auto p = gsl_sf_pow_int(s,2);
	  if(i)
	    {
	      vSumOfSquares.push_back(p);
	    }
	  SumOfSquares+=p;
	  //gsl_vector_set(SS.get(), i, p);
	}
      double VGcheck=0.0;
      for( std::size_t i=1 ; i < mut_keys.size() ; ++i )
	{
	  double pSS=vSumOfSquares[i]/SumOfSquares;
	  VGcollection.emplace_back(pop->mutations[mut_keys[i]].pos,
				    pop->mutations[mut_keys[i]].s,
				    pSS*VG,
				    vSumOfSquares[i],
				    pSS,
				    double(pop->mcounts[mut_keys[i]])/(2.0*double(pop->diploids.size())),
				    pop->mutations[mut_keys[i]].g,
				    generation
				    );
	  VGcheck+=pSS*VG;
	}
      std::cerr << VGcheck << ' ' << VG << ' ' << ' ' << SumOfSquares << '\n';
    }

    final_t final() const
    {
      return VGcollection;
    }

    additive_variance() : VGcollection(final_t())
    {
    }
  private:
    final_t VGcollection;
    template<typename pop_t>
    std::vector<std::size_t> get_mut_keys(const pop_t * pop)
    {
      std::vector<std::size_t> mut_keys; //array of keys for each segregating, non-neutral variant
      std::set<KTfwd::uint_t> ucounts; //the number of unique frequency bins...
      for(std::size_t i=0;i<pop->mutations.size();++i)
	{
	  //first check is to avoid extinct variants that fwdpp will recycle later.
	  if(pop->mcounts[i] && !pop->mutations[i].neutral)
	    {
	      mut_keys.push_back(i);
	      ucounts.insert(pop->mcounts[i]);
	    }
	}
      //Now, I need to sort based on frequency, descending order
      std::sort(mut_keys.begin(),mut_keys.end(),
		[&pop](std::size_t a, std::size_t b) { return pop->mcounts[a] > pop->mcounts[b]; });
      //Within each frequency class, sort in descending order via |effect size|...
      for(auto uc : ucounts)
	{
	  auto itr_b = std::find_if(mut_keys.begin(),mut_keys.end(),
				    [&pop,uc](std::size_t a) { return pop->mcounts[a]==uc; });
	  auto itr_e = std::find_if(itr_b+1,mut_keys.end(),
				    [&pop,uc](std::size_t a) { return pop->mcounts[a]!=uc; });
	  std::sort(itr_b,itr_e,[&pop](std::size_t a, std::size_t b){ return std::fabs(pop->mutations[a].s) > std::fabs(pop->mutations[b].s); });
	}
      return mut_keys;
    }
    
    
    template<typename pop_t>
    gsl_matrix_ptr_t make_variant_matrix(const pop_t * pop, const std::vector<std::size_t> & mut_keys)
    {
      return make_variant_matrix_details(pop,mut_keys);
    }

    template<typename pop_t, typename diploid_t>
    void update_matrix_counts_details(gsl_matrix_ptr_t & m,
				      const pop_t * pop,
				      const std::vector<std::size_t> & mut_keys,
				      const diploid_t & dip,
				      const size_t row)
    {
      for( const auto k : pop->gametes[dip.first].smutations )
	{
	  auto i = std::find(mut_keys.begin(),mut_keys.end(),k);
	  std::size_t col = std::distance(mut_keys.begin(),i);
	  auto mp = gsl_matrix_ptr(m.get(),row,col+1);
	  *mp += 1.0; //update counts
	}
      for( const auto k : pop->gametes[dip.second].smutations )
	{
	  auto i = std::find(mut_keys.begin(),mut_keys.end(),k);
	  std::size_t col = std::distance(mut_keys.begin(),i);
	  auto mp = gsl_matrix_ptr(m.get(),row,col+1);
	  *mp += 1.0; //update counts
	}
    }

    template<typename pop_t>
    void update_matrix_counts(const pop_t * pop,
			      const std::vector<std::size_t> & mut_keys,
			      gsl_matrix_ptr_t & rv)
    {
      return rv;
      //Fill the matrix
      std::size_t row=0;
      for( const auto & dip : pop->diploids )
	{
	  gsl_matrix_set(rv.get(),row,0,1.0); //set column 0 to a value of 1.0
	  update_matrix_counts_details(rv,pop,mut_keys,dip,row);
	  update_matrix_counts_details(rv,pop,mut_keys,dip,row);
	  row++;
	}
      return rv;
    }
    
    template<typename pop_t>
    gsl_matrix_ptr_t make_variant_matrix_details(const pop_t * pop, const std::vector<std::size_t> & mut_keys) //single pop -- specialized for multilocus pop
    /*!
      Return a 0,1,2 matrix of counts of causative alleles in each diploid.
    */
    {
      gsl_matrix_ptr_t rv(gsl_matrix_alloc(pop->diploids.size(),1+mut_keys.size()));
      gsl_matrix_set_zero(rv.get()); //set all values to 0.

      update_matrix_counts(pop,mut_keys,rv);
      return rv;
    };
    
    template<typename pop_t>
    gsl_vector_ptr_t fillG(const pop_t * pop, double  * VG) //single-pop...
    /*!
      Returns vector of genetic values of each diploid.
    */
    {
      gsl_vector_ptr_t g(gsl_vector_alloc(pop->diploids.size()));
      std::size_t i=0;
      for(const auto & dip : pop->diploids)
	{
	  gsl_vector_set(g.get(),i++,dip.g);
	}
      *VG = gsl_stats_variance(g->data,1,pop->diploids.size());
      return g;
    }
  };

  template<>
  inline
  gsl_vector_ptr_t additive_variance::fillG<multilocus_t>(const multilocus_t * pop,double *VG)
  {
    gsl_vector_ptr_t g(gsl_vector_alloc(pop->diploids.size()));
    std::size_t i=0;
    for(const auto & dip : pop->diploids)
      {
	gsl_vector_set(g.get(),i++,dip[0].g);
      }
    *VG = gsl_stats_variance(g->data,1,pop->diploids.size());
    return g;
  }

  template<>
  inline void
  additive_variance::update_matrix_counts<multilocus_t>(const multilocus_t * pop, const std::vector<std::size_t> & mut_keys,
							gsl_matrix_ptr_t & rv) 
  /*!
    Return a 0,1,2 matrix of counts of causative alleles in each diploid.

    Specialization for fwdpy::multilocus_t
  */
  {
    //Fill the matrix
    std::size_t row=0;
    for( const auto & dip : pop->diploids )
      {
	for( const auto & locus : dip)
	  {
	    update_matrix_counts_details(rv,pop,mut_keys,locus,row);
	    update_matrix_counts_details(rv,pop,mut_keys,locus,row);
	  }
	row++;
      }
  };
}

#endif
