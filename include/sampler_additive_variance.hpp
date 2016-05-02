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
#include <vector>
#include <set>
#include <algorithm>
#include <memory>
#include <iostream>
#include "types.hpp"

namespace fwdpy
{
  struct VAcum
  /*!
    Cumulative additive genetic varianace as a function
    of frequency over time.
   */
  {
    double freq,cumsum;
    unsigned generation,N;
    VAcum(double f, double c, unsigned g,unsigned n) : freq(f),
						       cumsum(c),
						       generation(g),
						       N(n)
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

  struct regression_results
  {
    gsl_vector_ptr_t sums;
    std::vector<std::size_t> ucol_labels;
    std::size_t DF;
    regression_results(gsl_vector_ptr_t && s, std::vector<std::size_t> && p, std::size_t df) :
      sums(std::move(s)),
      ucol_labels(std::move(p)),
      DF(df)
    {
    }
  };

  struct additive_variance
  {
    using final_t = std::vector<VAcum>;

    template<typename pop_t>
    inline void operator()(const pop_t * pop,
			   unsigned generation)
    {
      auto mut_keys = get_mut_keys(pop);
      if(mut_keys.empty())
	{
	  //TODO: add an empty thing to the return value.
	  return;
	}
      //Genetic values for each diploid
      double VG;
      auto G = fillG(pop,&VG);
      //0,1,2 count of mutations affecting fitness
      auto genotypes = make_variant_matrix(pop,mut_keys);
      //Get a vector of the mcounts corresponding to mut_kets
      std::vector<KTfwd::uint_t> mut_key_counts;
      for( const auto i : mut_keys ) mut_key_counts.emplace_back(pop->mcounts[i]);

      auto results = regression_details(G,genotypes,mut_keys,mut_key_counts);
      double SumOfSquares=0.0;
      std::vector<double> vSumOfSquares;
      for(std::size_t i = 1; i < genotypes->size1 ; i++)
      	{
      	  auto s = gsl_vector_get(results.sums.get(),i);
      	  auto p = gsl_sf_pow_int(s,2);
      	  if(i)
      	    {
      	      vSumOfSquares.push_back(p);
      	      SumOfSquares+=p;
      	    }
      	}
      //residual sum of squares
      //double RSS = std::accumulate(results.sums->data+(results.sums->size-results.DF),results.sums->data+results.sums->size,
      //0.,[](double a,double b) { return a + gsl_sf_pow_int(b,2); });

      //Assign VA to each frequency bin.
      std::set<size_t> ucounts;
      for(std::size_t i=0;i<results.ucol_labels.size();++i)
	{
	  if(results.ucol_labels[i]) ucounts.insert(mut_key_counts[i]);
	}

      double VGcheck=0.0,adjrsqcheck=0.0;;
      for( const auto ui : ucounts ) //this starts w/rares
	{
	  //Have to count up the number of markers at this count that were used in the regression
	  unsigned cui=0;
	  std::vector<std::size_t> uindexes;
	  for(std::size_t i=0;i<results.ucol_labels.size();++i) //Filling this is the problem!!!!
	    {
	      if(results.ucol_labels[i] && mut_key_counts[i] == ui)
		{
		  ++cui;
		  uindexes.push_back(i);
		}
	    }
	  //Get SS for mutations at this count and those not at this count
	  double SSui=0.0,SSnotui=0.0;
	  for(std::size_t i=0;i<vSumOfSquares.size();++i)
	    {
	      if( std::binary_search(uindexes.begin(),uindexes.end(),i) ) SSui += vSumOfSquares[i];
	      else SSnotui += vSumOfSquares[i];
	    }
	  double adj_rsq = SSui/SumOfSquares;
	  adjrsqcheck += adj_rsq;
	  VGcollection.emplace_back(VAcum(double(ui)/(2.0*double(pop->diploids.size())),adjrsqcheck,generation,pop->diploids.size()));
	  VGcheck += adj_rsq*VG;
	}
      std::cerr << generation << ' ' << adjrsqcheck << ' ' << VGcheck << ' ' << VG << '\n';
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

    regression_results regression_details(const gsl_vector_ptr_t & G,
					  const gsl_matrix_ptr_t & genotypes,
					  const std::vector<size_t> & mut_keys,
					  const std::vector<KTfwd::uint_t> & mut_key_counts)
    /*!
      Handles ugly details of the regression:
      1. Reduces genotypes just to the set of unique columns
      2. Does the QR regression
      3. Returns stuff

      \note The indexing in this code is tricky, and could be improved via a GSL matrix view skipping 1st column of genotypes.
     */
    {
      std::vector<size_t> column_labels(mut_keys.size(),1);
      unsigned identical = 0;
      for( std::size_t col = 1 ; col < genotypes->size2-1 ; ++col ) //skip column 0...
	{
	  //gsl_vector_view col_i = gsl_matrix_column(genotypes,col);
	  if(column_labels[col-1])
	    {
	      for( std::size_t col2 = col+1 ; col2 < genotypes->size2 ; ++col2 )
		{
		  //columns can only be identical if mutations have same frequency!
		  //The -1 is b/c genotypes has 1 extra columns
		  if(column_labels[col2-1] && mut_key_counts[col-1]==mut_key_counts[col2-1])
		    {
		      unsigned ndiff=0;
		      for( std::size_t row = 0 ; row < genotypes->size1 ; ++row )
			{
			  if(gsl_matrix_get(genotypes.get(),row,col) != gsl_matrix_get(genotypes.get(),row,col2))
			    ++ndiff;
			}
		      if(!ndiff)
			{
			  ++identical;
			  //column_labels[col2-1]=col-1; //This marks col2 as non-unique.
			  column_labels[col2-1]=0;
			}
		    }
		}
	    }
	}
      std::size_t NROW = genotypes->size1;
      std::size_t NCOL = genotypes->size2-identical;
      //Allocate placeholder variables
      gsl_vector_ptr_t tau(gsl_vector_alloc(NCOL));
      gsl_vector_ptr_t sums(gsl_vector_alloc(NROW));
      gsl_matrix_ptr_t Q(gsl_matrix_alloc(NROW,NROW));
      gsl_matrix_ptr_t R(gsl_matrix_alloc(NROW,NCOL));

      //FILE * fp1 = fopen("trait_values.txt","w");
      //FILE * fp2 = fopen("matrix.txt","w");
      if (identical)
	{
	  //Make new matrix of unique columns.
	  gsl_matrix_ptr_t ugeno(gsl_matrix_alloc(genotypes->size1,genotypes->size2-identical));
	  for(std::size_t i=0,j=0;i<genotypes->size2;++i)
	    {
	      //The -1 is b/c of the different lengths, as above...
	      if(i==0 ||(column_labels[i-1])) //Column is either column 0 or the column is unique
		{
		  for(std::size_t k = 0 ; k < genotypes->size1 ; ++k )
		    {
		      gsl_matrix_set(ugeno.get(),k,j,gsl_matrix_get(genotypes.get(),k,i));
		    }
		  ++j;
		}
	    }
	  // std::cerr << "type 1\n";
	  // for(std::size_t i=0;i<G->size;++i)
	  //   {
	  //     fprintf(fp1,"%lf\n",gsl_vector_get(G.get(),i));
	  //   }
	  // for(std::size_t i=0;i<ugeno->size1;++i)
	  //   {
	  //     for(std::size_t j=0;j<ugeno->size2;++j)
	  // 	{
	  // 	  fprintf(fp2,"%lf\t",gsl_matrix_get(ugeno.get(),i,j));
	  // 	}
	  //     fprintf(fp2,"\n");
	  //   }
	  // fclose(fp1);
	  // fclose(fp2);
	  // exit(1);
	  //QR decomposition...
	  gsl_linalg_QR_decomp(ugeno.get(), tau.get());
	  //Get the Q and R matrix separately
	  gsl_linalg_QR_unpack(ugeno.get(), tau.get(), Q.get(), R.get());
	  //Multiply t(Q) %*% b and store in sums
	  gsl_blas_dgemv(CblasTrans, 1.0, Q.get(), G.get(), 0.0, sums.get());
	  return regression_results(std::move(sums),
				    std::move(column_labels),
				    //This next value is the # of elements at end of sums
				    //from which RSS can be calculated.
				    NCOL-(2-identical));
	}
      // std::cerr << "type 2\n";
      // for(std::size_t i=0;i<G->size;++i)
      // 	{
      // 	  fprintf(fp1,"%lf\n",gsl_vector_get(G.get(),i));
      // 	}
      // for(std::size_t i=0;i<genotypes->size1;++i)
      // 	{
      // 	  for(std::size_t j=0;j<genotypes->size2;++j)
      // 	    {
      // 	      fprintf(fp2,"%lf\t",gsl_matrix_get(genotypes.get(),i,j));
      // 	    }
      // 	  fprintf(fp2,"\n");
      // 	}
      // fclose(fp1);
      // fclose(fp2);
      // exit(1);
      //QR decomposition...
      gsl_linalg_QR_decomp(genotypes.get(), tau.get());
      //Get the Q and R matrix separately
      gsl_linalg_QR_unpack(genotypes.get(), tau.get(), Q.get(), R.get());
      //Multiply t(Q) %*% b and store in sums
      gsl_blas_dgemv(CblasTrans, 1.0, Q.get(), G.get(), 0.0, sums.get());
      return regression_results(std::move(sums),
				std::move(column_labels),
				//This next value is the # of elements at end of sums
				//from which RSS can be calculated.
				NCOL-(2-identical));
    }

    template<typename pop_t>
    std::vector<std::size_t> get_mut_keys(const pop_t * pop)
    {
      std::vector<std::size_t> mut_keys; //array of keys for each segregating, non-neutral variant
      std::set<KTfwd::uint_t> ucounts; //the number of unique frequency bins...
      for(std::size_t i=0;i<pop->mutations.size();++i)
	{
	  //first check is to avoid extinct variants that fwdpp will recycle later.
	  //The second avoids fixed variants
	  if(pop->mcounts[i] && (pop->mcounts[i] < 2*pop->diploids.size()) && !pop->mutations[i].neutral)
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
      //auto x = pop->gametes[dip.first].smutations.size() + pop->gametes[dip.second].smutations.size();
      for( const auto k : pop->gametes[dip.first].smutations )
	{
	  if( pop->mcounts[k] < 2*pop->N ) //skip fixations!!!
	    {
	      auto i = std::find(mut_keys.begin(),mut_keys.end(),k);
	      std::size_t col = std::distance(mut_keys.begin(),i);
	      auto mp = gsl_matrix_ptr(m.get(),row,col+1);
	      *mp += 1.0; //update counts
	    }
	}
      for( const auto k : pop->gametes[dip.second].smutations )
	{
	  if( pop->mcounts[k] < 2*pop->N ) //skip fixations!!!
	    {
	      auto i = std::find(mut_keys.begin(),mut_keys.end(),k);
	      std::size_t col = std::distance(mut_keys.begin(),i);
	      auto mp = gsl_matrix_ptr(m.get(),row,col+1);
	      *mp += 1.0; //update counts
	    }
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
	gsl_matrix_set(rv.get(),row,0,1.0); //set column 0 to a value of 1.0
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
