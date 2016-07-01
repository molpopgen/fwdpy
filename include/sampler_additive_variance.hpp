/*!
  \file sampler_additive_variance.hpp

  Estimates cumulative contribution of mutations to
  V(A).

  Big thanks to Jaleal Sanjak for help with the QR
  decomposition code.
*/
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
#include "types.hpp"
#include "sampler_base.hpp"

namespace fwdpy
{
  struct VAcum
  /*!
    Cumulative additive genetic varianace as a function
    of frequency over time.
  */
  {
    double freq,pss;
    unsigned generation,N;
    VAcum(double f, double c, unsigned g,unsigned n) : freq(f),
						       pss(c),
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
    regression_results(gsl_vector_ptr_t && s, std::vector<std::size_t> && p) :
      sums(std::move(s)),
      ucol_labels(std::move(p))
    {
    }
  };

  struct additive_variance : public sampler_base
  {
    using final_t = std::vector<VAcum>;

    virtual void operator()(const singlepop_t * pop, const unsigned generation)
    {
      call_operator_details(pop,generation);
    }
    virtual void operator()(const multilocus_t * pop, const unsigned generation)
    {
      call_operator_details(pop,generation);
    }

    final_t final() const
    {
      return VGcollection;
    }

    additive_variance() : buffer(std::vector<double>()),
			  Gbuffer(std::vector<double>()),
			  taubuffer(std::vector<double>()),
			  Qbuffer(std::vector<double>()),
			  Rbuffer(std::vector<double>()),
			  sumsbuffer(std::vector<double>()),
			  VGcollection(final_t())
    {
      buffer.reserve(10000);
    }
  private:
    std::vector<double> buffer,Gbuffer,taubuffer,Qbuffer,Rbuffer,sumsbuffer;
    final_t VGcollection;

    template<typename pop_t>
    inline void call_operator_details(const pop_t * pop,
				      unsigned generation)
    {
      if(pop->N != pop->diploids.size()) throw std::runtime_error("pop->N != pop->diploids.size(), " +
								  std::string(__FILE__) + ", line " +
								  std::to_string(__LINE__));
      auto mut_keys = get_mut_keys(pop);
      if(mut_keys.empty())
	{
	  //TODO: add an empty thing to the return value.
	  return;
	}
      //Genetic values for each diploid
      double VG;
      fillG(pop,Gbuffer,&VG);
      auto Gview = gsl_vector_view_array(Gbuffer.data(),pop->N);
      auto G = &Gview.vector; //this is a gsl_vector *
      
      //Get a vector of the mcounts corresponding to mut_kets
      std::vector<KTfwd::uint_t> mut_key_counts;
      for( const auto i : mut_keys ) mut_key_counts.emplace_back(pop->mcounts[i]);

      //Check if we need to reallocate
      if( std::size_t(pop->N)*(mut_keys.size()+1) > buffer.size() )
	{
	  buffer.resize(std::size_t(pop->N)*(mut_keys.size()+1));
	}
      std::size_t tda = buffer.size()/pop->N;
      auto genotypes_view = gsl_matrix_view_array_with_tda(buffer.data(),pop->N,mut_keys.size()+1,tda);
      auto genotypes = &genotypes_view.matrix;
      gsl_matrix_set_zero(genotypes);
      update_matrix_counts(pop,mut_keys,genotypes);
      auto ucol_labels = prune_matrix(genotypes,mut_keys,mut_key_counts);      
      auto DF = std::count(ucol_labels.begin(),ucol_labels.end(),1);

      //Now, do the regression
      taubuffer.resize(std::min(genotypes->size1,genotypes->size2));
      sumsbuffer.resize(genotypes->size1);
      Qbuffer.resize(genotypes->size1*genotypes->size1);
      Rbuffer.resize(genotypes->size1*genotypes->size2);

      auto tau = gsl_vector_view_array(taubuffer.data(),taubuffer.size());
      auto sums = gsl_vector_view_array(sumsbuffer.data(),genotypes->size1);
      auto Q = gsl_matrix_view_array_with_tda(Qbuffer.data(),genotypes->size1,genotypes->size1,genotypes->size1);
      auto R = gsl_matrix_view_array_with_tda(Rbuffer.data(),genotypes->size1,genotypes->size2,genotypes->size2);
      //QR decomposition...
      gsl_linalg_QR_decomp(genotypes, &tau.vector);
      //Get the Q and R matrix separately
      gsl_linalg_QR_unpack(genotypes, &tau.vector,&Q.matrix,&R.matrix);
      //Multiply t(Q) %*% b and store in sums
      gsl_blas_dgemv(CblasTrans, 1.0, &Q.matrix, G, 0.0, &sums.vector);
      //Now, remove elements corresponding to columns not used in regression
      for( auto i = ucol_labels.rbegin() ; i != ucol_labels.rend() ; ++i )
	{
	  if(!*i)
	    {
	      auto d = std::distance(ucol_labels.begin(),i.base())-1;
	      mut_keys.erase(mut_keys.begin()+d);
	      mut_key_counts.erase(mut_key_counts.begin()+d);
	    }
	}
      if(mut_keys.size() != std::size_t(DF)) throw std::runtime_error("removal error: " +
								      std::string(__FILE__) +
								      ", " +
								      std::to_string(__LINE__));
      double SumOfSquares=0.0;
      std::vector<double> vSumOfSquares;
      /*
	j=1 b/c first value in sums is for the origin,
	not the first column.  GSL inside baseball.
      */
      for(std::size_t i=0,j=1;i<ucol_labels.size();++i)
	{
	  if(ucol_labels[i])
	    {
	      auto s = gsl_vector_get(&sums.vector,j++);
	      auto p = gsl_sf_pow_int(s,2);
	      vSumOfSquares.push_back(p);
	      SumOfSquares+=p;
	    }
	}
      //residual sum of squares
      double RSS = std::accumulate(sums.vector.data+DF+1,sums.vector.data+sums.vector.size,
				   0.,[](double a,double b) { return a + gsl_sf_pow_int(b,2); });
      SumOfSquares += RSS; //add in the RSS.
      std::set<std::size_t> ucounts({mut_key_counts.begin(),mut_key_counts.end()});
      double ttl_rsq=0.0,ttl_adj_rsq=0.0;
      double DFsum = double(genotypes->size1-1);
      double DFn=double(genotypes->size1-1-DF);
      for(auto uc : ucounts)
	{
	  //Get all mutations in regression with this frequency
	  auto er = std::equal_range(mut_key_counts.begin(),mut_key_counts.end(),uc,
				     [](const std::size_t a,const std::size_t b)
				     {
				       return a>b;
				     });
	  //Get the sum of squares for this frequency bin
	  auto d1 = std::distance(mut_key_counts.begin(),er.first);
	  auto d2 = std::distance(mut_key_counts.begin(),er.second);
	  double ssuc = std::accumulate(vSumOfSquares.begin() + d1,
					vSumOfSquares.begin() + d2,
					0.0);
	  //This is r^2 for this frequency bin
	  double rsq = ssuc/SumOfSquares;
	  //This is adjusted r^2 for this bin
	  double a = (RSS+SumOfSquares-ssuc)/SumOfSquares;
	  double b = DFsum/(DFn + double(mut_key_counts.size()-std::distance(er.first,er.second)));
	  double adj_rsq = 1.0 - a*b;
	  VGcollection.emplace_back(VAcum(double(uc)/(2.0*double(pop->diploids.size())),rsq,generation,pop->diploids.size()));
	}
    }

    //regression_results regression_details(const gsl_vector_ptr_t & G,
    std::vector<std::size_t> prune_matrix(gsl_matrix * genotypes,
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
      if(mut_keys.size()!=mut_key_counts.size()) throw std::runtime_error("key sizes unequal");
      std::vector<size_t> column_labels(mut_keys.size(),1);
      unsigned identical = 0;
      for( std::size_t col = 1 ; col < genotypes->size2-1 ; ++col ) //skip column 0...
	{
	  if(column_labels[col-1])
	    {
	      auto c1 = gsl_matrix_const_column(genotypes,col);
	      for( std::size_t col2 = col+1 ; col2 < genotypes->size2 ; ++col2 )
		{
		  //columns can only be identical if mutations have same frequency!
		  //The -1 is b/c genotypes has 1 extra column
		  if(column_labels[col2-1] && mut_key_counts[col-1]==mut_key_counts[col2-1])
		    {
		      bool ndiff=false;
		      auto c2 = gsl_matrix_const_column(genotypes,col2);
		      for(std::size_t i=0;!ndiff&&i<c1.vector.size;++i)
			{
			  if(gsl_vector_get(&c1.vector,i)!=gsl_vector_get(&c2.vector,i))
			    ndiff=true;
			}
		      if(!ndiff)
			{
			  ++identical;
			  column_labels[col2-1]=0;
			}
		    }
		}
	    }
	}
      std::size_t NROW = genotypes->size1;
      std::size_t NCOL = genotypes->size2-identical;
      if(NCOL-1 != std::count(column_labels.begin(),column_labels.end(),1))
	{
	  throw std::runtime_error("NCOL incorrect "+
				   std::to_string(NCOL) + " " +
				   std::to_string(std::count(column_labels.begin(),column_labels.end(),1)) + " " +
				   std::to_string(std::count(column_labels.begin(),column_labels.end(),0)) + " " +
				   std::to_string(genotypes->size2));
	}
      if (identical)
	{
	  std::size_t n = 0;
	  std::vector<std::size_t> offsets;
	  for(std::size_t i=0;i<column_labels.size();++i)
	    {
	      if(column_labels[i])
		{
		  offsets.push_back(n);
		}
	      else
		{
		  offsets.push_back(0);
		  ++n;
		}
	    }
	  for(std::size_t i=0;i<offsets.size();++i)
	    {
	      if(offsets[i])
		{
		  auto c = gsl_matrix_const_column(genotypes,i+1);
		  gsl_matrix_set_col(genotypes,i+1-offsets[i],&c.vector);
		}
	    }
	  genotypes->size2 -= identical;
      	}
      return column_labels;
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
      gsl_matrix_ptr_t rv(gsl_matrix_alloc(pop->diploids.size(),1+mut_keys.size()));
      gsl_matrix_set_zero(rv.get()); //set all values to 0.

      update_matrix_counts(pop,mut_keys,rv.get());
      return rv;
    }

    template<typename pop_t>
    void update_row_details(gsl_matrix * m,
			    const typename pop_t::gamete_t & g,
			    const pop_t * pop,
			    const std::vector<std::size_t> & mut_keys,
			    const size_t row)
    {
      for( auto && k : g.smutations )
	{
	  if( pop->mcounts[k] < 2*pop->N ) //skip fixations!!!
	    {
	      if(!pop->mcounts[k]) throw std::runtime_error("extinct mutation encountered: " + std::string(__FILE__) + ", " + std::to_string(__LINE__));
	      auto i = std::find(mut_keys.begin(),mut_keys.end(),k);
	      if(i==mut_keys.end()) throw std::runtime_error("mutation key not found: " +
							     std::string(__FILE__) + ", "
							     + std::to_string(__LINE__) +", "
							     + "mcount = " + std::to_string(pop->mcounts[k]));
	      std::size_t col = std::distance(mut_keys.begin(),i);
	      if( col + 1 >= m->size2 ) throw std::runtime_error("second dimension out of range: " + std::string(__FILE__) + ", " + std::to_string(__LINE__));
	      auto mp = gsl_matrix_ptr(m,row,col+1);
	      *mp += 1.0; //update counts
	    }
	}
    }
    template<typename pop_t, typename diploid_t>
    void update_matrix_counts_details(gsl_matrix * m,
				      const pop_t * pop,
				      const std::vector<std::size_t> & mut_keys,
				      const diploid_t & dip,
				      const size_t row)
    {
      update_row_details(m,pop->gametes[dip.first],pop,mut_keys,row);
      update_row_details(m,pop->gametes[dip.second],pop,mut_keys,row);
    }

    template<typename pop_t>
    void update_matrix_counts(const pop_t * pop,
			      const std::vector<std::size_t> & mut_keys,
			      gsl_matrix * rv)
    {
      //Fill the matrix
      std::size_t row=0;
      for( const auto & dip : pop->diploids )
	{
	  gsl_matrix_set(rv,row,0,1.0); //set column 0 to a value of 1.0
	  update_matrix_counts_details(rv,pop,mut_keys,dip,row);
	  row++;
	}
    }

    template<typename pop_t>
    void fillG(const pop_t * pop,
	       std::vector<double> & Gbuffer,
	       double  * VG) //single-pop...
    /*!
      Returns vector of genetic values of each diploid.
    */
    {
      Gbuffer.clear();
      for(const auto & dip : pop->diploids)
	{
	  Gbuffer.push_back(dip.g);
	}
    }
  };

  template<>
  inline
  void additive_variance::fillG<multilocus_t>(const multilocus_t * pop,
					      std::vector<double> & Gbuffer,
					      double *VG)
  {
    Gbuffer.clear();
    for(const auto & dip : pop->diploids)
      {
	Gbuffer.push_back(dip[0].g);
      }
  }

  template<>
  inline void
  additive_variance::update_matrix_counts<multilocus_t>(const multilocus_t * pop,
							const std::vector<std::size_t> & mut_keys,
							gsl_matrix * rv)
  /*!
    Return a 0,1,2 matrix of counts of causative alleles in each diploid.

    Specialization for fwdpy::multilocus_t
  */
  {
    //Fill the matrix
    std::size_t row=0;
    for( const auto & dip : pop->diploids )
      {
	gsl_matrix_set(rv,row,0,1.0); //set column 0 to a value of 1.0
	for( const auto & locus : dip)
	  {
	    update_matrix_counts_details(rv,pop,mut_keys,locus,row);
	  }
	row++;
      }
  };
}

#endif
