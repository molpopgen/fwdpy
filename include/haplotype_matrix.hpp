#ifndef FWDPY_HAPLOTYPE_MATRIX_HPP
#define FWDPY_HAPLOTYPE_MATRIX_HPP

#include <unordered_set>
#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>
#include <string>
#include <stdexcept>
#include "types.hpp"

namespace fwdpy
{
  struct haplotype_matrix
  {
    //!Indexes of neutral markers in matrix
    std::vector<std::size_t> n;
    //!Indexes of selected markers in matrix
    std::vector<std::size_t> s;
    //!Positions of neutral markers
    std::vector<double> np;
    //!Positions of selected markers
    std::vector<double> sp;
    //! Genetic value
    std::vector<double> G;
    //! Random value
    std::vector<double> E;
    //! fitness
    std::vector<double> w;
    std::size_t nrow,ncol_n,ncol_s;
  };

  template<typename mcont_t>
  std::pair<std::vector<std::size_t>,
	    std::vector<std::size_t>>
    get_mut_keys_common(const mcont_t & mutations,
			const std::unordered_set<std::size_t> & n,
			const std::unordered_set<std::size_t> & s,
			haplotype_matrix & hm)
  {
    std::pair<std::vector<std::size_t>,
	      std::vector<std::size_t>> rv;
    std::copy(n.begin(),n.end(),std::back_inserter(rv.first));
    std::copy(s.begin(),s.end(),std::back_inserter(rv.second));

    std::sort(rv.first.begin(),rv.first.end(),[&mutations](const std::size_t i, const std::size_t j) {
	return mutations[i].pos<mutations[j].pos;
      });
    std::sort(rv.second.begin(),rv.second.end(),[&mutations](const std::size_t i, const std::size_t j) {
	return mutations[i].pos<mutations[j].pos;
      });
    for( auto & ni : rv.first ) hm.np.push_back(mutations[ni].pos);
    for( auto & si : rv.second ) hm.sp.push_back(mutations[si].pos);
    return rv;
  }
				  
  template<typename dipvec_t,
	   typename gcont_t,
	   typename mcont_t>
  std::pair<std::vector<std::size_t>,
	    std::vector<std::size_t>>
    get_mut_keys( const dipvec_t & diploids,
		  const gcont_t & gametes,
		  const mcont_t & mutations,
		  const std::vector<unsigned> & mcounts,
		  const std::vector<std::size_t> diploids_sample,
		  haplotype_matrix & hm)
  /*!
    Workhorse function for a single deme
  */
  {
    std::unordered_set<std::size_t> n,s;
    for( auto & dip : diploids_sample )
      {
	if(dip >= diploids.size()) throw std::out_of_range("diploid index out of range");
	for(auto & m : gametes[diploids[dip].first].mutations)
	  {
	    if(mcounts[m]<2*diploids.size()) //avoid fixed variants
	      n.insert(m);
	  }
	for(auto & m : gametes[diploids[dip].second].mutations)
	  {
	    if(mcounts[m]<2*diploids.size())
	      n.insert(m);
	  } 
	for(auto & m : gametes[diploids[dip].first].smutations)
	  {
	    if(mcounts[m]<2*diploids.size())
	      s.insert(m);
	  }
	for(auto & m : gametes[diploids[dip].second].smutations)
	  {
	    if(mcounts[m]<2*diploids.size())
	      s.insert(m);
	  }
      }
    auto rv = get_mut_keys_common(mutations,n,s,hm);
    hm.nrow=2*diploids_sample.size();
    hm.ncol_n=rv.first.size();
    hm.ncol_s=rv.second.size();
    return rv;
  }

  template<typename dipvec_t,
	   typename gcont_t,
	   typename mcont_t>
  std::pair<std::vector<std::size_t>,
	    std::vector<std::size_t>>
    get_mut_keys_mloc( const dipvec_t & diploids,
		       const gcont_t & gametes,
		       const mcont_t & mutations,
		       const std::vector<unsigned> & mcounts,
		       const std::vector<std::size_t> diploids_sample,
		       haplotype_matrix & hm)
  {
    std::unordered_set<std::size_t> n,s;
    for( auto & dip : diploids_sample )
      {
	if(dip >= diploids.size()) throw std::out_of_range("diploid index out of range");
	for( auto & locus : diploids[dip] )
	  {
	    for(auto & m : gametes[locus.first].mutations)
	      {
		if(mcounts[m]<2*diploids.size()) //avoid fixed variants
		  n.insert(m);
	      }
	    for(auto & m : gametes[locus.second].mutations)
	      {
		if(mcounts[m]<2*diploids.size())
		  n.insert(m);
	      } 
	    for(auto & m : gametes[locus.first].smutations)
	      {
		if(mcounts[m]<2*diploids.size())
		  s.insert(m);
	      }
	    for(auto & m : gametes[locus.second].smutations)
	      {
		if(mcounts[m]<2*diploids.size())
		  s.insert(m);
	      }
	  }
      }
    auto rv=get_mut_keys_common(mutations,n,s,hm);
    hm.nrow=2*diploids_sample.size();
    hm.ncol_n=rv.first.size();
    hm.ncol_s=rv.second.size();
    return rv;
  }
      
  inline void update_matrix(const std::size_t key,
			    const std::size_t row,
			    const std::vector<std::size_t> & keys,
			    std::vector<std::size_t> & indexes)
  {
    auto i = std::find(keys.begin(),keys.end(),key);
    if(i==keys.end()) throw std::runtime_error("fatal error: mutation key not found: " + std::to_string(key));
    indexes.emplace_back(std::size_t(row*keys.size() + std::distance(keys.begin(),i)));
  }

  template<typename dipvec_t,
	   typename gcont_t,
	   typename mcont_t>
  haplotype_matrix make_haplotype_matrix_single_deme_details(const dipvec_t & diploids,
							     const gcont_t & gametes,
							     const mcont_t & mutations,
							     const std::vector<unsigned> & mcounts,
							     const std::vector<std::size_t> diploids_sample)
  {
    haplotype_matrix rv;
    //Step 1: get mutation keys/positions, no. rows/columns
    auto keys = get_mut_keys(diploids,gametes,mutations,mcounts,diploids_sample,rv);

    //Step 2: fill matrices
    std::size_t row = 0;
    for( auto & dip : diploids_sample )
      {
	rv.G.push_back(diploids[dip].g);
	rv.E.push_back(diploids[dip].e);
	rv.w.push_back(diploids[dip].w);
	for(auto & m : gametes[diploids[dip].first].mutations)
	  {
	    if(mcounts[m]<2*diploids.size())
	      update_matrix(m,row,keys.first,rv.n);
	  }
	for(auto & m : gametes[diploids[dip].first].smutations)
	  {
	    if(mcounts[m]<2*diploids.size())
	      update_matrix(m,row,keys.second,rv.s);
	  }
	++row;
	for(auto & m : gametes[diploids[dip].second].mutations)
	  {
	    if(mcounts[m]<2*diploids.size())
	      update_matrix(m,row,keys.first,rv.n);
	  }
	for(auto & m : gametes[diploids[dip].second].smutations)
	  {
	    if(mcounts[m]<2*diploids.size())
	      update_matrix(m,row,keys.second,rv.s);
	  }
	++row;
      }
    if (row!=rv.nrow) throw std::runtime_error("fatal error: row != rv.nrow");
    return rv;
  }
  haplotype_matrix make_haplotype_matrix(const singlepop_t * pop,
					 const std::vector<std::size_t> & diploids)
  /*!
    Create a haplotype matrix from a single deme from a specified set of individuals
  */
  {
    return make_haplotype_matrix_single_deme_details(pop->diploids,
						     pop->gametes,
						     pop->mutations,
						     pop->mcounts,
						     diploids);
  }

  haplotype_matrix make_haplotype_matrix(const metapop_t * pop,
					 const std::vector<std::size_t> & diploids,
					 const std::size_t deme)
  {
    if(deme>pop->diploids.size())throw std::out_of_range("deme index out of range");
    return make_haplotype_matrix_single_deme_details(pop->diploids[deme],
						     pop->gametes,
						     pop->mutations,
						     pop->mcounts,
						     diploids);
  }
  haplotype_matrix make_haplotype_matrix(const multilocus_t * pop,
					 const std::vector<std::size_t> & diploids)
  /*!
    Create a haplotype matrix from a single deme (multi-locus) from a specified set of individuals
  */
  {
    haplotype_matrix rv;
    auto keys = get_mut_keys_mloc(pop->diploids,pop->gametes,pop->mutations,pop->mcounts,diploids,rv);
    std::size_t row=0;
    for (auto & dip : diploids)
      {
	rv.G.push_back(pop->diploids[dip][0].g);
	rv.E.push_back(pop->diploids[dip][0].e);
	rv.w.push_back(pop->diploids[dip][0].w);
	for( auto & locus : pop->diploids[dip] )
	  {
	    for( auto & m : pop->gametes[locus.first].mutations)
	      {
		if(pop->mcounts[m]<2*pop->diploids.size())
		  update_matrix(m,row,keys.first,rv.n);
	      }
	    for( auto & m : pop->gametes[locus.first].smutations)
	      {
		if(pop->mcounts[m]<2*pop->diploids.size())
		  update_matrix(m,row,keys.second,rv.s);
	      }
	    ++row;
	    for( auto & m : pop->gametes[locus.second].mutations)
	      {
		if(pop->mcounts[m]<2*pop->diploids.size())
		  update_matrix(m,row,keys.first,rv.n);
	      }
	    for( auto & m : pop->gametes[locus.second].smutations)
	      {
		if(pop->mcounts[m]<2*pop->diploids.size())
		  update_matrix(m,row,keys.second,rv.s);
	      }
	    ++row;
	  }
      }
    return rv;
  }
}

#endif
