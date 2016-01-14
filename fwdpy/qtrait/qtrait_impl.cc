/* 
   Models of quantitative traits
   Trait values are additive over 1, 1+hs, 1+2s, where s is a Gaussian deviate

   The infinitely-many sites stuff is an Cython/fwdpp-based re-implementation of the 
   code used to generate preliminary data for R01GM115564.
*/

#include <types.hpp>
#include <internal/internal.hpp>
#include <qtrait/rules.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/experimental/sample_diploid.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/popgenmut.hpp>

#include <qtrait/details.hpp>
#include <thread>
#include <algorithm>
#include <memory>
#include <limits>

#include <gsl/gsl_statistics_double.h>

using namespace std;

using poptype = fwdpy::singlepop_t;

namespace fwdpy
{
  namespace qtrait
  {
    void evolve_qtraits_t( GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,
			   const unsigned * Nvector,
			   const size_t Nvector_length,
			   const double mu_neutral,
			   const double mu_selected,
			   const double littler,
			   const double f,
			   const double sigmaE,
			   const double optimum,
			   const double VS,
			   const int track,			   
			   const fwdpy::internal::region_manager * rm)
    {
      std::vector<std::thread> threads(pops->size());
      for(unsigned i=0;i<pops->size();++i)
	{
	  threads[i]=std::thread(fwdpy::qtrait::qtrait_sim_details_t<qtrait_model_rules>,
				 gsl_rng_get(rng->get()),
				 pops->operator[](i).get(),
				 Nvector,Nvector_length,
				 mu_neutral,mu_selected,littler,f,sigmaE,optimum,track,
				 std::move(KTfwd::extensions::discrete_mut_model(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks)),
				 std::move(KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw)),
				 std::move(qtrait_model_rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length))));
	}
      for(unsigned i=0;i<threads.size();++i) threads[i].join();
    }

    //Get properties out from the population
    std::map<string,double> qtrait_pop_props( const fwdpy::singlepop_t * pop )
    {
      std::vector<double> VG,VE,wbar;
      
      //Get data from the diploids
      for( auto itr = pop->diploids.cbegin() ; itr != pop->diploids.cend() ; ++itr )
	{
	  VG.push_back( itr->g );
	  VE.push_back( itr->e );
	  wbar.push_back( itr->w );
	}

      //Find the "leading factor"
      double twoN = 2.*double(pop->diploids.size());
      //std::size_t max_exp = std::numeric_limits<std::size_t>::max();
      double mvexpl = 0.,
	leading_e=std::numeric_limits<double>::quiet_NaN(),
	leading_f=std::numeric_limits<double>::quiet_NaN(); 
      for(std::size_t i = 0 ; i < pop->mcounts.size() ; ++i )
	{
	  if(pop->mcounts[i])
	    {
	      auto n = pop->mcounts[i];
	      double p1=double(n)/twoN;
	      if (2.0*p1*(1.-p1)*std::pow(pop->mutations[i].s,2.) > mvexpl)
		{
		  mvexpl=2.0*p1*(1.-p1);
		  leading_e = pop->mutations[i].s;
		  leading_f = p1;
		}
	    }
	}
      // auto itr = std::max_element(pop->mutations.cbegin(),pop->mutations.cend(),
      // 				  [&twoN]( const poptype::mutation_t & m1,
      // 					   const poptype::mutation_t & m2 ) {
      // 				    double p1 = double(m1.n)/twoN,p2=double(m2.n)/twoN;
      // 				    return p1*(1.-p1)*std::pow(m1.s,2.) < p2*(1.-p2)*std::pow(m2.s,2.);
      // 				  });
      
      // double mvexpl = std::numeric_limits<double>::quiet_NaN(),
      // 	leading_e=std::numeric_limits<double>::quiet_NaN(),
      // 	leading_f=std::numeric_limits<double>::quiet_NaN();
      // if(itr != pop->mutations.end())
      // 	{
      // 	  mvexpl = 2.*(double(itr->n)/twoN)*(1.-(double(itr->n)/twoN))*std::pow(itr->s,2.);
      // 	  leading_e = itr->s;
      // 	  leading_f = double(itr->n)/twoN;
      // 	}
      map<string,double> rv;
      rv["VG"] = gsl_stats_variance(&VG[0],1,VG.size());
      rv["VE"] = gsl_stats_variance(&VE[0],1,VE.size());
      rv["H2"] = rv["VG"]/(rv["VG"]+rv["VE"]);
      rv["wbar"] = gsl_stats_mean(&wbar[0],1,wbar.size());
      rv["max_expl"] = mvexpl;
      rv["leading_e"] = leading_e;
      rv["leading_q"] = leading_f;
      return rv;
    }
  } //ns qtrait
} //ns fwdpy

  
