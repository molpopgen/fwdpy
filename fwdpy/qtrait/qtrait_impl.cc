/* 
   Models of quantitative traits
   Trait values are additive over 1, 1+hs, 1+2s, where s is a Gaussian deviate

   The infinitely-many sites stuff is an Cython/fwdpp-based re-implementation of the 
   code used to generate preliminary data for R01GM115564.
*/

#include <types.hpp>
#include <fwdpy_internal.hpp>
#include <qtrait_rules.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/experimental/sample_diploid.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/popgenmut.hpp>

#include <qtrait_details.hpp>
#include <thread>
#include <algorithm>
#include <memory>
#include <limits>

#include <gsl/gsl_statistics_double.h>

//Libsequence
#include <Sequence/PolySIM.hpp>
//Note: requires libseq >= 1.8.5!!!!
#include <Sequence/PolyTableSlice.hpp>

using namespace std;
using namespace Sequence;

using poptype = fwdpy::singlepop_t;

namespace fwdpy
{
  namespace qtrait
  {
    /*
      Summary statistics -- these are copied from the qRHH repo that 
      was used for prelim. data for proposal.
    */

    struct GarudStats
    //http://arxiv.org/abs/1303.0906
    {
      double H1,H12,H2H1;
      GarudStats() : H1(1.),
		     H12(std::numeric_limits<double>::quiet_NaN()),
		     H2H1(std::numeric_limits<double>::quiet_NaN())
      {
      }
    };

    /*
      Garud et al. DOI: 10.1371/journal.pgen.1005004
      Messer & Petrov DOI: 10.1016/j.tree.2013.08.003
      Note that H1 = 1 - haplotype homozygosity, e.g. Depaulis and Veuille's "H"
    */
    GarudStats H1H12(const SimData & d)
    {
      if( d.empty() ) return GarudStats();
      set<string> uhaps(d.begin(),d.end());
      vector<double> hapcounts;
      GarudStats G;
      G.H1 = 0.;
      for_each(uhaps.cbegin(),uhaps.cend(),
	       [&](const string & hap) {
		 unsigned hapcount = count(d.begin(),d.end(),hap);
		 //Unbiased calc. of homozygosity in finite sample
		 G.H1 += double(hapcount)*double(hapcount-1)/(double(d.size())*double(d.size()-1));
		 hapcounts.push_back(double(hapcount));
	       });
      sort(hapcounts.begin(),hapcounts.end(),
	   std::bind(greater<double>(),std::placeholders::_1,std::placeholders::_2));
      G.H12 = G.H1 + 2.*hapcounts[0]*hapcounts[1]/pow(double(d.size()),2.);
      G.H2H1 = (G.H1-double(hapcounts[0]*(hapcounts[0]-1))/double(d.size()*(d.size()-1)))/G.H1;
      return G;
    }

    /*
      Mechanics of the nSL statistic
    */
    pair<double,double> __nlSsum(const unsigned & core,
				 const SimData & d,
				 const vector<size_t> & coretype,
				 const double * gmap = nullptr)
    {
      double s = 0.,s2=0.;
      unsigned nc=0u;
      for( unsigned i = 0 ; i < coretype.size() ; ++i )
	{
	  for( unsigned j = i+1 ; j < coretype.size() ; ++j )
	    {
	      auto right = mismatch(d[coretype[i]].cbegin()+core,d[coretype[i]].cend(),
				    d[coretype[j]].cbegin()+core);
	      string::const_reverse_iterator ri1(d[coretype[i]].cbegin()+core),
		ri2(d[coretype[j]].cbegin()+core);
	      auto left = mismatch(ri1,d[coretype[i]].crend(),ri2);
	      if(left.first != d[coretype[i]].rend() && right.first != d[coretype[i]].end())
		{
		  s += double(distance(left.first.base(),right.first) + 1);
		  s2 += (gmap == nullptr) ? fabs(d.position(distance(d[coretype[i]].cbegin(),right.first)-1) - 
						 d.position(distance(d[coretype[i]].cbegin(),left.first.base())))
		    : fabs( gmap[distance(d[coretype[i]].cbegin(),right.first)] -
			    gmap[distance(d[coretype[i]].cbegin(),left.first.base())] );
		  ++nc;
		}
	    }
	}
      return make_pair( s/double(nc),s2/double(nc) );
    }

    /*
      The nSL statistic of doi:10.1093/molbev/msu077
    */
    pair<double,double> nSL(const unsigned & core,
			    const SimData & d)
    {
      std::vector<size_t> der,anc;
      for(unsigned i=0;i<d.size();++i)
	{
	  if( d[i][core] == '1' ) der.push_back(i);
	  else anc.push_back(i);
	}
      pair<double,double> A = __nlSsum(core,d,anc),
	D = __nlSsum(core,d,der);
      return make_pair(log(A.first)-log(D.first),
		       log(A.second)-log(D.second));
    }

    /*
      Return max. abs value of unstandardized nSL and iHS, with the latter as defined by Ferrer-Admetella et al.
    */
    pair<double,double> mnSL(const SimData & d,const double minfreq)
    {
      if(d.empty()) return make_pair(std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN());
      vector<polymorphicSite> filtered;
      for_each(d.sbegin(),d.send(),
	       [&](const polymorphicSite & p) {
		 double f = double(count(p.second.begin(),p.second.end(),'1'))/double(d.size());
		 //if( min(f,1.-f) >= minfreq) {
		 if( f >= minfreq ) {
		   filtered.push_back(p);
		 }
	       });
      if(filtered.empty()) return make_pair(std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN());
      SimData __filtered( filtered.begin(), filtered.end() );
      double rv = 0.,rv2=0.;
      for(unsigned i=0;i<__filtered.numsites();++i) 
	{
	  pair<double,double> rvi = nSL(i,__filtered);
	  if( isfinite(rvi.first) && fabs(rvi.first) > fabs(rv) ) rv = rvi.first;
	  if( isfinite(rvi.second) && fabs(rvi.second) > fabs(rv2) ) rv2 = rvi.second;
	}
      return make_pair(rv,rv2);
    } 

    /*
      Return max. abs value of standardized nSL and iHS, with the latter as defined by Ferrer-Admetella et al.
    */
    pair<double,double> snSL(const SimData & d,const double minfreq, const double & binsize)
    {
      if(d.empty()) return make_pair(std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN());
      vector<polymorphicSite> filtered;
      for_each(d.sbegin(),d.send(),
	       [&](const polymorphicSite & p) {
		 double f = double(count(p.second.begin(),p.second.end(),'1'))/double(d.size());
		 if( min(f,1.-f) >= minfreq) {
		   filtered.push_back(p);
		 }
	       });
      if(filtered.empty()) return make_pair(std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN());
      SimData __filtered( filtered.begin(), filtered.end() );
      //Associate the stats with their DAFs
      vector< pair<double, pair<double,double> > > binning;
      for(unsigned i=0;i<__filtered.numsites();++i) 
	{
	  pair<double,double> rvi = nSL(i,__filtered);
	  unsigned dcount = count((__filtered.sbegin()+i)->second.begin(),
				  (__filtered.sbegin()+i)->second.end(),'1');
	  binning.push_back( make_pair( double(dcount)/double(__filtered.size()),
					rvi));
	}
      double rv = std::numeric_limits<double>::quiet_NaN(), rv2=std::numeric_limits<double>::quiet_NaN();
      //Now, bin, standardise, and move on...
      for( double l = minfreq ; l < 1. ; l += binsize )
	{
	  vector< pair<double, pair<double,double> > > thisbin;
	  double s1=0.,s2=0.,sq1=0.,sq2=0.;
	  copy_if( binning.begin(),binning.end(),
		   back_inserter(thisbin),[&]( const pair<double, pair<double,double> > & data )
		   {
		     s1 += (isfinite(data.second.first)) ? data.second.first : 0.;
		     s2 += (isfinite(data.second.second)) ? data.second.second : 0.;
		     sq1 += (isfinite(data.second.first)) ? pow(data.second.first,2.) : 0.;
		     sq2 += (isfinite(data.second.second)) ? pow(data.second.second,2.) : 0.;
		     return isfinite(data.second.first) && data.first >= l && data.first < l+binsize;
		   });
	  if( thisbin.size() > 1 ) //otherwise SD = 0, so there's nothing to standardize
	    {
	      double mean1 = s1/double(thisbin.size()),
		mean2 = s2/double(thisbin.size()),
		sd1 = pow(sq1/double(thisbin.size()) - pow(mean1,2.),0.5),
		sd2 = pow(sq2/double(thisbin.size()) - pow(mean2,2.),0.5);
	      for_each( thisbin.begin(), 
			thisbin.end(),
			[&](const pair<double, pair<double,double> > & data ) {
			  double z1 = (data.second.first-mean1)/sd1,
			    z2 = (data.second.second-mean2)/sd2;
			  //Get max abs val of each stat
			  if( isfinite(z1) && (!isfinite(rv) || fabs(z1) > fabs(rv)) ) rv = z1;
			  if( isfinite(z2) && (!isfinite(rv2) || fabs(z2) > fabs(rv2)) ) rv2 = z2;
			});
	    }
	}
      return make_pair(rv,rv2);
    }
    
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
      const KTfwd::extensions::discrete_mut_model m(rm->nb,rm->ne,rm->nw,rm->sb,rm->se,rm->sw,rm->callbacks);
      auto recmap = KTfwd::extensions::discrete_rec_model(rm->rb,rm->rw,rm->rw);
      std::vector<GSLrng_t> rngs;
      std::vector<qtrait_model_rules> rules;
      for(unsigned i=0;i<pops->size();++i)
	{
	  //Give each thread a new RNG + seed
	  rngs.emplace_back(GSLrng_t(gsl_rng_get(rng->get())) );
	  rules.emplace_back(qtrait_model_rules(sigmaE,optimum,VS,*std::max_element(Nvector,Nvector+Nvector_length)));
	}
      std::vector<std::thread> threads(pops->size());
      for(unsigned i=0;i<pops->size();++i)
	{
	  threads[i]=std::thread(fwdpy::qtrait::qtrait_sim_details_t<qtrait_model_rules>,
				 rngs[i].get(),
				 pops->operator[](i).get(),
				 Nvector,Nvector_length,
				 mu_neutral,mu_selected,littler,f,sigmaE,optimum,track,
				 std::cref(m),std::cref(recmap),
				 std::ref(rules[i]));
	}
      for(unsigned i=0;i<threads.size();++i) threads[i].join();
    }

    //Get properties out from the population
    std::map<string,double> qtrait_pop_props( const fwdpy::singlepop_t * pop )
    {
      std::vector<double> VG,VE,wbar,esize;
      
      //Get data from the diploids
      for( auto itr = pop->diploids.cbegin() ; itr != pop->diploids.cend() ; ++itr )
	{
	  VG.push_back( itr->g );
	  VE.push_back( itr->e );
	  wbar.push_back( itr->w );
	}
      
      //Get data from the mutations
      for( auto itr = pop->mutations.cbegin() ; itr != pop->mutations.cend() ; ++itr )
	{
	  esize.push_back( itr->s );
	}
      
      //Find the "leading factor"
      double twoN = 2.*double(pop->diploids.size());
      auto itr = std::max_element(pop->mutations.cbegin(),pop->mutations.cend(),
				  [&twoN]( const poptype::mutation_t & m1,
					   const poptype::mutation_t & m2 ) {
				    double p1 = double(m1.n)/twoN,p2=double(m2.n)/twoN;
				    return p1*(1.-p1)*std::pow(m1.s,2.) < p2*(1.-p2)*std::pow(m2.s,2.);
				  });
      
      double mvexpl = std::numeric_limits<double>::quiet_NaN(),leading_e=std::numeric_limits<double>::quiet_NaN(),leading_f=std::numeric_limits<double>::quiet_NaN();
      if(itr != pop->mutations.end())
	{
	  mvexpl = 2.*(double(itr->n)/twoN)*(1.-(double(itr->n)/twoN))*std::pow(itr->s,2.);
	  leading_e = itr->s;
	  leading_f = double(itr->n)/twoN;
	}
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

    map<string,vector<double> > get_qtrait_traj(const singlepop_t *pop,const unsigned minsojourn,const double minfreq)
    {
      std::vector<double> pos,freq,s;
      std::vector<double> generations;
      /*
	Key is origin, (pos,s), trajectories
      */
      //using trajtype = std::map< std::pair<unsigned,std::pair<double,double> >, std::vector<double> >;
      for( poptype::trajtype::const_iterator itr = pop->trajectories.begin() ;
	   itr != pop->trajectories.end() ; ++itr )
	{
	  double maxfreq = *std::max_element(itr->second.cbegin(),itr->second.cend());
	  if( itr->second.size() >= minsojourn && maxfreq >= minfreq )
	    {
	      std::vector<unsigned> times(itr->second.size());
	      unsigned itime = itr->first.first;
	      std::generate(times.begin(),times.end(),[&itime]{ return itime++; });
	      generations.insert(generations.end(),times.begin(),times.end());
	      std::fill_n(std::back_inserter(pos),itr->second.size(),itr->first.second.first);
	      std::fill_n(std::back_inserter(s),itr->second.size(),itr->first.second.second);
	      std::copy(itr->second.begin(),itr->second.end(),std::back_inserter(freq));
	    }
	}
      map<string,vector<double>> rv;
      rv["pos"]=std::move(pos);
      rv["freq"]=std::move(freq);
      rv["generation"]=std::move(generations);
      rv["esize"]=std::move(s);
      return rv;
    }

    map<string,vector<double> > qtrait_esize_freq(const singlepop_t * pop)
    {
      double twoN = 2.*double(pop->diploids.size());
      std::vector<double> esize,freq;
      for( const auto & __m : pop->mutations )
	{
	  esize.push_back(__m.s);
	  freq.push_back(double(__m.n)/twoN);
	}
      map<string,vector<double>>rv;
      rv["esize"]=std::move(esize);
      rv["freq"]=std::move(freq);
      return rv;
    }
  } //ns qtrait
} //ns fwdpy

  
