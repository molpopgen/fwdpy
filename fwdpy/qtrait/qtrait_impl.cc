/* 
   Models of quantitative traits
   Trait values are additive over 1, 1+hs, 1+2s, where s is a Gaussian deviate

   The infinitely-many sites stuff is an Cython/fwdpp-based re-implementation of the 
   code used to generate preliminary data for R01GM115564.
*/

#include <types.hpp>
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

//Boost stuff--will have to add autoconf tests for later
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>

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
			   const double optimum,
			   const int track,
			   const std::vector<double> & nbegs,
			   const std::vector<double> & nends,
			   const std::vector<double> & nweights,
			   const std::vector<double> & sbegs,
			   const std::vector<double> & sends,
			   const std::vector<double> & sweights,
			   const std::vector<KTfwd::extensions::shmodel> * callbacks,
			   const std::vector<double> & rbeg,
			   const std::vector<double> & rend,
			   const std::vector<double> & rweight,
			   const char * fitness)
    {
    }
  
  } //ns qtrait
} //ns fwdpy



//'Get mutation trajectories from House-of-Cards simulation
//'@param ppop External pointer to simuated population
//'@param minsojourn Minimum sojourn time for a mutation. Don't include trajectories whose length is less than this value.
//'@param minfreq Don't include mutations that never get to at least this value.
//'@examples
//'sregions = list( new('s.gaussian',sd=0.1) )
//'recregions = list(new('region'))
//'rng = init.gsl.rng(101)
//'N=500
//'pop = qtrait.evolve(rng,rep(N,10*N),0,0.1/(2*N),10/(4*N),0.1,list(),sregions,recregions,track=TRUE)
//'pop.traj = qtrait.get.traj(pop[[1]],5)
// [[Rcpp::export("qtrait.get.traj")]]
// Rcpp::DataFrame get_qtrait_traj(SEXP ppop,const unsigned minsojourn = 0,const double minfreq = 0.)
// {
//   Rcpp::XPtr<poptype> pop(ppop);
//   std::vector<double> pos,freq,s;
//   std::vector<unsigned> generations;
//   /*
//     Key is origin, (pos,s), trajectories
//   */
//   //using trajtype = std::map< std::pair<unsigned,std::pair<double,double> >, std::vector<double> >;
//   for( poptype::trajtype::const_iterator itr = pop->trajectories.begin() ;
//        itr != pop->trajectories.end() ; ++itr )
//     {
//       double maxfreq = *std::max_element(itr->second.cbegin(),itr->second.cend());
//       if( itr->second.size() >= minsojourn && maxfreq >= minfreq )
// 	{
// 	  std::vector<unsigned> times(itr->second.size());
// 	  unsigned itime = itr->first.first;
// 	  std::generate(times.begin(),times.end(),[&itime]{ return itime++; });
// 	  generations.insert(generations.end(),times.begin(),times.end());
// 	  std::fill_n(std::back_inserter(pos),itr->second.size(),itr->first.second.first);
// 	  std::fill_n(std::back_inserter(s),itr->second.size(),itr->first.second.second);
// 	  std::copy(itr->second.begin(),itr->second.end(),std::back_inserter(freq));
// 	}
//     }
//   return Rcpp::DataFrame::create( Rcpp::Named("pos")=pos,
// 				  Rcpp::Named("generation") = generations,
// 				  Rcpp::Named("esize")=s,
// 				  Rcpp::Named("freq")=freq );				  
// }


//' Obtain a sample from a infinitely-many sites House-of-Cards simulation
//' @param ppop Pointer to a population returned by evolve.qtrait or evolve.qtrait.more
//' @param prng Pointer to a GSL rng
//' @param nsam The sample size
//' @param remove_fixed Remove variants where the derived state is present in all samples
// [[Rcpp::export("qtrait.sample")]]
// SEXP sample_qtrait( SEXP ppop,SEXP prng, const unsigned & nsam, const bool remove_fixed = true )
// {
//    Rcpp::XPtr<poptype> pop(ppop);
//    Rcpp::XPtr<foRward::GSLrng> rng(prng);
//    Rcpp::IntegerVector diplist(nsam);
//    for(unsigned i=0;i<nsam;++i)
//      {
//        diplist[i] = unsigned(gsl_ran_flat(rng->get(),0,pop->diploids.size()));
//      }
//    Rcpp::List view = foRward::build_diploid_view_details(*pop,diplist);
//    Rcpp::DataFrame neutral(view["neutral"]),selected(view["selected"]);
//    Rcpp::IntegerVector indlist_neutral = neutral["ind"],
//      indlist_selected = selected["ind"];
//    Rcpp::NumericVector G(indlist_neutral.size()),E(indlist_neutral.size()),w(indlist_neutral.size()),
//      G_s(indlist_selected.size()),E_s(indlist_selected.size()),w_s(indlist_selected.size()); //pheno data and fitness for this diploids
//    for( unsigned i = 0 ; i < indlist_neutral.size() ; ++i )
//      {
//        G[i] = pop->diploids[ diplist[indlist_neutral[i]-1] ].g;
//        E[i] = pop->diploids[ diplist[indlist_neutral[i]-1] ].e;
//        w[i] = pop->diploids[ diplist[indlist_neutral[i]-1] ].w;
//      }
//    for( unsigned i = 0 ; i < indlist_selected.size() ; ++i )
//      {
//        G_s[i] = pop->diploids[ diplist[indlist_selected[i]-1] ].g;
//        E_s[i] = pop->diploids[ diplist[indlist_selected[i]-1] ].e;
//        w_s[i] = pop->diploids[ diplist[indlist_selected[i]-1] ].w;
//      }
//    return Rcpp::List::create(Rcpp::Named("neutral") = Rcpp::DataFrame::create(Rcpp::Named("ind") = neutral["ind"],
// 									      Rcpp::Named("chrom") = neutral["chrom"],
// 									      Rcpp::Named("pos") = neutral["pos"],
// 									      Rcpp::Named("s") = neutral["s"],
// 									      Rcpp::Named("h") = neutral["h"],
// 									      Rcpp::Named("n") = neutral["n"],
// 									      Rcpp::Named("G") = G,
// 									      Rcpp::Named("E") = E,
// 									      Rcpp::Named("w") = w),
// 			     Rcpp::Named("selected") = Rcpp::DataFrame::create(Rcpp::Named("ind") = selected["ind"],
// 									       Rcpp::Named("chrom") = selected["chrom"],
// 									       Rcpp::Named("pos") = selected["pos"],
// 									       Rcpp::Named("s") = selected["s"],
// 									       Rcpp::Named("h") = selected["h"],
// 									       Rcpp::Named("n") = selected["n"],
// 									       Rcpp::Named("G") = G_s,
// 									       Rcpp::Named("E") = E_s,
// 									       Rcpp::Named("w") = w_s)
// 			     );
// }

// Summary stats of a population
// [[Rcpp::export("qtrait.pop.properties")]]
// SEXP qtrait_pop_props( SEXP ppop )
// {
//   Rcpp::XPtr<poptype> pop(ppop);

//   using namespace boost::accumulators;
//   accumulator_set< double, boost::accumulators::stats<tag::variance> > VG,VE;
//   accumulator_set< double, boost::accumulators::stats<tag::mean> > wbar,esize;

//   //Get data from the diploids
//   for( auto itr = pop->diploids.cbegin() ; itr != pop->diploids.cend() ; ++itr )
//     {
//       VG( itr->g );
//       VE( itr->e );
//       wbar( itr->w );
//     }

//   //Get data from the mutations
//   for( auto itr = pop->mutations.cbegin() ; itr != pop->mutations.cend() ; ++itr )
//     {
//       esize( itr->s );
//     }

//   //Find the "leading factor"
//   double twoN = 2.*double(pop->diploids.size());
//   auto itr = std::max_element(pop->mutations.cbegin(),pop->mutations.cend(),
// 			      [&twoN]( const poptype::mutation_t & m1,
// 				       const poptype::mutation_t & m2 ) {
// 				double p1 = double(m1.n)/twoN,p2=double(m2.n)/twoN;
// 				return p1*(1.-p1)*std::pow(m1.s,2.) < p2*(1.-p2)*std::pow(m2.s,2.);
// 			      });
  
//   double mvexpl = std::numeric_limits<double>::quiet_NaN(),leading_e=std::numeric_limits<double>::quiet_NaN(),leading_f=std::numeric_limits<double>::quiet_NaN();
//   if(itr != pop->mutations.end())
//     {
//       mvexpl = 2.*(double(itr->n)/twoN)*(1.-(double(itr->n)/twoN))*std::pow(itr->s,2.);
//       leading_e = itr->s;
//       leading_f = double(itr->n)/twoN;
//     }

//   return Rcpp::DataFrame::create( Rcpp::Named("generation") = pop->generation,
// 				  Rcpp::Named("VG") = variance(VG),
// 				  Rcpp::Named("VE") = variance(VE),
// 				  Rcpp::Named("wbar") = boost::accumulators::mean(wbar),
// 				  Rcpp::Named("mean.esize") =boost::accumulators::mean(esize),
// 				  Rcpp::Named("max.var.expl") = mvexpl,
// 				  Rcpp::Named("leading.esize") = leading_e,
// 				  Rcpp::Named("leading.freq") =  leading_f
// 				  );
// }

//' Effect size vs frequency
//' @param ppop External pointer to population.
//' @return A data.frame with all mutation effect sizes and frequencies.
// [[Rcpp::export("qtrait.esize.freq")]]
// SEXP qtrait_esize_freq(SEXP ppop)
// {
//   Rcpp::XPtr<poptype> pop(ppop);
//   double twoN = 2.*double(pop->diploids.size());
//   std::vector<double> esize,freq;
//   for( const auto & __m : pop->mutations )
//     {
//       esize.push_back(__m.s);
//       freq.push_back(double(__m.n)/twoN);
//     }
//   return Rcpp::DataFrame::create( Rcpp::Named("esize") = esize,
// 				  Rcpp::Named("freq") = freq );
// }
  
