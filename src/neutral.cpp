#include <neutral.hpp>
#include <map>

using namespace std;
using namespace KTfwd;

namespace fwdpy {

  void evolve_pop(GSLrng_t * rng,
		  singlepop_t * pop,
		  const unsigned & ngens,
		  const double & theta,
		  const double & rho)
  {
    double mu = theta/(4.*double(pop->N)),littler=rho/(4.*double(pop->N));

    std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng->get());

    for( unsigned generation = 0; generation < ngens; ++generation )
      {
	//Iterate the population through 1 generation
	double wbar = sample_diploid(rng->get(),
					    &pop->gametes,  //non-const pointer to gametes
					    &pop->diploids, //non-const pointer to diploids
					    &pop->mutations, //non-const pointer to mutations
					    pop->N,     //current pop size, remains constant
					    mu,    //mutation rate per gamete
					    std::bind(infsites(),rng->get(),&pop->mut_lookup,generation,
						      mu,0.,[&rng](){return gsl_rng_uniform(rng->get());},[](){return 0.;},[](){return 0.;}),
					    //The recombination policy includes the uniform crossover rate
					    std::bind(genetics101(),std::placeholders::_1,std::placeholders::_2,
						      std::placeholders::_3,
						      //Pass as reference
						      std::ref(pop->neutral),std::ref(pop->selected),
						      &pop->gametes,
						      littler,
						      rng->get(),
						      recmap),
					    std::bind(insert_at_end<singlepop_t::mutation_t,singlepop_t::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					    std::bind(insert_at_end<singlepop_t::gamete_t,singlepop_t::glist_t>,std::placeholders::_1,std::placeholders::_2),
					    std::bind(multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
					    std::bind(mutation_remover(),std::placeholders::_1,0,2*pop->N));
	remove_fixed_lost(&pop->mutations,&pop->fixations,&pop->fixation_times,&pop->mut_lookup,generation,2*pop->N);
	assert(check_sum(pop->gametes,2*pop->N));
      }
  }

  std::vector<int> sfs_from_sample(GSLrng_t * rng,const singlepop_t * pop,const unsigned & nsam)
  {
    map<double,unsigned> mutfreqs;
    unsigned twoN = 2*pop->N;
    for( unsigned i = 0 ; i < nsam ; ++i )
      {
	//pick a random chrom (w/replacement...)
	unsigned chrom = unsigned(gsl_ran_flat(rng->get(),0.,double(twoN)));
	//get pointer to that chrom from the individual
	auto gamete = (chrom%2==0.) ? pop->diploids[chrom/2].first : pop->diploids[chrom/2].second;
	//In this example, there are only neutral mutations, so that's what we'll iterate over
	for( auto m = gamete->mutations.begin() ; m != gamete->mutations.end() ; ++m )
	  {
	    auto pos_itr = mutfreqs.find( (*m)->pos );
	    if( pos_itr == mutfreqs.end() )
	      {
		mutfreqs.insert(std::make_pair((*m)->pos,1));
	      }
	    else
	      {
		pos_itr->second++;
	      }
	  }
      }
    //Now, fill in the SFS, omitting positions that are fixed in the sample
    std::vector<int> __rv(nsam-1,0u);
    for( const auto & __x : mutfreqs )
      {
	if (__x.second < nsam) __rv[__x.second-1]++;
      }
    return __rv;
  }
}
