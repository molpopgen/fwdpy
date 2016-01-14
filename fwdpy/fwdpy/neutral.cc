#include <neutral.hpp>
#include <map>
#include <thread>
using namespace std;
using namespace KTfwd;

namespace fwdpy {

  void evolve_pop_details(gsl_rng * rng,
			  singlepop_t * pop,
			  const std::vector<unsigned> & nlist,
			  const double & theta,
			  const double & rho)
  {
    double mu = theta/(4.*double(pop->N)),littler=rho/(4.*double(pop->N));

    std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng);

    for( unsigned generation = 0; generation < nlist.size() ; ++generation,++pop->generation )
      {
	//Iterate the population through 1 generation
	double wbar = sample_diploid(rng,
				     pop->gametes,  //non-const pointer to gametes
				     pop->diploids, //non-const pointer to diploids
				     pop->mutations, //non-const pointer to mutations
				     pop->mcounts,
				     pop->N,     //current pop size
				     nlist[generation], //next pop size,
				     mu,    //mutation rate per gamete
				     std::bind(infsites(),std::placeholders::_1,std::placeholders::_2,
					       rng,std::ref(pop->mut_lookup),generation,
					       mu,0.,[&rng](){return gsl_rng_uniform(rng);},[](){return 0.;},[](){return 0.;}),
				     //The recombination policy includes the uniform crossover rate
				     std::bind(KTfwd::poisson_xover(),rng,littler,0.,1.,
					       std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
				     std::bind(multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,2.),
				     pop->neutral,pop->selected);
	//update pop size in data structure
	pop->N = nlist[generation];
	update_mutations(pop->mutations,pop->fixations,pop->fixation_times,pop->mut_lookup,pop->mcounts,generation,2*pop->N);
	assert(check_sum(pop->gametes,2*pop->N));
      }
  }

  void evolve_pop(GSLrng_t * rng, std::vector<std::shared_ptr<singlepop_t> > * pops,const std::vector<unsigned> & nlist,const double & theta, const double & rho)
  {
    vector<thread> threads;
    std::vector<GSLrng_t> rngs;
    for( unsigned i = 0 ; i < pops->size() ; ++i )
      {
	rngs.emplace_back(GSLrng_t(gsl_rng_get(rng->get())));
      }
    for( unsigned i = 0 ; i < pops->size() ; ++i )
      {
	threads.push_back( thread(evolve_pop_details,rngs[i].get(),pops->operator[](i).get(),nlist,theta,rho) );
      }
    for(unsigned i=0;i<threads.size();++i) threads[i].join();
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
	//auto gamete = (chrom%2==0.) ? pop->diploids[chrom/2].first : pop->diploids[chrom/2].second;
	const auto & gamete = (chrom%2==0.) ? pop->gametes[pop->diploids[chrom/2].first] :
	  pop->gametes[pop->diploids[chrom/2].second];
	//In this example, there are only neutral mutations, so that's what we'll iterate over
	//for( auto m = gamete->mutations.begin() ; m != gamete->mutations.end() ; ++m )
	for(const auto & m : gamete.mutations )
	  {
	    auto pos_itr = mutfreqs.find( pop->mutations[m].pos );
	    if( pos_itr == mutfreqs.end() )
	      {
		mutfreqs.insert(std::make_pair(pop->mutations[m].pos,1));
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
