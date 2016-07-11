/*!
  \file fwdpp_features.hpp

  fwdpy serves as a test-bed for fwdpp.  Sometimes, I realize that fwdpp
  is missing a feature.  Those features may first appear here before getting
  moved over to fwdpp.
*/

namespace fwdpy
{
  /*!
    Label all fixed neutral variant and all extinct variants for recycling. Copy fixations and fixation times
    for neutral mutations into containers.

    Fixed, non-neutral variants get copied into fixations and fixation times, so that fixation times
    can get records.

    This function is identical in name and interface to the current fwdpp function in fwdpp/util.hpp

    \note: lookup must be compatible with lookup->erase(lookup->find(double))
  */
  template<typename mcont_t,
	   typename fixation_container_t,
	   typename fixation_time_container_t,
	   typename mutation_lookup_table>
  void update_mutations_n( mcont_t & mutations,
			   fixation_container_t & fixations,
			   fixation_time_container_t & fixation_times,
			   mutation_lookup_table & lookup,
			   std::vector<uint_t> & mcounts,
			   const unsigned & generation,const unsigned & twoN )
  {
    static_assert( typename traits::is_mutation_t<typename mcont_t::value_type>::type(),
		   "mutation_type must be derived from KTfwd::mutation_base" );
    assert(mcounts.size()==mutations.size());
    for(unsigned i=0;i<mcounts.size();++i)
      {
	assert(mcounts[i] <= twoN);
	if(mcounts[i]==twoN)
	  {
	    if(mutations[i].neutral)
	      {
		fixations.push_back(mutations[i]);
		fixation_times.push_back(generation);
		mcounts[i]=0; //set count to zero to mark mutation as "recyclable"
		lookup.erase(mutations[i].pos);
	      }
	    else
	      {
		fixations.push_back(mutations[i]);
		fixation_times.push_back(generation);
	      }
	  }
	if(!mcounts[i]) lookup.erase(mutations[i].pos);
      }
  }
}
