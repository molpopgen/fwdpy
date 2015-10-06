#include <algorithm>
#include <memory>
#include <trajectories.hpp>

using namespace std;

namespace fwdpy
{
  map<string,vector<double> > get_singlepop_traj(const singlepop_t *pop,const unsigned minsojourn,const double minfreq)
  {
    vector<double> pos,freq,s;
    vector<double> generations;
    /*
      Key is origin, (pos,s), trajectories
    */
    //using trajtype = map< pair<unsigned,pair<double,double> >, vector<double> >;
    for( singlepop_t::trajtype::const_iterator itr = pop->trajectories.begin() ;
	 itr != pop->trajectories.end() ; ++itr )
      {
	double maxfreq = *max_element(itr->second.cbegin(),itr->second.cend());
	if( itr->second.size() >= minsojourn && maxfreq >= minfreq )
	  {
	    vector<unsigned> times(itr->second.size());
	    unsigned itime = itr->first.first;
	    generate(times.begin(),times.end(),[&itime]{ return itime++; });
	    generations.insert(generations.end(),times.begin(),times.end());
	    fill_n(back_inserter(pos),itr->second.size(),itr->first.second.first);
	    fill_n(back_inserter(s),itr->second.size(),itr->first.second.second);
	    copy(itr->second.begin(),itr->second.end(),back_inserter(freq));
	  }
      }
    map<string,vector<double>> rv;
    rv["pos"]=move(pos);
    rv["freq"]=move(freq);
    rv["generation"]=move(generations);
    rv["esize"]=move(s);
    return rv;
  }
}
