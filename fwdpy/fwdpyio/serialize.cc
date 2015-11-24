#include <fwdpyio/serialize.hpp>

using namespace std;

namespace fwdpy
{
  namespace serialize
  {
    void serialize_trajectories( const singlepop_t::trajtype & t, ostream & o )
    {
      size_t ntraj = t.size();
      o.write( reinterpret_cast<char*>(&ntraj),sizeof(size_t) );
      for( auto i = t.cbegin(); i != t.cend() ; ++i )
	{
	  o.write( reinterpret_cast<const char *>(&(i->first.first)),(sizeof(unsigned)));
	  o.write( reinterpret_cast<const char *>(&(i->first.second.first)),(sizeof(double)) );
	  o.write( reinterpret_cast<const char *>(&(i->first.second.second)),(sizeof(double)) );
	  ntraj = i->second.size();
	  o.write( reinterpret_cast<char*>(&ntraj),sizeof(unsigned) );
	  o.write( reinterpret_cast<const char *>(&(i->second[0])),ntraj*sizeof(double) );
	}
    }

    void deserialize_trajectories(istream & in, singlepop_t::trajtype * t )
    {
      size_t ntraj = t->size(),nfreqs;
      unsigned a;
      double rest[2];
      in.read( reinterpret_cast<char*>(&ntraj),sizeof(size_t) );
      for(unsigned i=0;i<ntraj;++i)
	{
	  in.read(reinterpret_cast<char*>(&a),sizeof(unsigned));
	  in.read(reinterpret_cast<char*>(&rest[0]),2*sizeof(double));
	  in.read(reinterpret_cast<char*>(&nfreqs),sizeof(unsigned));
	  vector<double> temp(nfreqs);
	  if(nfreqs)
	    {
	      in.read(reinterpret_cast<char*>(&temp[0]),nfreqs*sizeof(double));
	    }
	  t->insert( std::make_pair( std::make_pair(a,std::make_pair(rest[0],rest[1])), std::move(temp) ) );
	}
    }
    
    string serialize_singlepop(const fwdpy::singlepop_t * pop)
    {
      KTfwd::serialize rv;
      rv.buffer.write(reinterpret_cast<const char *>(&(pop->generation)),sizeof(unsigned));
      serialize_trajectories(pop->trajectories,rv.buffer);
      rv(*pop,KTfwd::mutation_writer(),fwdpy::diploid_writer());
      return rv.buffer.str();
    }

    vector<shared_ptr<singlepop_t> > deserialize_singlepop(const vector<string> & strings)
    {
      vector<shared_ptr<singlepop_t> > rv;
      for(unsigned i=0;i<strings.size();++i)
	{
	  KTfwd::serialize st;
	  st.buffer.str(strings[i]);
	  st.buffer.seekg(0);
	  singlepop_t pop(0);
	  st.buffer.read(reinterpret_cast<char*>(&pop.generation),sizeof(unsigned));
	  deserialize_trajectories(st.buffer,&pop.trajectories);
	  KTfwd::deserialize d;
	  d(pop,st,KTfwd::mutation_reader<KTfwd::popgenmut>(),fwdpy::diploid_reader());
	  rv.emplace_back( shared_ptr<singlepop_t>(new singlepop_t(move(pop))) );
	}
      return rv;
    }

    string serialize_metapop(const fwdpy::metapop_t * pop)
    {
      KTfwd::serialize rv;
      rv.buffer.write(reinterpret_cast<const char *>(&(pop->generation)),sizeof(unsigned));
      rv(*pop,KTfwd::mutation_writer(),fwdpy::diploid_writer());
      return rv.buffer.str();
    }

    vector<shared_ptr<metapop_t> > deserialize_metapop(const vector<string> & strings)
    {
      vector<shared_ptr<metapop_t> > rv;
      for(unsigned i=0;i<strings.size();++i)
	{
	  KTfwd::serialize st;
	  st.buffer.str(strings[i]);
	  st.buffer.seekg(0);
	  metapop_t pop({0});
	  st.buffer.read(reinterpret_cast<char*>(&pop.generation),sizeof(unsigned));
	  KTfwd::deserialize d;
	  d(pop,st,KTfwd::mutation_reader<KTfwd::popgenmut>(),fwdpy::diploid_reader());
	  rv.emplace_back( shared_ptr<metapop_t>(new metapop_t(move(pop))) );
	}
      return rv;
    }
  }
}

