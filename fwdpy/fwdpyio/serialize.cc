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
	  o.write( reinterpret_cast<const char*>(&(std::get<static_cast<std::size_t>(traj_key_values::deme) >(i->first))),
		   sizeof(decltype(std::get<static_cast<std::size_t>(traj_key_values::deme) >(i->first))) );
	  o.write( reinterpret_cast<const char*>(&(std::get<static_cast<std::size_t>(traj_key_values::origin) >(i->first))),
		   sizeof(decltype(std::get<static_cast<std::size_t>(traj_key_values::origin) >(i->first))) );
	  o.write( reinterpret_cast<const char*>(&(std::get<static_cast<std::size_t>(traj_key_values::pos) >(i->first))),
		   sizeof(decltype(std::get<static_cast<std::size_t>(traj_key_values::pos) >(i->first))) );
	  o.write( reinterpret_cast<const char*>(&(std::get<static_cast<std::size_t>(traj_key_values::esize) >(i->first))),
		   sizeof(decltype(std::get<static_cast<std::size_t>(traj_key_values::esize) >(i->first))) );
	  ntraj = i->second.size();
	  o.write( reinterpret_cast<char*>(&ntraj),sizeof(unsigned) );
	  o.write( reinterpret_cast<const char *>(i->second.data()),ntraj*sizeof(double) );
	}
    }

    void deserialize_trajectories(istream & in, singlepop_t::trajtype * t )
    {
      using qvec_t = singlepop_t::trajtype::value_type::second_type;
      std::size_t ntraj,nfreqs;
      in.read( reinterpret_cast<char*>(&ntraj),sizeof(decltype(ntraj)) );

      unsigned d,o;
      double q,s;
      for( unsigned i=0;i<ntraj;++i )
	{
	  in.read( reinterpret_cast<char*>(&d),sizeof(decltype(d)) );
	  in.read( reinterpret_cast<char*>(&o),sizeof(decltype(o)) );
	  in.read( reinterpret_cast<char*>(&q),sizeof(decltype(q)) );
	  in.read( reinterpret_cast<char*>(&s),sizeof(decltype(s)) );
	  auto key = std::make_tuple(d,o,q,s);
	  in.read( reinterpret_cast<char*>(&nfreqs),sizeof(decltype(nfreqs)));
	  qvec_t qvec(nfreqs);
	  in.read( reinterpret_cast<char*>(qvec.data()),nfreqs*sizeof(qvec_t::value_type));
	  t->emplace( std::move(key),std::move(qvec) );
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

