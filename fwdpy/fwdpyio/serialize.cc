#include "fwdpyio_serialize.hpp"

using namespace std;

namespace fwdpy
{
  namespace serialize
  {
    string serialize_singlepop(const fwdpy::singlepop_t * pop)
    {
      KTfwd::serialize rv;
      rv.buffer.write(reinterpret_cast<const char *>(&(pop->generation)),sizeof(unsigned));
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

