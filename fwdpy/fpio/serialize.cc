#include <serialize.hpp>
#include <iostream>
namespace fwdpy
{
  namespace serialize
  {
    std::string serialize_singlepop(const fwdpy::singlepop_t * pop)
    {
      KTfwd::serialize rv;
      rv(*pop,KTfwd::mutation_writer(),fwdpy::diploid_writer());
      return rv.buffer.str();
    }

    void deserialize_singlepop(const std::vector<std::string> & strings,
			       std::vector<std::shared_ptr<singlepop_t> > * pops)
    {
      for(unsigned i=0;i<strings.size();++i)
	{
	  KTfwd::serialize st;
	  st.buffer.str(strings[i]);
	  st.buffer.seekg(0);
	  //std::cerr << st.buffer.str() << '\n';
	  KTfwd::deserialize t;
	  t(*pops->operator[](i).get(),st,KTfwd::mutation_reader<KTfwd::popgenmut>(),fwdpy::diploid_reader());
	  std::cerr << pops->operator[](i).get()->mutations.size() << ' '
		    << pops->operator[](i).get()->diploids.size() << ' '
		    << pops->operator[](i).get()->gametes.size() << '\n';
	}
    }
  }
}

