#include "fwdpyio_serialize.hpp"

using namespace std;

namespace fwdpy {
namespace serialize {

template<typename poptype>
string serialize_details(const poptype * pop) {
    KTfwd::serialize rv;
    rv.buffer.write(reinterpret_cast<const char *>(&(pop->generation)),sizeof(unsigned));
    rv(*pop,KTfwd::mutation_writer(),fwdpy::diploid_writer());
    return rv.buffer.str();
}

string serialize_singlepop(const fwdpy::singlepop_t * pop) {
    return serialize_details(pop);
}

string serialize_metapop(const fwdpy::metapop_t * pop) {
    return serialize_details(pop);
}

string serialize_multilocus(const fwdpy::multilocus_t * pop) {
    return serialize_details(pop);
}

template<typename poptype>
struct deserialize_details {
    template<typename ...constructor_data>
    vector<shared_ptr<poptype> > operator()(const vector<string> & strings,
                                            constructor_data... cdata) {
        vector<shared_ptr<poptype> > rv;
        for(unsigned i=0; i<strings.size(); ++i) {
            KTfwd::serialize st;
            st.buffer.str(strings[i]);
            st.buffer.seekg(0);
            poptype pop(cdata...);
            st.buffer.read(reinterpret_cast<char*>(&pop.generation),sizeof(unsigned));
            KTfwd::deserialize d;
            d(pop,st,KTfwd::mutation_reader<KTfwd::popgenmut>(),fwdpy::diploid_reader());
            rv.emplace_back(shared_ptr<poptype>(new poptype(move(pop))));
        }
        return rv;
    }
};

vector<shared_ptr<singlepop_t> > deserialize_singlepop(const vector<string> & strings) {
    return deserialize_details<singlepop_t>()(strings,0u);
}

vector<shared_ptr<metapop_t> > deserialize_metapop(const vector<string> & strings) {
    return deserialize_details<metapop_t>()(strings,0u);
}

vector<shared_ptr<multilocus_t> > deserialize_multilocus(const vector<string> & strings) {
    return deserialize_details<multilocus_t>()(strings,0u,0u);
}
}
}
