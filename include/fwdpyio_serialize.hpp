#ifndef __FWDPY_SERIALIZE_HPP__
#define __FWDPY_SERIALIZE_HPP__

#include <memory>
#include <vector>
#include <string>
#include <fwdpp/sugar/serialization.hpp>
#include "types.hpp"

namespace fwdpy {
namespace serialize {

template<typename poptype>
std::string serialize_details(const poptype * pop) {
    KTfwd::serialize rv;
    rv.buffer.write(reinterpret_cast<const char *>(&(pop->generation)),sizeof(unsigned));
    rv(*pop,KTfwd::mutation_writer(),fwdpy::diploid_writer());
    return rv.buffer.str();
}

template<typename poptype>
struct deserialize_details {
    template<typename ...constructor_data>
    std::vector<std::shared_ptr<poptype> > operator()(const std::vector<std::string> & strings,
            constructor_data... cdata) {
        std::vector<std::shared_ptr<poptype> > rv;
        for(unsigned i=0; i<strings.size(); ++i) {
            KTfwd::serialize st;
            st.buffer.str(strings[i]);
            st.buffer.seekg(0);
            poptype pop(cdata...);
            st.buffer.read(reinterpret_cast<char*>(&pop.generation),sizeof(unsigned));
            KTfwd::deserialize d;
            d(pop,st,KTfwd::mutation_reader<KTfwd::popgenmut>(),fwdpy::diploid_reader());
            rv.emplace_back(std::shared_ptr<poptype>(new poptype(std::move(pop))));
        }
        return rv;
    }
};

std::string serialize_singlepop(const singlepop_t * spop);
std::vector<std::shared_ptr<singlepop_t> > deserialize_singlepop(const std::vector<std::string> & strings);
std::string serialize_metapop(const fwdpy::metapop_t * pop);
std::vector<std::shared_ptr<metapop_t> > deserialize_metapop(const std::vector<std::string> & strings);
std::string serialize_multilocus(const fwdpy::multilocus_t * pop);
std::vector<std::shared_ptr<multilocus_t> > deserialize_multilocus(const std::vector<std::string> & strings);
}
}


#endif
