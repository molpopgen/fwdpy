/*! \file fwdpy_serialization.hpp
 * \brief Helper functions for object-level serialization
 */
#ifndef FWDPY_SERIALIZATION_HPP
#define FWDPY_SERIALIZATION_HPP
#include <fwdpp/sugar/serialization.hpp>

namespace fwdpy {
namespace serialize {

template<typename poptype,typename mwriter_t,typename dipwriter_t>
int gzserialize_details(const poptype & pop,
                        const mwriter_t & mwriter,
                        const dipwriter_t & dipwriter,
                        const char * filename, bool append) {
    gzFile f;
    if(append) {
        f=gzopen(filename,"ab");
    } else {
        f=gzopen(filename,"wb");
    }
    auto rv = gzwrite(f,reinterpret_cast<const char*>(&pop.generation),
                      sizeof(decltype(pop.generation)));
    KTfwd::gzserialize s;
    rv += s(f,pop,mwriter,dipwriter);
    gzclose(f);
    return rv;
}
}
}
#endif
