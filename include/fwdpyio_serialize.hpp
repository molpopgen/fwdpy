#ifndef __FWDPY_SERIALIZE_HPP__
#define __FWDPY_SERIALIZE_HPP__

#include "serialization_common.hpp"
#include "types.hpp"
#include <fwdpp/sugar/serialization.hpp>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace fwdpy
{
    namespace serialize
    {

        template <typename poptype> struct deserialize_details
        {
            template <typename... constructor_data>
            std::vector<std::shared_ptr<poptype>>
            operator()(const std::vector<std::string> &strings,
                       constructor_data... cdata)
            {
                std::vector<std::shared_ptr<poptype>> rv;
                for (unsigned i = 0; i < strings.size(); ++i)
                    {
                        std::istringstream buffer;
                        buffer.str(strings[i]);
                        buffer.seekg(0);
                        poptype pop(cdata...);
						pop.deserialize(buffer.str());
                        rv.emplace_back(std::shared_ptr<poptype>(
                            new poptype(std::move(pop))));
                    }
                return rv;
            }
        };

        std::string serialize_singlepop(const singlepop_t *spop);
        std::vector<std::shared_ptr<singlepop_t>>
        deserialize_singlepop(const std::vector<std::string> &strings);
        std::string serialize_metapop(const fwdpy::metapop_t *pop);
        std::vector<std::shared_ptr<metapop_t>>
        deserialize_metapop(const std::vector<std::string> &strings);
        std::string serialize_multilocus(const fwdpy::multilocus_t *pop);
        std::vector<std::shared_ptr<multilocus_t>>
        deserialize_multilocus(const std::vector<std::string> &strings);
    }
}

#endif
