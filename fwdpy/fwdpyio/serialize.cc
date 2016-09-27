#include "fwdpyio_serialize.hpp"

using namespace std;

namespace fwdpy
{
    namespace serialize
    {

        string
        serialize_singlepop(const fwdpy::singlepop_t *pop)
        {
			return pop->serialize();
        }

        string
        serialize_metapop(const fwdpy::metapop_t *pop)
        {
			return pop->serialize();
        }

        string
        serialize_multilocus(const fwdpy::multilocus_t *pop)
        {
			return pop->serialize();
        }

        vector<shared_ptr<singlepop_t>>
        deserialize_singlepop(const vector<string> &strings)
        {
            return deserialize_details<singlepop_t>()(strings, 0u);
        }

        vector<shared_ptr<metapop_t>>
        deserialize_metapop(const vector<string> &strings)
        {
            return deserialize_details<metapop_t>()(strings, std::vector<unsigned>(0u));
        }

        vector<shared_ptr<multilocus_t>>
        deserialize_multilocus(const vector<string> &strings)
        {
            return deserialize_details<multilocus_t>()(strings, 0u, 0u);
        }
    }
}
