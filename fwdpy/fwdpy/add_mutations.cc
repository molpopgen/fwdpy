#include "fwdpy_add_mutations.hpp"
#include <fwdpp/sugar/add_mutation.hpp>
#include <stdexcept>

namespace fwdpy
{
    std::size_t
    add_mutation_cpp(singlepop_t *pop, const std::vector<std::size_t> &indlist,
                     const std::vector<short> &clist, const double pos,
                     const double s, const double h)
    {
        // If mutation position already exists, raise exception
        if (pop->mut_lookup.find(pos) != pop->mut_lookup.end())
            {
                throw std::runtime_error("attempt to create a mutation at "
                                         "position that already exists");
            }
        /*
          This function is a variadic template,
        */
        auto key = KTfwd::add_mutation(*pop, indlist, clist, pos, s, h,
                                       pop->generation);
        // update the lookup
        pop->mut_lookup.insert(pop->mutations[key].pos);
        return key;
    }
}
