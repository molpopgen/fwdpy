#include "metapop.hpp"
#include <exception>
#include <fwdpp/sugar/demography.hpp>
namespace fwdpy
{
    size_t
    migpop(const size_t &source_pop, const gsl_rng *r, const double &mig_prob)
    {
    }

    void
    re_init_mpop(metapop_t *mpop, const singlepop_t *pop)
    {
    }

    void
    copy_deme(metapop_t *mpop, const std::size_t i)
    {
        auto rv = KTfwd::copy_pop(*mpop, i);
        if (rv == -1)
            throw std::out_of_range("deme index out of range");
    }

    void
    remove_deme(metapop_t *mpop, const std::size_t i)
    {
        auto rv = KTfwd::remove_pop(*mpop, i);
        if (rv == -1)
            throw std::out_of_range("deme index out of range");
    }

    void
    merge_demes(metapop_t *mpop, const std::size_t i, const std::size_t j)
    {
        auto rv = KTfwd::merge_pops(*mpop, i, j);
        if (rv == -1)
            throw std::out_of_range("deme index out of range");
        if (rv == 1)
            throw std::runtime_error("cannot merge population into itself");
    }

    void
    split_deme(const gsl_rng *r, metapop_t *mpop, const std::size_t i,
               const unsigned N_new, const bool replacement)
    {
        auto rv = KTfwd::split_pop(r, *mpop, i, N_new, replacement);
        if (rv == -1)
            throw std::out_of_range("deme index out of range");
        if (rv == 1)
            throw std::runtime_error("N_new larger than parental deme size");
    }

    void
    admix_demes(const gsl_rng *r, metapop_t *mpop, const std::size_t i,
                const std::size_t j, const double prop_i, const unsigned N_new,
                const bool replacement)
    {
        auto rv
            = KTfwd::admix_pops(r, *mpop, i, j, prop_i, N_new, replacement);
        if (rv == -1)
            throw std::out_of_range("deme index out of range");
        if (rv == 1)
            throw std::runtime_error("invalid admixture proportion");
    }

    void
    swap_demes(metapop_t *mpop, const std::size_t i, const size_t j)
    {
        auto rv = KTfwd::swap_pops(*mpop, i, j);
        if (rv == -1)
            throw std::out_of_range("deme index out of range");
        // We "silently" ignore return values of 1
    }
}
