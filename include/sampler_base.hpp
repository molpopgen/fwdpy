#ifndef FWDPY_SAMPLER_BASE_HPP
#define FWDPY_SAMPLER_BASE_HPP

#include "types.hpp"
#include <stdexcept>
#include <thread>
#include <vector>

namespace fwdpy
{
    struct sampler_base
    /*!
      Base class for a temporal sampler.

      This class cannot be used as a sampler. Attempting to build a sampler
      around this will result in std::runtime_error being thrown from
      operator()
      during simulations.
     */
    {
        virtual void
        operator()(const singlepop_t *, const unsigned)
        {
            throw std::runtime_error(
                "sampler type not implemented for single deme simulations");
        };
        virtual void
        operator()(const multilocus_t *, const unsigned)
        {
            throw std::runtime_error(
                "sampler type not implemented for multi-locus simulations");
        };
        virtual void
        operator()(const metapop_t *, const unsigned)
        {
            throw std::runtime_error(
                "sampler type not implemented for metapopulation simulations");
        };
        virtual void
        cleanup()
        /*!
          Evolve functions will call this before returning.
          This function gives samplers a chance to do things like clean
          up any RAM that they may have allocated during the course of
          a simulation.

          \note: the stuff to return to Python should not be cleaned up by this
          function!
         */
        {
        }
        virtual ~sampler_base() {}
    };

    inline void
    clear_samplers(std::vector<std::unique_ptr<sampler_base>> &v)
    {
        v.clear();
        std::vector<std::unique_ptr<sampler_base>> t;
        v = std::move(t);
    }

    template <typename pop_t>
    inline void
    apply_sampler_wrapper(sampler_base *s, const pop_t *pop)
    /*!
      Thin wrapper function for apply_sampler_cpp.  This funcion is needed to
      avoid slicing the pointer to a sampler down to a pointer to the base
      class.
     */
    {
        s->operator()(pop, pop->generation);
    }

    template <typename T>
    inline void
    apply_sampler_cpp(
        const std::vector<std::shared_ptr<T>> &popvec,
        const std::vector<std::unique_ptr<sampler_base>> &samplers)
    /*!
      Apply the i-th ampler to the i-th pop using std::thread.

      Throws runtime_error if popvec.size()!=samplers.size().
     */
    {
        if (popvec.size() != samplers.size())
            throw std::runtime_error("Containers of populations and samplers "
                                     "must be equal in length");
        std::vector<std::thread> threads;
        for (std::size_t i = 0; i < popvec.size(); ++i)
            {
                threads.emplace_back(apply_sampler_wrapper<T>,
                                     samplers[i].get(), popvec[i].get());
            }
        for (auto &t : threads)
            t.join();
    }

    template <typename T>
    inline void
    apply_sampler_single_cpp(
        const T *pop,
        const std::vector<std::unique_ptr<sampler_base>> &samplers)
    {
        apply_sampler_wrapper(samplers[0].get(), pop);
    }

    template <typename final_t> struct custom_sampler : public sampler_base
    /*!
      The basis of a custom sampler that requires no extra data
      in order to implement its call operator.
    */
    {
        using singlepop_call_operator
            = void (*)(const singlepop_t *, const unsigned, final_t &);
        using multilocus_call_operator
            = void (*)(const multilocus_t *, const unsigned, final_t &);
        using metapop_call_operator
            = void (*)(const metapop_t *, const unsigned, final_t &);
        using cleanup_fxn = void (*)(final_t &);

        final_t f;
        singlepop_call_operator scall;
        multilocus_call_operator mcall;
        metapop_call_operator metacall;
        cleanup_fxn cleanupcall;

        custom_sampler(singlepop_call_operator s)
            : f(final_t()), scall(s), mcall(nullptr), metacall(nullptr),
              cleanupcall(nullptr)

        {
        }

        custom_sampler(multilocus_call_operator s)
            : f(final_t()), scall(nullptr), mcall(s), metacall(nullptr),
              cleanupcall(nullptr)
        {
        }

        custom_sampler(metapop_call_operator s)
            : f(final_t()), scall(nullptr), mcall(nullptr), metacall(s),
              cleanupcall(nullptr)
        {
        }

        void
        operator()(const singlepop_t *pop, const unsigned generation)
        {
            if (scall != nullptr)
                {
                    scall(pop, generation, f);
                }
            else
                sampler_base::operator()(pop, generation);
        }

        void
        operator()(const multilocus_t *pop, const unsigned generation)
        {
            if (mcall != nullptr)
                {
                    mcall(pop, generation, f);
                }
            else
                sampler_base::operator()(pop, generation);
        }

        void
        operator()(const metapop_t *pop, const unsigned generation)
        {
            if (metacall != nullptr)
                {
                    metacall(pop, generation, f);
                }
            else
                sampler_base::operator()(pop, generation);
        }

        void
        register_callback(singlepop_call_operator c)
        {
            scall = c;
        }

        void
        register_callback(multilocus_call_operator c)
        {
            mcall = c;
        }

        void
        register_callback(metapop_call_operator c)
        {
            metacall = c;
        }

        void
        register_callback(cleanup_fxn c)
        {
            cleanupcall = c;
        }

        void
        cleanup()
        {
            if (cleanupcall != nullptr)
                {
                    cleanupcall(f);
                }
            else
                sampler_base::cleanup();
        }

        final_t
        final()
        {
            return f;
        }
    };

    template <typename final_t, typename data_t>
    struct custom_sampler_data : public sampler_base
    /*!
      The basis of a custom sampler requires extra data
      in order to implement its call operator.
    */
    {
        using singlepop_call_operator
            = void (*)(const singlepop_t *, const unsigned, final_t &,
                       data_t &);
        using multilocus_call_operator
            = void (*)(const multilocus_t *, const unsigned, final_t &,
                       data_t &);
        using metapop_call_operator
            = void (*)(const metapop_t *, const unsigned, final_t &, data_t &);
        using cleanup_fxn = void (*)(final_t &, data_t &);
        final_t f;
        data_t data;

        singlepop_call_operator scall;
        multilocus_call_operator mcall;
        metapop_call_operator metacall;
        cleanup_fxn cleanupcall;

        custom_sampler_data(singlepop_call_operator s, const data_t &d)
            : f(final_t()), data(d), scall(s), mcall(nullptr),
              metacall(nullptr), cleanupcall(nullptr)
        {
        }

        custom_sampler_data(multilocus_call_operator s, const data_t &d)
            : f(final_t()), data(d), scall(nullptr), mcall(s),
              metacall(nullptr), cleanupcall(nullptr)
        {
        }

        custom_sampler_data(metapop_call_operator s, const data_t &d)
            : f(final_t()), data(d), scall(nullptr), mcall(nullptr),
              metacall(s), cleanupcall(nullptr)
        {
        }

        void
        operator()(const singlepop_t *pop, const unsigned generation)
        {
            if (scall != nullptr)
                {
                    scall(pop, generation, f, data);
                }
            else
                sampler_base::operator()(pop, generation);
        }

        void
        operator()(const multilocus_t *pop, const unsigned generation)
        {
            if (mcall != nullptr)
                {
                    mcall(pop, generation, f, data);
                }
            else
                sampler_base::operator()(pop, generation);
        }

        void
        operator()(const metapop_t *pop, const unsigned generation)
        {
            if (metacall != nullptr)
                {
                    metacall(pop, generation, f, data);
                }
            else
                sampler_base::operator()(pop, generation);
        }

        void
        register_callback(singlepop_call_operator c)
        {
            scall = c;
        }

        void
        register_callback(multilocus_call_operator c)
        {
            mcall = c;
        }

        void
        register_callback(metapop_call_operator c)
        {
            metacall = c;
        }

        void
        register_callback(cleanup_fxn c)
        {
            cleanupcall = c;
        }

        void
        cleanup()
        {
            if (cleanupcall != nullptr)
                {
                    cleanupcall(f, data);
                }
            else
                sampler_base::cleanup();
        }

        final_t
        final()
        {
            return f;
        }
    };
}

#endif
