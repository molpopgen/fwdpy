#ifndef FWDPY_GET_SELECTED_MUT_DATA_HPP
#define FWDPY_GET_SELECTED_MUT_DATA_HPP
#include "sampler_base.hpp"
#include "types.hpp"
#include <limits>
#include <unordered_map>
namespace fwdpy
{
    class selected_mut_tracker : public sampler_base
    /*!
      \brief A "sampler" for recording frequency trajectories of selected
      mutations.
      \ingroup samplers
    */
    {
      public:
        /*!
          \brief Internal representation of mutation frequencies during a
          simulation

          \note Used in fwdpy::selected_mut_tracker
        */
        using trajVec = std::vector<std::pair<unsigned, double>>;
        using posEsize = std::pair<double, double>;
        using innerMap = std::map<posEsize, trajVec>;
        using final_t = std::unordered_map<KTfwd::uint_t, innerMap>;

        virtual void
        operator()(const singlepop_t *pop, const unsigned generation)
        {
            call_operator_details(pop, generation);
        }
        virtual void
        operator()(const multilocus_t *pop, const unsigned generation)
        {
            call_operator_details(pop, generation);
        }

        final_t
        final() const
        {
            return trajectories;
        }

        explicit selected_mut_tracker() noexcept : trajectories(final_t())
        {
            trajectories.reserve(1000000);
        }

      private:
        final_t trajectories;
        template <typename pop_t>
        inline void
        call_operator_details(const pop_t *pop, const unsigned generation)
        {
            for (std::size_t i = 0; i < pop->mcounts.size(); ++i)
                {
                    if (pop->mcounts[i])
                        { // if mutation is not extinct
                            const auto &__m = pop->mutations[i];
                            if (!__m.neutral)
                                {
                                    const auto freq
                                        = double(pop->mcounts[i])
                                          / double(2 * pop->diploids.size());
                                    auto __itr = trajectories.find(__m.g);
                                    if (__itr == trajectories.end())
                                        {
                                            final_t::mapped_type m;
                                            trajVec v(1, { generation, freq });
                                            v.reserve(50);
                                            m[{ __m.pos, __m.s }]
                                                = std::move(v);
                                            trajectories.emplace(__m.g,
                                                                 std::move(m));
                                        }
                                    else
                                        {
                                            auto __ps = __itr->second.find(
                                                { __m.pos, __m.s });
                                            if (__ps == __itr->second.end())
                                                {
                                                    // update the data
                                                    __itr->second.emplace(
                                                        std::make_pair(__m.pos,
                                                                       __m.s),
                                                        std::
                                                            vector<std::
                                                                       pair<unsigned,
                                                                            double>>(
                                                                1,
                                                                std::make_pair(
                                                                    generation,
                                                                    freq)));
                                                }
                                            else
                                                {
                                                    // Don't keep updating for
                                                    // fixed
                                                    // variants
                                                    if (__ps->second.back()
                                                            .second
                                                        < 1.)
                                                        {
                                                            __ps->second
                                                                .emplace_back(
                                                                    generation,
                                                                    freq);
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
    };
}

#endif
