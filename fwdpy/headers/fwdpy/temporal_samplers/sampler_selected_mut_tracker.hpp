#ifndef FWDPY_GET_SELECTED_MUT_DATA_HPP
#define FWDPY_GET_SELECTED_MUT_DATA_HPP
#include "fwdpy/temporal_samplers/sampler_base.hpp"
#include "fwdpy/types.hpp"
#include <mutex>
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

        final_t::const_iterator
        begin() const
        {
            return this->trajectories.cbegin();
        }
        final_t::const_iterator
        end() const
        {
            return this->trajectories.cend();
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

    using origin_filter_fxn = bool (*)(const unsigned);
    using pos_esize_filter_fxn = bool (*)(const std::pair<double, double> &);
    using freq_filter_fxn
        = bool (*)(const std::vector<std::pair<KTfwd::uint_t, double>> &);
    bool all_origins_pass(const unsigned);
    bool all_pos_esize_pass(const std::pair<double, double> &);
    bool all_freqs_pass(const std::vector<std::pair<KTfwd::uint_t, double>> &);
    struct trajFilter
    {
        origin_filter_fxn origin_filter;
        pos_esize_filter_fxn pos_esize_filter;
        freq_filter_fxn freq_filter;
        trajFilter()
            : origin_filter(&all_origins_pass),
              pos_esize_filter(&all_pos_esize_pass),
              freq_filter(&all_freqs_pass)
        {
        }
        virtual bool
        apply_origin_filter(const unsigned origin) const
        {
            return origin_filter(origin);
        }
        virtual bool
        apply_pos_esize_filter(const std::pair<double, double> &pe) const
        {
            return pos_esize_filter(pe);
        }
        virtual bool
        apply_freq_filter(
            const std::vector<std::pair<unsigned, double>> &freqs) const
        {
            return freq_filter(freqs);
        }
    };

    template <typename T> class trajFilterData : public trajFilter
    {
      public:
        using origin_filter_fxn_T = bool (*)(const unsigned, const T &);
        using pos_esize_filter_fxn_T
            = bool (*)(const std::pair<double, double> &, const T &);
        using freq_filter_fxn_T
            = bool (*)(const std::vector<std::pair<KTfwd::uint_t, double>> &,
                       const T &);

      private:
        T data;
        origin_filter_fxn_T origin_filter;
        pos_esize_filter_fxn_T pos_esize_filter;
        freq_filter_fxn_T freq_filter;

      public:
        trajFilterData(const T &data_)
            : data(data_), origin_filter(nullptr), pos_esize_filter(nullptr),
              freq_filter(nullptr)
        {
        }
        void
        register_callback(origin_filter_fxn_T o)
        {
            origin_filter = o;
        }
        void
        register_callback(pos_esize_filter_fxn_T p)
        {
            pos_esize_filter = p;
        }
        void
        register_callback(freq_filter_fxn_T f)
        {
            freq_filter = f;
        }
        bool
        apply_origin_filter(const unsigned origin) const final
        {
            if (origin_filter == nullptr)
                {
                    return trajFilter::apply_origin_filter(origin);
                }
            return origin_filter(origin, data);
        }
        bool
        apply_pos_esize_filter(const std::pair<double, double> &pe) const final
        {
            if (pos_esize_filter == nullptr)
                {
                    return trajFilter::apply_pos_esize_filter(pe);
                }
            return pos_esize_filter(pe,data);
        }
        bool
        apply_freq_filter(
            const std::vector<std::pair<unsigned, double>> &freqs) const final
        {
            if (freq_filter == nullptr)
                {
                    return trajFilter::apply_freq_filter(freqs);
                }
            return freq_filter(freqs, data);
        }
    };
    void
    traj2sql(const std::vector<std::unique_ptr<fwdpy::sampler_base>> &samplers,
             const std::shared_ptr<std::mutex> &dblock, const trajFilter *tf,
             const std::string &dbname, unsigned threshold,
             const unsigned label, const bool onedb, const bool append);
}

#endif
