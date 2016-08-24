#ifndef FWDPY_GET_SELECTED_MUT_DATA_HPP
#define FWDPY_GET_SELECTED_MUT_DATA_HPP

#include <limits>
#include "types.hpp"
#include "sampler_base.hpp"
namespace fwdpy {
struct selected_mut_data {
    double pos,esize;
    unsigned origin;
    using label_t =decltype(KTfwd::mutation_base::xtra);
    label_t label;
    selected_mut_data(unsigned g,double p,double e,label_t l) :
        pos(p),esize(e),origin(g),label(l) {
    }
    selected_mut_data() : pos(std::numeric_limits<double>::quiet_NaN()),
        esize(std::numeric_limits<double>::quiet_NaN()),
        origin(std::numeric_limits<unsigned>::max()),
        label(std::numeric_limits<label_t>::max())
        /*!
          This constructor assigns NaN or "max_int"
          values to members.
         */
    {
    }
    inline bool operator==(const selected_mut_data & rhs) const noexcept {
        return this->origin == rhs.origin &&
               this->pos == rhs.pos &&
               this->esize == rhs.esize
               && this->label == rhs.label;

    }
	inline bool operator<(const selected_mut_data & rhs) const noexcept {
		if(this->origin < rhs.origin) return true;
		if(rhs.origin < this->origin) return false;
		if(this->pos < rhs.pos) return true;
		if(rhs.pos < this->pos) return false;
		if(this->esize < rhs.esize) return true;
		if(rhs.esize < this->esize) return false;
		if(this->label < rhs.label) return true;
		if(rhs.label < this->label) return false;
		return false;
	}
};

struct selected_mut_data_tidy
/*!
  \brief Convenience type for conversion to dict -> pandas.DataFrame
 */
{
    double pos,esize,freq;
    unsigned origin,generation;
    using label_t=decltype(KTfwd::mutation_base::xtra);
    label_t label;
    selected_mut_data_tidy(unsigned o,unsigned g,double p,double q,double e,label_t l) :
        pos(p),esize(e),freq(q),origin(o),generation(g),label(l) {
    }
};

//non-inline!  This is part of fwdpy's main module.
std::vector<selected_mut_data_tidy> tidy_trajectory_info( const std::vector<std::pair<selected_mut_data,std::vector<std::pair<unsigned,double>>>> & trajectories,
        const unsigned min_sojourn ,
        const double min_freq);

//! Used internally to convert C++11 types to something Cython will understand
enum class traj_key_values : std::size_t { deme,origin,pos,esize,label };

/*!
  \brief Internal representation of mutation frequencies during a simulation

  \note Used in fwdpy::selected_mut_tracker
*/
using trajectories_t = std::map<selected_mut_data,std::vector<std::pair<unsigned,double> >>;

class selected_mut_tracker : public sampler_base
/*!
  \brief A "sampler" for recording frequency trajectories of selected mutations.
  \ingroup samplers
*/
{
  public:
    using final_t = std::vector< std::pair<selected_mut_data, std::vector<std::pair<unsigned,double> > > >;

    virtual void operator()(const singlepop_t * pop, const unsigned generation) {
        call_operator_details(pop,generation);
    }
    virtual void operator()(const multilocus_t * pop, const unsigned generation) {
        call_operator_details(pop,generation);
    }

    final_t final() const {
        return final_t(trajectories.begin(),trajectories.end());
    }

    explicit selected_mut_tracker() noexcept : trajectories(trajectories_t()) {
    }
  private:
    trajectories_t trajectories;
    template<typename pop_t>
    inline void call_operator_details(const pop_t * pop, const unsigned generation) {
        for(std::size_t i = 0 ; i < pop->mcounts.size() ; ++i ) {
            if(pop->mcounts[i]) { //if mutation is not extinct
                const auto & __m = pop->mutations[i];
                if( !__m.neutral ) {
                    const auto freq = double(pop->mcounts[i])/double(2*pop->diploids.size());
					selected_mut_data __p(__m.g,__m.pos,__m.s,__m.xtra);
                    auto __itr = trajectories.find(__p);
                    if(__itr == trajectories.end()) {
                        trajectories[__p] = std::vector<std::pair<unsigned,double> >(1,std::make_pair(generation,freq));
                    } else {
                        //Don't keep updating for fixed variants
                        if( __itr->second.back().second < 1.) {
                            __itr->second.emplace_back(generation,freq);
                        }
                    }
                }
            }
        }
    }
};
}

#endif
