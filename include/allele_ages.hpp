#ifndef FWDPY_ALLELE_AGES_HPP
#define FWDPY_ALLELE_AGES_HPP

#include <limits>
#include "sampler_selected_mut_tracker.hpp"

//NAMESPACE POLLUTION
struct allele_age_data_t
/*!
  \brief Structure representing allele ages.

  Intended use is conversion to Python dictionary.

  Returned to Python via Cython making a call to fwdpy::allele_ages_details
*/
{
  double esize,max_freq,last_freq;
  unsigned origin,tlen;
  allele_age_data_t( double e,
		     double mf,
		     double lf,
		     unsigned o,
		     unsigned t ) : esize(e),
				    max_freq(mf),
				    last_freq(lf),
				    origin(o),
				    tlen(t)
  {
  }
  allele_age_data_t() : esize(std::numeric_limits<double>::quiet_NaN()),
			max_freq(std::numeric_limits<double>::quiet_NaN()),
			last_freq(std::numeric_limits<double>::quiet_NaN()),
			origin(std::numeric_limits<unsigned>::max()),
			tlen(std::numeric_limits<unsigned>::max())
  {
  }
};

namespace fwdpy
{

  /*!
    \brief Take an object of type fwdpy::selected_mut_tracker and return allele age data.

    The parameters exclude mutations whose max frequency are < minfreq or whose sojourn time
    in simulation is < minsojourn.
  */
  std::vector< allele_age_data_t > allele_ages_details( const selected_mut_tracker::final_t & trajectories,
							const double minfreq, const unsigned minsojourn );

  /*
    \brief Merge containers of mutation trajectories.

    \note traj1 is a "source" that will "sink" into a copy via std::move, hence the
    expensive-looking prototype :)
  */
  selected_mut_tracker::final_t merge_trajectories_details( selected_mut_tracker::final_t traj1,
							    const selected_mut_tracker::final_t & traj2 );
}

#endif
