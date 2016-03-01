#ifndef FWDPY_ALLELE_AGES_HPP
#define FWDPY_ALLELE_AGES_HPP

#include <limits>
#include <get_selected_mut_data.hpp>

//NAMESPACE POLLUTION
struct allele_age_data_t
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
  std::vector< allele_age_data_t > allele_ages_details( const get_selected_mut_data::final_t & trajectories,
							const double minfreq, const unsigned minsojourn );

  /*
    merge traj1 and traj2.

    traj1 is a "source" that will "sink" into a copy via std::move, hence the
    expensive-looking prototype :)
  */
  get_selected_mut_data::final_t merge_trajectories_details( get_selected_mut_data::final_t traj1,
							     const get_selected_mut_data::final_t & traj2 );
}

#endif
