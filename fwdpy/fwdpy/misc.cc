#include <types.hpp>

namespace fwdpy
{
  std::vector<qtrait_stats_cython> convert_qtrait_stats( const qtrait_stats_t & qts )
  {
    std::vector<qtrait_stats_cython> rv;
    for( const auto & i : qts )
      {
	rv.emplace_back( "VG",i[std::size_t(qtrait_stat_names::VG)],i[std::size_t(qtrait_stat_names::GEN)] );
	rv.emplace_back( "VE",i[std::size_t(qtrait_stat_names::VE)],i[std::size_t(qtrait_stat_names::GEN)] );
	rv.emplace_back( "leading_q",i[std::size_t(qtrait_stat_names::PLF)],i[std::size_t(qtrait_stat_names::GEN)] );
	rv.emplace_back( "leading_e",i[std::size_t(qtrait_stat_names::LE)],i[std::size_t(qtrait_stat_names::GEN)] );
	rv.emplace_back( "max_expl",i[std::size_t(qtrait_stat_names::MAXEXP)],i[std::size_t(qtrait_stat_names::GEN)] );
	rv.emplace_back( "ebar",i[std::size_t(qtrait_stat_names::EBAR)],i[std::size_t(qtrait_stat_names::GEN)] );
	rv.emplace_back( "wbar",i[std::size_t(qtrait_stat_names::WBAR)],i[std::size_t(qtrait_stat_names::GEN)] );
      }
    return rv;
  }
}
