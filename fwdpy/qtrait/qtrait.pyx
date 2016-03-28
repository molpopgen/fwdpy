# distutils: language = c++
# distutils: sources = fwdpy/qtrait/qtrait_impl.cc fwdpy/qtrait/ew2010.cc

from cython.operator cimport dereference as deref,preincrement as inc
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.utility cimport pair
from libcpp.string cimport string
from fwdpy.fwdpy cimport *
from fwdpy.internal.internal cimport shwrappervec
from fwdpy.fwdpp cimport sep_sample_t
import fwdpy.internal as internal
import pandas

cdef extern from "qtrait_evolve_rules.hpp" namespace "fwdpy::qtrait" nogil:
    cdef cppclass qtrait_model_rules:
        qtrait_model_rules(const double & sigmaE, const double & optimum, const double & VS, const unsigned maxN)

cdef extern from "qtrait_pleiotropic.hpp" namespace "fwdpy::qtrait" nogil:    
    cdef struct ew_mut_details:
        double s
        double e
        double p
        
    map[double,ew_mut_details] ew2010_assign_effects(GSLrng_t * rng, const singlepop_t * pop, const double tau, const double sigma) except +
    vector[double] ew2010_traits_cpp(const singlepop_t * pop, const map[double,ew_mut_details] & effects) except +

cdef extern from "qtrait_evolve.hpp" namespace "fwdpy::qtrait" nogil:
    void evolve_qtrait_no_sampling_async( GSLrng_t * rng,
                                          vector[shared_ptr[singlepop_t]] * pops,
				          const unsigned * Nvector,
				          const size_t Nvector_length,
				          const double mu_neutral,
				          const double mu_selected,
				          const double littler,
				          const double f,
				          const double sigmaE,
				          const double optimum,
				          const double VS,
				          const region_manager * rm)

    
    vector[vector[pair[uint,detailed_deme_sample]]] evolve_qtrait_sample_async( GSLrng_t * rng,
                                                                                vector[shared_ptr[singlepop_t]] * pops,
				                                                const unsigned * Nvector,
				                                                const size_t Nvector_length,
				                                                const double mu_neutral,
				                                                const double mu_selected,
				                                                const double littler,
				                                                const double f,
				                                                const double sigmaE,
				                                                const double optimum,
				                                                const double VS,
				                                                const int sample,
				                                                const unsigned nsam,
				                                                const region_manager * rm)

    vector[vector[qtrait_stats_cython]] evolve_qtrait_popstats_async( GSLrng_t * rng,
                                                                      vector[shared_ptr[singlepop_t]] * pops,
				                                      const unsigned * Nvector,
				                                      const size_t Nvector_length,
				                                      const double mu_neutral,
				                                      const double mu_selected,
				                                      const double littler,
				                                      const double f,
				                                      const double sigmaE,
				                                      const double optimum,
				                                      const double VS,
				                                      const int sample,
				                                      const region_manager * rm)

    vector[vector[pair[selected_mut_data,vector[double]]]] evolve_qtrait_track_async( GSLrng_t * rng,
                                                                                              vector[shared_ptr[singlepop_t]] * pops,
				                                                              const unsigned * Nvector,
				                                                              const size_t Nvector_length,
				                                                              const double mu_neutral,
				                                                              const double mu_selected,
				                                                              const double littler,
				                                                              const double f,
				                                                              const double sigmaE,
				                                                              const double optimum,
				                                                              const double VS,
				                                                              const int track,
				                                                              const region_manager * rm)

    void evolve_gbr_no_sampling_async( GSLrng_t * rng,
                                       vector[shared_ptr[singlepop_t]] * pops,
				       const unsigned * Nvector,
				       const size_t Nvector_length,
				       const double mu_neutral,
				       const double mu_selected,
				       const double littler,
				       const double f,
				       const double sigmaE,
				       const double optimum,
				       const double VS,
				       const region_manager * rm)
    
    
    vector[vector[pair[uint,detailed_deme_sample]]] evolve_gbr_sample_async( GSLrng_t * rng,
                                                                             vector[shared_ptr[singlepop_t]] * pops,
				                                             const unsigned * Nvector,
				                                             const size_t Nvector_length,
				                                             const double mu_neutral,
				                                             const double mu_selected,
				                                             const double littler,
				                                             const double f,
				                                             const double sigmaE,
				                                             const double optimum,
				                                             const double VS,
				                                             const int sample,
				                                             const unsigned nsam,
				                                             const region_manager * rm)

    vector[vector[qtrait_stats_cython]] evolve_gbr_popstats_async( GSLrng_t * rng,
                                                                   vector[shared_ptr[singlepop_t]] * pops,
				                                   const unsigned * Nvector,
				                                   const size_t Nvector_length,
				                                   const double mu_neutral,
				                                   const double mu_selected,
				                                   const double littler,
				                                   const double f,
				                                   const double sigmaE,
				                                   const double optimum,
				                                   const double VS,
				                                   const int sample,
				                                   const region_manager * rm)

    vector[vector[pair[selected_mut_data,vector[double]]]] evolve_gbr_track_async( GSLrng_t * rng,
                                                                                           vector[shared_ptr[singlepop_t]] * pops,
				                                                           const unsigned * Nvector,
				                                                           const size_t Nvector_length,
				                                                           const double mu_neutral,
				                                                           const double mu_selected,
				                                                           const double littler,
				                                                           const double f,
				                                                           const double sigmaE,
				                                                           const double optimum,
                                                                                           const double VS,
				                                                           const int track,
				                                                           const region_manager * rm)
    
include "evolve_qtraits.pyx"
include "ew2010.pyx"
include "misc.pyx"



