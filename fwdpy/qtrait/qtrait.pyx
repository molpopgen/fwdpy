# distutils: language = c++
# distutils: sources = fwdpy/qtrait/qtrait_impl.cc fwdpy/qtrait/ew2010.cc fwdpy/qtrait/ewvw.cc
from cython.operator cimport dereference as deref,preincrement as inc
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.utility cimport pair
from fwdpy.fwdpy cimport *
from fwdpy.internal.internal cimport shwrappervec
import fwdpy.internal as internal
import pandas

cdef extern from "qtraits.hpp" namespace "fwdpy::qtrait":
    void evolve_qtraits_t( GSLrng_t * rng, vector[shared_ptr[singlepop_t] ] * pops,
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
        const vector[double] & nbegs,
        const vector[double] & nends,
        const vector[double] & nweights,
        const vector[double] & sbegs,
        const vector[double] & sends,
        const vector[double] & sweights,
        const vector[shmodel] * callbacks,
        const vector[double] & rbeg,
        const vector[double] & rend,
        const vector[double] & rweight)
    map[string,double] qtrait_pop_props( const singlepop_t * pop );
    map[string,vector[double]] get_qtrait_traj(const singlepop_t *pop,const unsigned minsojourn,const double minfreq)
    map[string,vector[double]] qtrait_esize_freq(const singlepop_t * pop)
    cdef struct ew_mut_details:
       double s
       double e
       double p
    map[double,ew_mut_details] ew2010_assign_effects(GSLrng_t * rng, const singlepop_t * pop, const double tau, const double sigma) except +
    vector[double] ew2010_traits_cpp(const singlepop_t * pop, const map[double,ew_mut_details] & effects) except +

cdef extern from "ewvw.hpp" namespace "fwdpy::qtrait":
    void evolve_ewvw_t( GSLrng_t * rng,
        vector[shared_ptr[singlepop_t]] * pops,
        const unsigned * Nvector,
        const size_t Nvector_length,
        const double mu_neutral,
        const double mu_selected,
        const double littler,
        const double f,
        const double sigmaE,
        const double VS_total,
        const double optimum,
        const int track,
        const vector[double] & nbegs,
        const vector[double] & nends,
        const vector[double] & nweights,
        const vector[double] & sbegs,
        const vector[double] & sends,
        const vector[double] & sweights,
        const vector[shmodel] * callbacks,
        const vector[double] & rbeg,
        const vector[double] & rend,
        const vector[double] & rweight)
    
include "evolve_qtraits.pyx"
include "ew2010.pyx"
include "misc.pyx"



