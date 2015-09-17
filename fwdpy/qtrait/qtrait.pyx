# distutils: language = c++
# distutils: sources = fwdpy/qtrait/qtrait_impl.cc
from libcpp.vector cimport vector
from fwdpy.fwdpy cimport *

cdef extern from "evolve_qtraits.hpp" namespace "fwdpy::qtrait":
    void evolve_qtraits_t( GSLrng_t * rng, vector[shared_ptr[singlepop_t] ] * pops,
        const unsigned * Nvector,
        const size_t Nvector_length,
        const double mu_neutral,
        const double mu_selected,
        const double littler,
        const double f,
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
        const vector[double] & rweight,
        const char * fitness)
