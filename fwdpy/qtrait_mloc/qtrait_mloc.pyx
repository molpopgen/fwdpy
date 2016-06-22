# distutils: language = c++
# distutils: sources = fwdpy/qtrait_mloc/qtrait_mloc_impl.cc

from cython.operator cimport dereference as deref,preincrement as inc
from libcpp.vector cimport vector
from fwdpy.fwdpp cimport shmodel
from fwdpy.fwdpy cimport *
from fwdpy.internal.internal cimport shwrappervec,shmodel
from fwdpy.fitness cimport multilocus_fitness

cdef extern from "qtrait_evolve_mlocus.hpp" namespace "fwdpy::qtrait" nogil:
    void evolve_qtrait_mloc_cpp( GSLrng_t * rng,
			         vector[shared_ptr[multilocus_t]] * pops,
			         vector[unique_ptr[sampler_base]] & samplers,
			         const unsigned * Nvector,
			         const size_t Nvector_length,
			         const vector[double] & neutral_mutation_rates,
			         const vector[double] & selected_mutation_rates,
                                 const vector[shmodel] & shmodels,
			         const vector[double] & within_region_rec_rates,
			         const vector[double] & between_region_rec_rates,
			         const double f,
			         const double sigmaE,
			         const double optimum,
			         const double VS,
                                 const int sample,
			         const multilocus_fitness & fitness ) except +
    
include "evolve_qtraits_mloc.pyx"
