# distutils: language = c++
# distutils: sources = fwdpy/qtrait_mloc/qtrait_mloc_impl.cc

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

cdef extern from "qtrait_evolve_mlocus.hpp" namespace "fwdpy::qtrait" nogil:
     vector[vector[pair[uint,vector[detailed_deme_sample]]]] evolve_qtrait_mloc_sample_async( GSLrng_t * rng,
				                                                              GSLrng_t * rng_sample,
				                                                              vector[shared_ptr[multilocus_t]] * pops,
				                                                              const unsigned * Nvector,
				                                                              const size_t Nvector_length,
				                                                              const vector[double] & neutral_mutation_rates,
				                                                              const vector[double] & selected_mutation_rates,
				                                                              const vector[double] & sigma_mus,
				                                                              const vector[double] & within_region_rec_rates,
				                                                              const vector[double] & between_region_rec_rates,
				                                                              const double f,
				                                                              const double sigmaE,
				                                                              const double optimum,
				                                                              const double VS,
				                                                              const int sample,
                                                                                              const unsigned nsam)

include "evolve_qtraits_mloc.pyx"
