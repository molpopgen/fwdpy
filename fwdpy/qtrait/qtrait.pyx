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
from fwdpy.structs cimport VAcum

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
    void evolve_regions_qtrait_cpp(GSLrng_t * rng,
				   vector[shared_ptr[singlepop_t]] * pops,
				   vector[unique_ptr[sampler_base]] & samplers,
				   const unsigned * Nvector,
				   const size_t Nvector_length,
				   const double neutral,
				   const double selected,
				   const double recrate,
				   const double f,
				   const double sigmaE,
				   const double optimum,
				   const double VS,
				   const int interval,
				   const region_manager * rm,
				   const singlepop_fitness & fitness) except +

include "evolve_qtraits.pyx"
include "ew2010.pyx"
include "misc.pyx"



