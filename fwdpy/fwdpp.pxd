"""
This file exposes the fwdpp libray to Cython.

It is currently not possible (and not necessarily even desirable) to expose all of fwdpp.

Cython currently doesn't seem to support all of the possible type of C++ template declaration necessary,
making this file a work in progress.

It is possible that some of these issues can be addresses via further evolution of the fwdpp 'sugar' layer,
simplifying the API considerably over the "raw" fwdpp API.
"""

from libcpp cimport bool
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.list cimport list as cpplist

cdef extern from "gsl/gsl_rng.h" nogil:
    ctypedef struct gsl_rng

##We will expose some low-level types from fwdpp:
cdef extern from "fwdpp/forward_types.hpp" namespace "KTfwd" nogil:
    cdef cppclass mutation_base:
        double pos
        unsigned n
        bool neutral

cdef extern from "fwdpp/sugar/popgenmut.hpp" namespace "KTfwd" nogil:
    cdef cppclass popgenmut(mutation_base):
        unsigned g
        double s
        double h

cdef extern from "fwdpp/forward_types.hpp" namespace "KTfwd" nogil:
    cdef cppclass gamete_base[popgenmut]:
        unsigned n
        vector[cpplist[popgenmut].iterator] mutations
        vector[cpplist[popgenmut].iterator] smutations

cdef extern from "fwdpp/sugar/sampling.hpp" namespace "KTfwd" nogil:
    ctypedef vector[pair[double,string]] sample_t
    ctypedef pair[sample_t,sample_t] sep_sample_t
    sample_t sample[POPTYPE](gsl_rng *,const POPTYPE &,const unsigned nsam , const bool removeFixed)
    sep_sample_t sample_separate[POPTYPE](gsl_rng *,const POPTYPE &,const unsigned nsam , const bool removeFixed)
    sep_sample_t sample_separate[POPTYPE](gsl_rng *,const POPTYPE &,const unsigned deme , const unsigned nsam , const bool removeFixed)

## fwdpp's extensions sub-library:    
cdef extern from "fwdpp/extensions/callbacks.hpp" namespace "KTfwd::extensions":
    cdef cppclass shmodel:
        shmodel()
    cdef cppclass constant:
        constant(double)
    cdef cppclass exponential:
        exponential(double)
    cdef cppclass uniform:
        uniform(double,double)
    cdef cppclass beta:
        beta(double,double,double)
    cdef cppclass gaussian:
        gaussian(double)
    cdef cppclass gamma:
        gamma(double,double)

cdef extern from "internal/callbacks.hpp" namespace "fwdpy::internal":
    void make_gamma_s(shmodel *, double,double)
    void make_constant_s(shmodel * s, const double scoeff);
    void make_uniform_s(shmodel * s, const double lo, const double hi);
    void make_exp_s(shmodel * s, const double mean);
    void make_gaussian_s(shmodel * s, const double sd);
    void make_constant_h(shmodel * s, const double h);

cdef extern from "internal/internal.hpp" namespace "fwdpy::internal":
    cdef cppclass region_manager:
        region_manager()
        vector[shmodel] callbacks
        vector[double] nb
        vector[double] ne
        vector[double] nw
        vector[double] sb
        vector[double] se
        vector[double] sw
        vector[double] rb
        vector[double] re
        vector[double] rw
