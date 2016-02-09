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
from fwdpy.gsl cimport gsl_rng


##We will expose some low-level types from fwdpp:
cdef extern from "fwdpp/forward_types.hpp" namespace "KTfwd" nogil:
    cdef cppclass mutation_base:
        double pos
        bool neutral

cdef extern from "fwdpp/sugar/popgenmut.hpp" namespace "KTfwd" nogil:
    cdef cppclass popgenmut(mutation_base):
        unsigned g
        double s
        double h

cdef extern from "fwdpp/sugar/generalmut.hpp" namespace "KTfwd" nogil:
    cdef cppclass generalmut_vec(mutation_base):
        vector[double] s
        vector[double] h
        unsigned g        

cdef extern from "fwdpp/forward_types.hpp" namespace "KTfwd" nogil:
    cdef cppclass gamete_base[T]:
        unsigned n
        vector[size_t] mutations
        vector[size_t] smutations

cdef extern from "fwdpp/sugar/sampling.hpp" namespace "KTfwd" nogil:
    ctypedef vector[pair[double,string]] sample_t
    ctypedef pair[sample_t,sample_t] sep_sample_t
    sample_t sample[POPTYPE](gsl_rng *,const POPTYPE &,const unsigned nsam , const bool removeFixed)
    sep_sample_t sample_separate[POPTYPE](gsl_rng *,const POPTYPE &,const unsigned nsam , const bool removeFixed)
    sep_sample_t sample_separate[POPTYPE](const POPTYPE &,const vector[unsigned] & individuals, const bool removeFixed) except +
    sep_sample_t sample_separate[POPTYPE](gsl_rng *,const POPTYPE &,const unsigned deme , const unsigned nsam , const bool removeFixed) except+

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
