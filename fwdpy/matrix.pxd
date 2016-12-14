from libcpp.vector cimport vector
from libcpp.utility cimport pair
from fwdpy cimport uint,mcont_t,ucont_t

cdef extern from "fwdpp/sugar/matrix.hpp" namespace "KTfwd" nogil:
    cdef struct data_matrix:
        vector[char] neutral,selected
        vector[double] neutral_positions,selected_positions,neutral_popfreq,selected_popfreq
        size_t nrow

    pair[vector[pair[size_t,uint]],vector[pair[size_t,uint]]] mutation_keys[POPTYPE](const POPTYPE & pop,
            const vector[size_t] & inividuals,
            const bint include_neutral, 
            const bint include_selected, const size_t deme) except +

    data_matrix fwdpp_genotype_matrix"KTfwd::genotype_matrix"[POPTYPE](const POPTYPE & pop,
            const vector[size_t] & inividuals,
            const vector[pair[size_t,uint]] & neutral_keys,
            const vector[pair[size_t,uint]] & selected_keys, const size_t deme) except +

    data_matrix fwdpp_haplotype_matrix"KTfwd::haplotype_matrix"[POPTYPE](const POPTYPE & pop,
            const vector[size_t] & inividuals,
            const vector[pair[size_t,uint]] & neutral_keys,
            const vector[pair[size_t,uint]] & selected_keys, const size_t deme) except +

cdef class DataMatrix(object):
    cdef readonly object neutral
    """Data for neutral positions"""
    cdef readonly object selected
    """Data for selected positions"""
    cdef readonly object neutral_positions
    """Neutral positions"""
    cdef readonly object selected_positions
    """selected positions"""
    cdef readonly object neutral_popfreq
    """Neutral frequencies in simulation"""
    cdef readonly object selected_popfreq
    """selected frequencies in simulation"""
    cdef readonly int nrow
    """Number of rows"""
    cdef readonly int nn
    """Number of neutral mutations"""
    cdef readonly int ns
    """Number of selected nutations"""

cdef class GenotypeMatrix(DataMatrix):
    pass

cdef class HaplotypeMatrix(DataMatrix):
    pass

ctypedef pair[vector[pair[size_t,uint]],vector[pair[size_t,uint]]] key_pair

#Uses some python list tricks, and hence nogil
cdef sort_keys_by_position(key_pair & keys,const mcont_t & mutations)

#For the following functions, the argument keys is treated as if const.
#The argument is not declared const due to a current limitation with Cython
cdef key_pair remove_fixed_keys(key_pair & keys,const uint n) nogil
cdef key_pair apply_min_sample_daf(key_pair & keys,const double n,const double x) nogil
cdef key_pair apply_min_pop_daf(key_pair & keys,const unsigned twoN,const double q, const ucont_t & mcounts) nogil
