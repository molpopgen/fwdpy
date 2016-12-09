from libcpp.vector cimport vector
from libcpp.utility cimport pair
from fwdpy cimport uint

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
    """Matrix representation of a sample.
    Returned by :func:`fwdpy.matrix.haplotype_matrix` pr
    :func:`fwdpy.matrix.genotype_matrix`"""
    cdef readonly object neutral
    """Data for neutral positions"""
    cdef readonly object selected
    """Data for selected positions"""
    cdef readonly list neutral_positions
    """Neutral positions"""
    cdef readonly list selected_positions
    """selected positions"""
    cdef readonly list neutral_popfreq
    """Neutral frequencies in simulation"""
    cdef readonly list selected_popfreq
    """selected frequencies in simulation"""
    cdef readonly int nrow
    """Number of rows"""
    cdef readonly int nn
    """Number of neutral mutations"""
    cdef readonly int ns
    """Number of selected nutations"""
ctypedef pair[vector[pair[size_t,uint]],vector[pair[size_t,uint]]] key_pair

#For the following function, the argument keys is treated as if const.
#The argument is not declared const due to a current limitation with Cython
cdef key_pair remove_fixed_keys(key_pair & keys,const uint n) nogil
cdef key_pair apply_min_daf(key_pair & keys,const double n,const double x) nogil
