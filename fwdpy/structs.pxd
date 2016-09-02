#This file contains wrappers for structs used by
#various modules.  They are kept here b/c they
#may not be used by all modules and having them in
#fwdpy.pxd had a "whole world recompiles" issue
#when any header was touched.

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.utility cimport pair
from libc.stdint cimport uint16_t
from fwdpy.fwdpp cimport sep_sample_t

cdef extern from "sample.hpp" namespace "fwdpy" nogil:
    cdef struct popsample_details:
        vector[double] s, h, p
        vector[unsigned] age, dcount, generation, locus
        vector[uint16_t] label
    
cdef extern from "sampler_selected_mut_tracker.hpp" namespace "fwdpy" nogil:
    cdef struct selected_mut_data:
        double pos
        double esize
        unsigned origin

    cdef struct selected_mut_data_tidy:
        double pos
        double esize
        double freq
        unsigned origin
        unsigned generation
        
cdef extern from "sampler_pop_properties.hpp" namespace "fwdpy" nogil:
    cdef struct qtrait_stats_cython:
        string stat
        double value
        unsigned generation

cdef extern from "sampler_additive_variance.hpp" namespace "fwdpy" nogil:
    cdef struct VAcum:
        double freq
        double pss
        unsigned generation
        unsigned N

cdef extern from "allele_ages.hpp" namespace "fwdpy" nogil:
    cdef struct allele_age_data_t:
        double esize
        double max_freq
        double last_freq
        unsigned origin
        unsigned tlen

cdef extern from "haplotype_matrix.hpp" namespace "fwdpy" nogil:
    cdef struct haplotype_matrix:
        vector[size_t] n
        vector[size_t] s
        vector[double] np
        vector[double] sp
        vector[double] nf
        vector[double] sf
        vector[double] G
        vector[double] E
        vector[double] w
        vector[double] esizes
        vector[double] h
        size_t nrow
        size_t ncol_n
        size_t ncol_s
