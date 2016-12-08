from libcpp.vector cimport vector
from libcpp.utility cimport pair
from fwdpy cimport uint

cdef extern from "fwdpp/sugar/matrix.hpp" namespace "KTfwd" nogil:
    cdef struct data_matrix:
        vector[short] neutral,selected
        vector[double] neutral_positions,selected_positions,neutral_popfreq,selected_popfreq
        size_t nrow

    pair[vector[pair[size_t,uint]],vector[pair[size_t,uint]]] mutation_keys[POPTYPE](const POPTYPE & pop, 
            const vector[size_t] & inividuals,
            const bint include_neutral, const bint include_selected, const size_t deme) except +

    data_matrix genotype_matrix[POPTYPE](const POPTYPE & pop,
            const vector[size_t] & inividuals,
            const vector[pair[size_t,uint]] & neutral_keys,
            const vector[pair[size_t,uint]] & selected_keys, const size_t deme) except +

    data_matrix haplotype_matrix[POPTYPE](const POPTYPE & pop,
            const vector[size_t] & inividuals,
            const vector[pair[size_t,uint]] & neutral_keys,
            const vector[pair[size_t,uint]] & selected_keys, const size_t deme) except +

ctypedef pair[vector[pair[size_t,uint]],vector[pair[size_t,uint]]] key_pair

#For the following function, the argument keys is treated as if const.
#The argument is not declared const due to a current limitation with Cython
cdef key_pair remove_fixed_keys(key_pair & keys,const uint n) nogil
cdef key_pair apply_min_daf(key_pair & keys,const double n,const double x) nogil
