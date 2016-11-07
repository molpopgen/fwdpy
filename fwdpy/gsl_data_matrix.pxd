from libcpp.vector cimport vector
from libcpp.string cimport string as cppstring
from libcpp.utility cimport pair
from libcpp.memory cimport unique_ptr
from fwdpy.fwdpy cimport uint
from fwdpy.gsl cimport gsl_matrix_ptr_t 
from cython_gsl cimport gsl_matrix 

cdef extern from "gsl_data_matrix.hpp" namespace "fwdpy::gsl_data_matrix" nogil:
    vector[uint] get_mut_keys[POPTYPE](const POPTYPE * pop, const bint sort_freq, const bint sort_esizes) 
    void update_matrix_counts[POPTYPE] (const POPTYPE *pop, const vector[uint] &mut_keys, gsl_matrix * m)
    cdef cppclass geno_matrix:
        geno_matrix(size_t nrow,size_t ncol)
        vector[double] G
        gsl_matrix_ptr_t m
        size_t nrow,ncol
    void write_geno_matrix(const geno_matrix *m, const uint generation,
                      cppstring stub, const int repid,
                      const bint keep_origin)
    void emplace_move(vector[pair[uint,unique_ptr[geno_matrix]]] & v,pair[uint,unique_ptr[geno_matrix]] & p)


