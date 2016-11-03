from libcpp.memory cimport unique_ptr
from cython_gsl cimport gsl_vector,gsl_matrix

cdef extern from "gsl.hpp" namespace "fwdpy::gsl" nogil:
    cdef cppclass gsl_matrix_deleter:
        pass
    cdef cppclass gsl_vector_deleter:
        pass

    ctypedef unique_ptr[gsl_vector,gsl_vector_deleter] gsl_vector_ptr_t
    ctypedef unique_ptr[gsl_matrix,gsl_matrix_deleter] gsl_matrix_ptr_t
