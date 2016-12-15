from libcpp.memory cimport unique_ptr,shared_ptr
from cython_gsl cimport gsl_vector,gsl_matrix

cdef extern from "gsl.hpp" namespace "fwdpy::gsl" nogil:
    cdef cppclass gsl_matrix_deleter:
        pass
    cdef cppclass gsl_vector_deleter:
        pass

    ctypedef unique_ptr[gsl_vector] gsl_vector_shared_ptr_t
    ctypedef unique_ptr[gsl_matrix] gsl_matrix_shared_ptr_t
    ctypedef unique_ptr[gsl_vector,gsl_vector_deleter] gsl_vector_ptr_t
    ctypedef unique_ptr[gsl_matrix,gsl_matrix_deleter] gsl_matrix_ptr_t

    #Cython does not understand the constructors taking deleters,
    #So we define functions to make and return them
    gsl_matrix_shared_ptr_t make_gsl_matrix_shared_ptr_t(const size_t size1, const size_t size2)
    gsl_vector_shared_ptr_t make_gsl_vector_shared_ptr_t(const size_t size1)
