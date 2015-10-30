##This file is for "gory internal details" things from fwdpp
from fwdpy.gsl cimport gsl_ran_discrete_t
from fwdpy.cpp cimport unique_ptr
    
cdef extern from "fwdpp/internal/gsl_discrete.hpp" namespace "KTfwd::internal" nogil:
    cdef cppclass gsl_ran_discrete_t_deleter:
        pass
    ctypedef unique_ptr[gsl_ran_discrete_t,gsl_ran_discrete_t_deleter] gsl_ran_discrete_t_ptr 
