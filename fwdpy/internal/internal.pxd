from libcpp.vector cimport vector
from libcpp.list cimport list as cpplist 

from fwdpy.fwdpp cimport *
from fwdpy.fwdpy cimport singlepop_t,metapop_t

##These are the callback wrappers from fwdpp
cdef class shwrappervec:
    """
    Wrapper for a vector of callback objects from fwdpp's extension library.

    Users will not interact with this type directly.  Rather, it is used
    by other module functions to process user inputs.
    """
    cdef vector[shmodel] vec

##Quick wrapper for region_manager
cdef class region_manager_wrapper:
    cdef region_manager * thisptr



