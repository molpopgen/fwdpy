from libcpp.vector cimport vector

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

cdef extern from "callbacks.hpp" namespace "fwdpy":
    void make_gamma_s(shmodel *, double,double)
    void make_constant_s(shmodel * s, const double scoeff);
    void make_uniform_s(shmodel * s, const double lo, const double hi);
    void make_exp_s(shmodel * s, const double mean);
    void make_gaussian_s(shmodel * s, const double sd);
    void make_constant_h(shmodel * s, const double h);

##These are the callback wrappers from fwdpp
cdef class shwrappervec:
    """
    Wrapper for a vector of callback objects from fwdpp's extension library.

    Users will not interact with this type directly.  Rather, it is used
    by other module functions to process user inputs.
    """
    cdef vector[shmodel] vec
