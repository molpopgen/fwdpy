from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr

##Create hooks to C++ types

#Wrap the classes:
cdef extern from "types.hpp" namespace "fwdpy":
    cdef cppclass singlepop_t:
        singlepop_t(unsigned)
        const unsigned N
        const unsigned generation
        unsigned gen()
        unsigned popsize()
        int sane()
    cdef cppclass GSLrng_t:
        GSLrng_t(unsigned)

##Now, wrap the functions
cdef extern from "neutral.hpp" namespace "fwdpy":
    void evolve_pop(GSLrng_t * rng, vector[shared_ptr[singlepop_t]] * pops, const vector[unsigned] nlist, const double & theta, const double & rho)

cdef extern from "sample.hpp" namespace "fwdpy":
    vector[pair[double,string]] take_sample_from_pop(GSLrng_t * rng,const singlepop_t * pop,const unsigned & nsam)
    double tajd( const vector[pair[double,string]] & __data )
    void get_sh( const vector[pair[double,string]] & ms_sample, const singlepop_t * pop, vector[double] * s,vector[double] * h, vector[double] * p, vector[double] * a)

cdef extern from "deps.hpp" namespace "fwdpy":
    vector[string] dependencies()

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
