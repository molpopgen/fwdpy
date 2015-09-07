from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr

##Create hooks to C++ types

#Wrap the classes:
cdef extern from "types.hpp" namespace "fwdpy":
    cdef cppclass popvector:
        popvector(unsigned,unsigned)
        unsigned size()
        unsigned generation(unsigned)
        unsigned popsize(unsigned)
        int sane(unsigned)
    cdef cppclass singlepop_t:
        singlepop_t(unsigned)
        unsigned gen()
        unsigned popsize()
        int sane()
    cdef cppclass GSLrng_t:
        GSLrng_t(unsigned)

##Now, wrap the functions
cdef extern from "neutral.hpp" namespace "fwdpy":
    void evolve_pop(GSLrng_t * rng, vector[shared_ptr[singlepop_t]] * pops, const vector[unsigned] nlist, const double & theta, const double & rho)

cdef extern from "sample.hpp" namespace "fwdpy":
    vector[vector[pair[double,string]]] take_sample_from_pop(GSLrng_t * rng,const vector[shared_ptr[singlepop_t]] & pops,const unsigned & nsam)
    double tajd( const vector[pair[double,string]] & __data )
    void get_sh( const vector[vector[pair[double,string]]] & samples, const vector[shared_ptr[singlepop_t]] & pops, const unsigned i,	vector[double] * s,vector[double] * h, vector[double] * p, vector[double] * a)
