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
    vector[pair[double,string]] take_sample_from_pop(GSLrng_t * rng,const singlepop_t * pop,const unsigned nsam, const int remove_fixed)
    pair[vector[pair[double,string]],vector[pair[double,string]]] sample_specific_diploids(const singlepop_t * pop, const vector[unsigned] & indlist, const int remove_fixed)
    double tajd( const vector[pair[double,string]] & __data )
    void get_sh( const vector[pair[double,string]] & ms_sample, const singlepop_t * pop, vector[double] * s,vector[double] * h, vector[double] * p, vector[double] * a)

cdef extern from "deps.hpp" namespace "fwdpy":
    vector[string] fwdpy_dependencies()
    vector[string] fwdpy_version()

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

cdef extern from "evolve_regions.hpp" namespace "fwdpy":
    void evolve_regions_t( GSLrng_t * rng, vector[shared_ptr[singlepop_t]] * pops,
		       const unsigned * popsizes,
               const size_t popsizes_len,
		       const double mu_neutral,
		       const double mu_selected,
		       const double littler,
		       const double f,
		       const vector[double] & nbegs,
		       const vector[double] & nends,
		       const vector[double] & nweights,
		       const vector[double] & sbegs,
		       const vector[double] & sends,
		       const vector[double] & sweights,
		       const vector[shmodel] * callbacks,
		       const vector[double] & rbeg,
		       const vector[double] & rend,
		       const vector[double] & rweight,
		       const char * fitness)

cdef extern from "callbacks.hpp" namespace "fwdpy":
    void make_gamma_s(shmodel *, double,double)
    void make_constant_s(shmodel * s, const double scoeff);
    void make_uniform_s(shmodel * s, const double lo, const double hi);
    void make_exp_s(shmodel * s, const double mean);
    void make_gaussian_s(shmodel * s, const double sd);
    void make_constant_h(shmodel * s, const double h);
