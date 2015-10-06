from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from libcpp.map cimport map

from fwdpy.internal.internal cimport *

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
    cdef cppclass metapop_t:
        metapop_t(vector[unsigned])
        unsigned generation
        vector[unsigned] Ns
        int sane()
        int size()
    cdef cppclass GSLrng_t:
        GSLrng_t(unsigned)

#Now, provied definitions for classes in classes.pyx
cdef class poptype(object):
    """
    Empty base class for a population object.

    Example derived types include :class:`fwdpy.fwdpy.singlepop` and :class:`fwdpy.fwdpy.metapop`
    """
    pass

cdef class singlepop(poptype):
    cdef shared_ptr[singlepop_t] pop
    cpdef gen(self)
    cpdef popsize(self)
    cpdef sane(self)
    
cdef class metapop(poptype):
    cdef shared_ptr[metapop_t] mpop
    cpdef gen(self)
    cpdef popsizes(self)
    cpdef sane(self)
    
cdef class popcont(object):
    """
    Empty base class for containers of population objects.

    Example derived types include :class:`fwdpy.fwdpy.popvec` and :class:`fwdpy.fwdpy.mpopvec`
    """
    pass

cdef class popvec(popcont):
    cdef vector[shared_ptr[singlepop_t]] pops
    cdef public object pypops
    cpdef size(self)
    cdef reset(self,const vector[shared_ptr[singlepop_t]] newpops)
    
cdef class mpopvec(popcont):
    cdef vector[shared_ptr[metapop_t]] mpops
    cdef public object pympops
    cpdef size(self)
    cdef reset(self,const vector[shared_ptr[metapop_t]]  & mpops)
    
cdef class GSLrng:
    cdef GSLrng_t * thisptr


    
##Now, wrap the functions
cdef extern from "neutral.hpp" namespace "fwdpy":
    void evolve_pop(GSLrng_t * rng, vector[shared_ptr[singlepop_t]] * pops, const vector[unsigned] nlist, const double & theta, const double & rho)

cdef extern from "sample.hpp" namespace "fwdpy":
    vector[pair[double,string]] take_sample_from_pop(GSLrng_t * rng,const singlepop_t * pop,const unsigned nsam, const int remove_fixed)
    pair[vector[pair[double,string]],vector[pair[double,string]]] take_sample_from_pop_sep(GSLrng_t * rng,const singlepop_t * pop,const unsigned nsam, const int remove_fixed)
    pair[vector[pair[double,string]],vector[pair[double,string]]] take_sample_from_metapop_sep(GSLrng_t * rng,const metapop_t * mpop,const unsigned & nsam, const int remove_fixed, const int deme)
    pair[vector[pair[double,string]],vector[pair[double,string]]] sample_specific_diploids(const singlepop_t * pop, const vector[unsigned] & indlist, const int remove_fixed)
    void get_sh( const vector[pair[double,string]] & ms_sample, const singlepop_t * pop, vector[double] * s,vector[double] * h, vector[double] * p, vector[double] * a)
    void get_sh( const vector[pair[double,string]] & samples, const metapop_t * pop, vector[double] * s, vector[double] * h, vector[double] * p, vector[double] * a)
    map[string,vector[double]] diploid_view_cpp(const singlepop_t *pop, const size_t ind, const int remove_fixed) except +
    map[string,vector[double]] diploid_view_cpp(const metapop_t * pop, const size_t ind, const int remove_fixed, const int deme) except +
    
cdef extern from "deps.hpp" namespace "fwdpy":
    vector[string] fwdpy_dependencies()
    vector[string] fwdpy_version()

cdef extern from "metapop.hpp" namespace "fwdpy":
    void re_init_mpop( metapop_t * mpop, const singlepop_t * pop)
    void copy_deme( metapop_t * mpop, const size_t i, const int update_counts)

cdef extern from "evolve_regions.hpp" namespace "fwdpy":
    void evolve_regions_t( GSLrng_t * rng, vector[shared_ptr[singlepop_t]] * pops,
		       const unsigned * popsizes,
               const size_t popsizes_len,
		       const double mu_neutral,
		       const double mu_selected,
		       const double littler,
		       const double f,
		       const region_manager * rm,
		       const char * fitness)

    void split_and_evolve_t(GSLrng_t * rng,
                vector[shared_ptr[metapop_t]] * mpops,
                const unsigned * Nvector_A,
                const size_t Nvector_A_len,
                const unsigned * Nvector_B,
                const size_t Nvector_B_len,
                const double & neutral,
                const double & selected,
                const double & recrate,
                const vector[double] & fs,
                const region_manager * rm,
                const char * fitness)
