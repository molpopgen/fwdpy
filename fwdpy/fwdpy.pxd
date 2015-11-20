from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from libcpp.map cimport map

from fwdpy.internal.internal cimport *
from fwdpy.fwdpp cimport popgenmut,gamete_base
from fwdpy.gsl cimport gsl_rng

##Create hooks to C++ types

#Wrap the classes:
cdef extern from "types.hpp" namespace "fwdpy" nogil:
    # "Standard" popgen types
    ctypedef gamete_base[popgenmut] gamete_t
    ctypedef cpplist[gamete_t] glist_t
    ctypedef cpplist[popgenmut] mlist_t

    cdef cppclass diploid_t:
        cpplist[gamete_t].iterator first
        cpplist[gamete_t].iterator second
        double g,e,w

    ctypedef vector[diploid_t] dipvector_t

    cdef cppclass singlepop_t:
        singlepop_t(unsigned)
        const unsigned N
        const unsigned generation
        mlist_t mutations
        glist_t gametes
        dipvector_t diploids
        vector[popgenmut] fixations
        vector[unsigned] fixation_times
        unsigned gen()
        unsigned popsize()
        int sane()
        void clearTrajectories()

    cdef cppclass metapop_t:
        metapop_t(vector[unsigned])
        unsigned generation
        vector[unsigned] Ns
        mlist_t mutations
        glist_t gametes
        vector[dipvector_t] diploids
        vector[popgenmut] fixations
        vector[unsigned] fixation_times
        int sane()
        int size()

    # Types based around KTfwd::generalmut_vec
    ctypedef gamete_base[generalmut_vec] gamete_gm_vec_t
    ctypedef cpplist[gamete_gm_vec_t] glist_gm_vec_t
    ctypedef cpplist[generalmut_vec] mlist_gm_vec_t

    cdef cppclass diploid_gm_vec_t:
        cpplist[gamete_gm_vec_t].iterator first
        cpplist[gamete_gm_vec_t].iterator second
        double g,e,w

    ctypedef vector[diploid_gm_vec_t] dipvector_gm_vec_t

    cdef cppclass singlepop_gm_vec_t:
        singlepop_gm_vec_t(unsigned)
        const unsigned N
        const unsigned generation
        mlist_gm_vec_t mutations
        glist_gm_vec_t gametes
        dipvector_gm_vec_t diploids
        vector[generalmut_vec] fixations
        vector[unsigned] fixation_times
        unsigned gen()
        unsigned popsize()
        int sane()

    cdef cppclass GSLrng_t:
        GSLrng_t(unsigned)
        gsl_rng * get()

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
    cpdef clearTraj(self)
    
cdef class metapop(poptype):
    cdef shared_ptr[metapop_t] mpop
    cpdef gen(self)
    cpdef popsizes(self)
    cpdef sane(self)

cdef class singlepop_gm_vec(poptype):
    cdef shared_ptr[singlepop_gm_vec_t] pop
    cpdef gen(self)
    cpdef popsize(self)
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

cdef class popvec_gmv(popcont):
    cdef vector[shared_ptr[singlepop_gm_vec_t]] pops
    cdef public object pypops
    cpdef size(self)
    cdef reset(self,const vector[shared_ptr[singlepop_gm_vec_t]] newpops)

cdef class mpopvec(popcont):
    cdef vector[shared_ptr[metapop_t]] mpops
    cdef public object pympops
    cpdef size(self)
    cdef reset(self,const vector[shared_ptr[metapop_t]]  & mpops)
    
cdef class GSLrng:
    cdef GSLrng_t * thisptr

#Typedefs for convenience
ctypedef cpplist[popgenmut].iterator mlist_t_itr
ctypedef vector[mlist_t_itr] mut_container_t
ctypedef cpplist[gamete_t].iterator glist_t_itr
ctypedef vector[diploid_t].iterator dipvector_t_itr

##Define some low-level functions that may be useful for others
cdef get_mutation( const cpplist[popgenmut].iterator & )
cdef get_gamete( const cpplist[gamete_t].iterator & )
cdef get_diploid( const vector[diploid_t].iterator & itr )

##Now, wrap the functions.
##To whatever extent possible, we avoid cdef externs in favor of Cython fxns based on cpp types.
##Many of the functions below rely on templates or other things that are too complex for Cython to handle at the moment
cdef extern from "neutral.hpp" namespace "fwdpy" nogil:
    void evolve_pop(GSLrng_t * rng, vector[shared_ptr[singlepop_t]] * pops, const vector[unsigned] nlist, const double & theta, const double & rho)

cdef extern from "sample.hpp" namespace "fwdpy" nogil:
    void get_sh( const vector[pair[double,string]] & ms_sample, const singlepop_t * pop, vector[double] * s,vector[double] * h, vector[double] * p, vector[double] * a)
    void get_sh( const vector[pair[double,string]] & samples, const metapop_t * pop, vector[double] * s, vector[double] * h, vector[double] * p, vector[double] * a)
    
cdef extern from "deps.hpp" namespace "fwdpy" nogil:
    vector[string] fwdpy_dependencies()
    vector[string] fwdpy_version()

cdef extern from "metapop.hpp" namespace "fwdpy" nogil:
    void re_init_mpop( metapop_t * mpop, const singlepop_t * pop)
    void copy_deme( metapop_t * mpop, const size_t i, const int update_counts)

cdef extern from "evolve_regions.hpp" namespace "fwdpy" nogil:
    void evolve_regions_t( GSLrng_t * rng, vector[shared_ptr[singlepop_t]] * pops,
                           const unsigned * popsizes,
                           const size_t popsizes_len,
                           const double mu_neutral,
                           const double mu_selected,
                           const double littler,
                           const double f,
                           const int track,
                           const region_manager * rm,
                           const char * fitness)

    void evolve_regions_t( GSLrng_t * rng, shared_ptr[singlepop_t] pops,
                           const unsigned * popsizes,
                           const size_t popsizes_len,
                           const double mu_neutral,
                           const double mu_selected,
                           const double littler,
                           const double f,
                           const int track,
                           const region_manager * rm,
                           const char * fitness)

    void split_and_evolve_t(GSLrng_t * rng,
                vector[shared_ptr[metapop_t]] * mpops,
                const unsigned * Nvector_A,
                const size_t Nvector_A_len,
                const unsigned * Nvector_B,
                const size_t Nvector_B_len,
                const double neutral,
                const double selected,
                const double recrate,
                const vector[double] & fs,
                const region_manager * rm,
                const char * fitness)

    vector[shared_ptr[singlepop_t]]  evolve_regions_async(const unsigned npops,
							  GSLrng_t * rng, 
							  const unsigned * Nvector,
							  const size_t Nvector_len,
							  const double mu_neutral,
							  const double mu_selected,
							  const double littler,
							  const double f,
							  const int track,
							  const region_manager * rm,
							  const char * fitness)

cdef extern from "trajectories.hpp" namespace "fwdpy" nogil:
    map[string,vector[double] ] get_singlepop_traj(const singlepop_t *pop,const unsigned minsojourn,const double minfreq)
