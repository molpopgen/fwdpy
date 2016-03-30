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
    ctypedef gamete_base[void] gamete_t
    ctypedef vector[gamete_t] gcont_t
    ctypedef vector[popgenmut] mcont_t

    cdef cppclass diploid_t:
        size_t first,second;
        double g,e,w

    ctypedef vector[diploid_t] dipvector_t

    cdef cppclass singlepop_t:
        singlepop_t(unsigned)
        unsigned N
        unsigned generation
        mcont_t mutations
        vector[unsigned] mcounts
        gcont_t gametes
        dipvector_t diploids
        vector[popgenmut] fixations
        vector[unsigned] fixation_times
        unsigned gen()
        unsigned popsize()
        int sane()

    cdef cppclass metapop_t:
        metapop_t(vector[unsigned])
        unsigned generation
        vector[unsigned] Ns
        mcont_t mutations
        vector[unsigned] mcounts
        gcont_t gametes
        vector[dipvector_t] diploids
        vector[popgenmut] fixations
        vector[unsigned] fixation_times
        int sane()
        int size()

    # Types based around KTfwd::generalmut_vec
    ctypedef vector[generalmut_vec] mlist_gm_vec_t

    cdef cppclass singlepop_gm_vec_t:
        singlepop_gm_vec_t(unsigned)
        const unsigned N
        const unsigned generation
        mlist_gm_vec_t mutations
        vector[unsigned] mcounts
        gcont_t gametes
        vector[diploid_t] diploids
        vector[generalmut_vec] fixations
        vector[unsigned] fixation_times
        unsigned gen()
        unsigned popsize()
        int sane()

    cdef cppclass GSLrng_t:
        GSLrng_t(unsigned)
        gsl_rng * get()

#Now, provide definitions for classes in classes.pyx
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
    cpdef append(self,popvec p)

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
    cpdef append(self,mpopvec p)

cdef class GSLrng:
    cdef GSLrng_t * thisptr

#Typedefs for convenience
ctypedef vector[popgenmut].iterator mcont_t_itr
ctypedef vector[mcont_t_itr] mut_container_t
ctypedef vector[gamete_t].iterator gcont_t_itr
ctypedef vector[diploid_t].iterator dipvector_t_itr
#vector of mutation counts (replaces KTfwd::mutation_base::n in fwdpp >= 0.4.4)
ctypedef vector[unsigned] mcounts_cont_t

##Define some low-level functions that may be useful for others
cdef struct popgen_mut_data:
    double pos,s,h
    unsigned n,g
    bint neutral

cdef struct gamete_data:
    vector[popgen_mut_data] neutral,selected
    unsigned n

cdef struct diploid_data:
    gamete_data chrom0,chrom1
    double g,e,w,sh0,sh1
    int n0,n1

#cdef popgen_mut_data get_mutation( const vector[popgenmut].iterator & ) nogil
#cdef gamete_data get_gamete( const vector[gamete_t].iterator & ) nogil
#cdef diploid_data get_diploid( const vector[diploid_t].iterator & itr ) nogil
cdef popgen_mut_data get_mutation( const popgenmut & m, size_t n) nogil
cdef gamete_data get_gamete( const gamete_t & g, const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil
cdef diploid_data get_diploid( const diploid_t & dip, const gcont_t & gametes, const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil

##Now, wrap the functions.
##To whatever extent possible, we avoid cdef externs in favor of Cython fxns based on cpp types.
##Many of the functions below rely on templates or other things that are too complex for Cython to handle at the moment
cdef extern from "sample.hpp" namespace "fwdpy" nogil:
    void get_sh( const vector[pair[double,string]] & ms_sample, const singlepop_t * pop, vector[double] * s,vector[double] * h, vector[double] * p, vector[double] * a)
    void get_sh( const vector[pair[double,string]] & samples, const metapop_t * pop, vector[double] * s, vector[double] * h, vector[double] * p, vector[double] * a)

cdef extern from "deps.hpp" namespace "fwdpy" nogil:
    vector[string] fwdpy_version()
    void fwdpy_citation()
    
cdef extern from "metapop.hpp" namespace "fwdpy" nogil:
    void re_init_mpop( metapop_t * mpop, const singlepop_t * pop)
    void copy_deme( metapop_t * mpop, const size_t i ) except +
    void remove_deme( metapop_t * mpop, const size_t i ) except +
    void merge_demes(metapop_t  * mpop, const size_t i, const size_t j) except +
    void split_deme(const gsl_rng * r, metapop_t * mpop, const size_t i, const unsigned N_new, const bint replacement ) except +
    void admix_demes(const gsl_rng * r, metapop_t * mpop, const size_t i, const size_t j, const double prop_i, const bint replacement) except +
    
cdef extern from "evolve_regions.hpp" namespace "fwdpy" nogil:
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
    

cdef extern from "sampler_sample_n.hpp" nogil:
    cdef struct detailed_deme_sample:
        sep_sample_t genotypes
        vector[pair[double,double]] sh

cdef extern from "sampler_selected_mut_tracker.hpp" nogil:
    cdef struct selected_mut_data:
        double pos
        double esize
        unsigned origin

cdef extern from "sampler_pop_properties.hpp" nogil:
    cdef struct qtrait_stats_cython:
        string stat
        double value
        unsigned generation

cdef extern from "allele_ages.hpp" nogil:
    cdef struct allele_age_data_t:
        double esize
        double max_freq
        double last_freq
        unsigned origin
        unsigned tlen

cdef extern from "allele_ages.hpp" namespace "fwdpy" nogil:
    vector[allele_age_data_t] allele_ages_details( const vector[pair[selected_mut_data,vector[double]]] & trajectories,
						   const double minfreq, const unsigned minsojourn ) except +
    
    vector[pair[selected_mut_data,vector[double]]] merge_trajectories_details( vector[pair[selected_mut_data,vector[double]]] traj1,
                                                                               const vector[pair[selected_mut_data,vector[double]]] & traj2 )
    
ctypedef unsigned uint
cdef extern from "evolve_regions_sampler.hpp" namespace "fwdpy" nogil:
    void evolve_regions_no_sampling_async(GSLrng_t * rng,
                                          vector[shared_ptr[singlepop_t]] * pops,
                                          const unsigned * Nvector,
                                          const size_t Nvector_len,
                                          const double mu_neutral,
                                          const double mu_selected,
                                          const double littler,
                                          const double f,
                                          const region_manager * rm,
                                          const char * fitness)

    vector[vector[pair[uint,detailed_deme_sample]]] evolve_regions_sample_async(GSLrng_t * rng,
                                                                        vector[shared_ptr[singlepop_t]] * pops,
                                                                        const unsigned * Nvector,
                                                                        const size_t Nvector_len,
                                                                        const double mu_neutral,
                                                                        const double mu_selected,
                                                                        const double littler,
                                                                        const double f,
                                                                        const int sample,
                                                                        const unsigned nsam,
                                                                        const region_manager * rm,
                                                                        const char * fitness)

    vector[vector[pair[selected_mut_data,vector[double]]]] evolve_regions_track_async(GSLrng_t * rng,
                                                                                      vector[shared_ptr[singlepop_t]] * pops,
                                                                                      const unsigned * Nvector,
                                                                                      const size_t Nvector_len,
                                                                                      const double mu_neutral,
                                                                                      const double mu_selected,
                                                                                      const double littler,
                                                                                      const double f,
                                                                                      const int sample,
                                                                                      const region_manager * rm,
                                                                                      const char * fitness)

