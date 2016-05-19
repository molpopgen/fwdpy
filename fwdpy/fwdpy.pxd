from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from libcpp.map cimport map

from fwdpy.internal.internal cimport *
from fwdpy.fwdpp cimport popgenmut,gamete_base
from fwdpy.gsl cimport gsl_rng
from fwdpy.structs cimport detailed_deme_sample,selected_mut_data,selected_mut_data_tidy,qtrait_stats_cython,allele_age_data_t,haplotype_matrix

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
        metapop_t(const singlepop_t &)
        unsigned generation
        vector[unsigned] Ns
        mcont_t mutations
        vector[unsigned] mcounts
        gcont_t gametes
        vector[dipvector_t] diploids
        vector[popgenmut] fixations
        vector[unsigned] fixation_times
        vector[unsigned] popsizes()
        int sane()
        int size()

    cdef cppclass multilocus_t:
        multilocus_t(unsigned,unsigned)
        multilocus_t(const multilocus_t &)
        unsigned generation
        unsigned N
        mcont_t mutations
        vector[unsigned] mcounts
        gcont_t gametes
        #This has different interpretation than for a metapop--see fwdpp dox
        vector[dipvector_t] diploids
        vector[popgenmut] fixations
        vector[unsigned] fixation_times
        int gen()
        int sane()
        int popsize()

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
        const gsl_rng * get()

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
    cpdef from_singlepop(self,singlepop p)

cdef class singlepop_mloc(poptype):
    cdef shared_ptr[multilocus_t] pop
    cpdef gen(self)
    cpdef popsize(self)
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

cdef class popvec_mloc(popcont):
    cdef vector[shared_ptr[multilocus_t]] pops
    cdef public object pypops
    cpdef size(self)
    cdef reset(self,const vector[shared_ptr[multilocus_t]] newpops)
    cpdef append(self,popvec_mloc p)


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
    uint16_t label

cdef struct gamete_data:
    vector[popgen_mut_data] neutral,selected
    unsigned n

cdef struct diploid_data:
    gamete_data chrom0,chrom1
    double g,e,w,sh0,sh1
    int n0,n1

cdef struct diploid_mloc_data:
    vector[gamete_data] chrom0,chrom1
    double g,e,w
    vector[double] sh0,sh1
    vector[int] n0,n1

#cdef popgen_mut_data get_mutation( const vector[popgenmut].iterator & ) nogil
#cdef gamete_data get_gamete( const vector[gamete_t].iterator & ) nogil
#cdef diploid_data get_diploid( const vector[diploid_t].iterator & itr ) nogil
cdef popgen_mut_data get_mutation( const popgenmut & m, size_t n) nogil
cdef gamete_data get_gamete( const gamete_t & g, const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil
cdef diploid_data get_diploid( const diploid_t & dip, const gcont_t & gametes, const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil
cdef diploid_mloc_data get_diploid_mloc( const dipvector_t & dip, const gcont_t & gametes, const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil

##Now, wrap the functions.
##To whatever extent possible, we avoid cdef externs in favor of Cython fxns based on cpp types.
##Many of the functions below rely on templates or other things that are too complex for Cython to handle at the moment
cdef extern from "sample.hpp" namespace "fwdpy" nogil:
    void get_sh( const vector[pair[double,string]] & ms_sample,
                 const vector[popgenmut] & mutations,
                 const vector[unsigned] & mcounts,
                 const unsigned & ttlN,
                 const unsigned & generation,
                 vector[double] * s,vector[double] * h, vector[double] * p, vector[double] * a, vector[uint16_t] * l)

cdef extern from "deps.hpp" namespace "fwdpy" nogil:
    vector[string] fwdpy_version()
    void fwdpy_citation()

cdef extern from "sampler_selected_mut_tracker.hpp" namespace "fwdpy" nogil:
    vector[selected_mut_data_tidy] tidy_trajectory_info( const vector[pair[selected_mut_data,vector[double]]] & trajectories,
                                                         const unsigned min_sojourn, const double min_freq);

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
                                          const char * fitness) except +

    vector[vector[pair[uint,detailed_deme_sample]]] evolve_regions_sample_async(GSLrng_t * rng_evolve,GSLrng_t * rng_sample,
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
                                                                        const char * fitness) except +

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
                                                                                      const char * fitness) except +

cdef extern from "sampling_wrappers.hpp" namespace "fwdpy" nogil:
    sample_t sample_single[POPTYPE](gsl_rng * r,const POPTYPE & p, const unsigned nsam, const bool removeFixed)  except +
    sep_sample_t sample_sep_single[POPTYPE](gsl_rng * r,const POPTYPE & p, const unsigned nsam, const bool removeFixed)  except +
    vector[sample_t] sample_single_mloc[POPTYPE](gsl_rng * r,const POPTYPE & p, const unsigned nsam, const bool removeFixed)  except +
    vector[sep_sample_t] sample_sep_single_mloc[POPTYPE](gsl_rng * r,const POPTYPE & p, const unsigned nsam, const bool removeFixed)  except +

cdef extern from "haplotype_matrix.hpp" namespace "fwdpy" nogil:
    haplotype_matrix make_haplotype_matrix(const singlepop_t * pop, const vector[size_t] & diploids) except +
    haplotype_matrix make_haplotype_matrix(const metapop_t * pop, const vector[size_t] & diploids,const size_t deme) except +
    haplotype_matrix make_haplotype_matrix(const multilocus_t * pop, const vector[size_t] & diploids) except +
