from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr,unique_ptr

from libcpp.map cimport map

from fwdpy.internal.internal cimport *
from fwdpy.fwdpp cimport popgenmut,gamete_base
from fwdpy.cpp cimport hash
from libcpp.unordered_set cimport unordered_set
from cython_gsl cimport gsl_rng
from fwdpy.structs cimport selected_mut_data,selected_mut_data_tidy,qtrait_stats_cython,allele_age_data_t,VAcum,popsample_details
from fwdpy.fitness cimport singlepop_fitness

##Create hooks to C++ types
ctypedef vector[unsigned] ucont_t

#Wrap the classes:
cdef extern from "types.hpp" namespace "fwdpy" nogil:
    # "Standard" popgen types
    ctypedef gamete_base[void] gamete_t
    ctypedef vector[gamete_t] gcont_t
    ctypedef vector[popgenmut] mcont_t
    ctypedef unordered_set[double,equal_eps] lookup_t
    
    cdef cppclass diploid_t:
        size_t first,second,label;
        double g,e,w

    ctypedef vector[diploid_t] dipvector_t

    cdef cppclass singlepop_t:
        singlepop_t(unsigned)
        unsigned N
        unsigned generation
        mcont_t mutations
        ucont_t mcounts
        gcont_t gametes
        dipvector_t diploids
        mcont_t fixations
        ucont_t fixation_times
        lookup_t mut_lookup
        unsigned gen()
        unsigned popsize()
        int sane()
        string serialize() const
        void deserialize(const string &)
        int tofile(const char *,bint)
        void fromfile(const char *,size_t)
        void clear()

    cdef cppclass metapop_t:
        metapop_t(ucont_t)
        metapop_t(const singlepop_t &)
        unsigned generation
        ucont_t Ns
        mcont_t mutations
        ucont_t mcounts
        gcont_t gametes
        vector[dipvector_t] diploids
        mcont_t fixations
        ucont_t fixation_times
        lookup_t mut_lookup
        ucont_t popsizes()
        int sane()
        int size()
        string serialize() const
        void deserialize(const string &)
        int tofile(const char *,bint)
        void fromfile(const char *,size_t)
        void clear()

    cdef cppclass multilocus_t:
        multilocus_t(unsigned,unsigned)
        multilocus_t(const multilocus_t &)
        unsigned generation
        unsigned N
        mcont_t mutations
        ucont_t mcounts
        gcont_t gametes
        #This has different interpretation than for a metapop--see fwdpp dox
        vector[dipvector_t] diploids
        mcont_t fixations
        ucont_t fixation_times
        lookup_t mut_lookup
        int gen()
        int sane()
        int popsize()
        string serialize() const
        void deserialize(const string &)
        int tofile(const char *,bint)
        void fromfile(const char *,size_t)
        void clear()

    # Types based around KTfwd::generalmut_vec
    ctypedef vector[generalmut_vec] mlist_gm_vec_t

    cdef cppclass singlepop_gm_vec_t:
        singlepop_gm_vec_t(unsigned)
        const unsigned N
        const unsigned generation
        mlist_gm_vec_t mutations
        ucont_t mcounts
        gcont_t gametes
        vector[diploid_t] diploids
        vector[generalmut_vec] fixations
        ucont_t fixation_times
        lookup_t mut_lookup
        unsigned gen()
        unsigned popsize()
        int sane()
        string serialize() const
        void deserialize(const string &)
        int tofile(const char *,bint)
        void fromfile(const char *,size_t)
        void clear()

    cdef cppclass GSLrng_t:
        GSLrng_t(unsigned)
        const gsl_rng * get()

#Now, provide definitions for classes in classes.pyx
cdef class PopType(object):
    """
    Empty base class for a population object.

    Example derived types include :class:`fwdpy.fwdpy.singlepop` and :class:`fwdpy.fwdpy.metapop`
    """
    pass

cdef class Spop(PopType):
    cdef shared_ptr[singlepop_t] pop
    cpdef gen(self)
    cpdef popsize(self)
    cpdef sane(self)

cdef class MetaPop(PopType):
    cdef shared_ptr[metapop_t] mpop
    cpdef gen(self)
    cpdef popsizes(self)
    cpdef sane(self)
    cpdef from_Spop(self,Spop p)

cdef class MlocusPop(PopType):
    cdef shared_ptr[multilocus_t] pop
    cpdef gen(self)
    cpdef popsize(self)
    cpdef sane(self)

cdef class SpopGenMut(PopType):
    cdef shared_ptr[singlepop_gm_vec_t] pop
    cpdef gen(self)
    cpdef popsize(self)
    cpdef sane(self)

cdef class PopVec(object):
    """
    Empty base class for containers of population objects.

    Example derived types include :class:`fwdpy.fwdpy.SpopVec` and :class:`fwdpy.fwdpy.MetaPopVec`
    """
    pass

cdef class SpopVec(PopVec):
    cdef vector[shared_ptr[singlepop_t]] pops
    cpdef size(self)
    cdef reset(self,const vector[shared_ptr[singlepop_t]] newpops)
    cpdef append(self,SpopVec p)

cdef class SpopGenMutVec(PopVec):
    cdef vector[shared_ptr[singlepop_gm_vec_t]] pops
    cpdef size(self)
    cdef reset(self,const vector[shared_ptr[singlepop_gm_vec_t]] newpops)

cdef class MetaPopVec(PopVec):
    cdef vector[shared_ptr[metapop_t]] mpops
    cpdef size(self)
    cdef reset(self,const vector[shared_ptr[metapop_t]]  & mpops)
    cpdef append(self,MetaPopVec p)

cdef class MlocusPopVec(PopVec):
    cdef vector[shared_ptr[multilocus_t]] pops
    cpdef size(self)
    cdef reset(self,const vector[shared_ptr[multilocus_t]] newpops)
    cpdef append(self,MlocusPopVec p)


cdef class GSLrng:
    cdef GSLrng_t * thisptr

#Functions relating to built-in temporal sampling features

cdef extern from "sampler_base.hpp" namespace "fwdpy" nogil:
    cdef cppclass sampler_base:
        pass
    void clear_samplers(vector[unique_ptr[sampler_base]] & v)

    #This is a template for a custom temporal sampler
    #Each generation, a "generation_t" is recorded, and
    #are kept in a container alled FINALT, which is returned,
    #and represents a time series of observations.
    #The best policy is to keep FINALT as something coercible
    #to Python via Cython directly.
    cdef cppclass custom_sampler[FINALT](sampler_base):
        FINALT final()
        FINALT f
        #To make a valid sampler, you must construct and object with a callback.
        #Below are the constructors for singlepop_t, multilocus_t, and metapop_t,
        #respectively.
        custom_sampler(void(*)(const singlepop_t *, const unsigned, FINALT &))
        custom_sampler(void(*)(const multilocus_t *, const unsigned, FINALT &))
        custom_sampler(void(*)(const metapop_t *, const unsigned, FINALT &))
        #You may register additional callbacks as needed
        void register_callback(void(*)(const singlepop_t *, const unsigned, FINALT &) )
        void register_callback(void(*)(const multilocus_t *, const unsigned, FINALT &))
        void register_callback(void(*)(const metapop_t *, const unsigned, FINALT &))
        #Typically, you won't need to register a cleanup callback,
        #but the interface is provided anyways
        void register_callback(void(*)(FINALT &))

    #Template for custom temporal sampler that need additional data
    #The DATAT may be arbitrary--it could contain parameters used
    #for each call when sampling, or reusable buffers for intermediate
    #compuations, etc.
    cdef cppclass custom_sampler_data[FINALT,DATAT](sampler_base):
        DATAT data
        FINALT f
        FINALT final()
        #To make a valid sampler, you must construct and object with a callback.
        #Below are the constructors for singlepop_t, multilocus_t, and metapop_t,
        #respectively.  Each constructor also requires a DATAT, and DATAT must be copyable.
        custom_sampler_data(void(*)(const singlepop_t *, const unsigned, FINALT &, DATAT &), const DATAT &)
        custom_sampler_data(void(*)(const multilocus_t *, const unsigned, FINALT &, DATAT &), const DATAT &)
        custom_sampler_data(void(*)(const metapop_t *, const unsigned, FINALT &, DATAT & ), const DATAT &)
        #You may register additional callbacks as needed
        void register_callback(void(*)(const singlepop_t *, const unsigned, FINALT &, DATAT &))
        void register_callback(void(*)(const multilocus_t *, const unsigned, FINALT &, DATAT &))
        void register_callback(void(*)(const metapop_t *, const unsigned, FINALT &, DATAT &))
        #if DATAT can get large, you may register a callback to clean it up
        #when evolution functions exit
        void register_callback(void(*)(FINALT &,DATAT &))
                          
    void apply_sampler_cpp[T](const vector[shared_ptr[T]] & popvec,
			      const vector[unique_ptr[sampler_base]] & samplers)

    void apply_sampler_single_cpp[T](const T *pop,const vector[unique_ptr[sampler_base]] & samplers)

cdef extern from "sampler_no_sampling.hpp" namespace "fwdpy" nogil:
    cdef cppclass no_sampling(sampler_base):
        no_sampling()

cdef extern from "sampler_pop_properties.hpp" namespace "fwdpy" nogil:
    cdef cppclass pop_properties(sampler_base):
        pop_properties(double optimum)
        vector[qtrait_stats_cython] final() const

cdef extern from "sampler_additive_variance.hpp" namespace "fwdpy" nogil:
    cdef cppclass additive_variance(sampler_base):
        additive_variance()
        vector[VAcum] final()

ctypedef vector[pair[sep_sample_t,popsample_details]] popSampleData

cdef extern from "sampler_sample_n.hpp" namespace "fwdpy" nogil:
    cdef cppclass sample_n(sampler_base):
        sample_n(unsigned, const gsl_rng * r,
                const string & nfile,const string & sfile,
                bint removeFixed, bint recordSamples, bint recordDetails,
                const vector[pair[double,double]] & boundaries,const bint append)
        popSampleData rv
        popSampleData final() const

#The following typedefs help us with the
#frequency tracker API.
ctypedef pair[uint,double] genfreqPair
ctypedef vector[pair[selected_mut_data,vector[genfreqPair]]] freqTraj

cdef extern from "sampler_selected_mut_tracker.hpp" namespace "fwdpy" nogil:
    cdef cppclass selected_mut_tracker(sampler_base):
        selected_mut_tracker()
        freqTraj final() const

#Extension classes for temporal sampling
cdef class TemporalSampler:
    """
    Base class representing containers of functions
    to be applied to simluated populations at regular intervals.
    """
    cdef vector[unique_ptr[sampler_base]] vec
    cpdef size_t size(self)

cdef class NothingSampler(TemporalSampler):
    pass

cdef class QtraitStatsSampler(TemporalSampler):
    pass

cdef class PopSampler(TemporalSampler):
    pass

cdef class VASampler(TemporalSampler):
    pass

cdef class FreqSampler(TemporalSampler):
    pass





#Typedefs for convenience
ctypedef vector[popgenmut].iterator mcont_t_itr
ctypedef vector[mcont_t_itr] mut_container_t
ctypedef vector[gamete_t].iterator gcont_t_itr
ctypedef vector[diploid_t].iterator dipvector_t_itr
#vector of mutation counts (replaces KTfwd::mutation_base::n in fwdpp >= 0.4.4)
ctypedef ucont_t mcounts_cont_t

##Define some low-level functions that may be useful for others
cdef struct popgen_mut_data:
    double pos,s,h
    unsigned n,g,ftime
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

cdef popgen_mut_data get_mutation( const popgenmut & m, size_t n) nogil
cdef gamete_data get_gamete( const gamete_t & g, const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil
cdef diploid_data get_diploid( const diploid_t & dip, const gcont_t & gametes, const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil
cdef diploid_mloc_data get_diploid_mloc( const dipvector_t & dip, const gcont_t & gametes, const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil

##Now, wrap the functions.
##To whatever extent possible, we avoid cdef externs in favor of Cython fxns based on cpp types.
##Many of the functions below rely on templates or other things that are too complex for Cython to handle at the moment
cdef extern from "sample.hpp" namespace "fwdpy" nogil:
    popsample_details get_sh( const vector[pair[double,string]] & ms_sample,
                 const mcont_t & mutations,
                 const mcont_t & fixations,
                 const ucont_t & fixation_times,
                 const ucont_t & mcounts,
                 const unsigned & ttlN,
                 const unsigned & generation,
                 const unsigned & locusID)

cdef extern from "deps.hpp" namespace "fwdpy" nogil:
    vector[string] fwdpy_version()
    void fwdpy_citation()


cdef extern from "sampler_selected_mut_tracker.hpp" namespace "fwdpy" nogil:
    vector[selected_mut_data_tidy] tidy_trajectory_info( const freqTraj & trajectories,
                                                         const unsigned min_sojourn, const double min_freq,
                                                         const unsigned remove_gone_before,
                                                         const unsigned remove_arose_after)
cdef extern from "allele_ages.hpp" namespace "fwdpy" nogil:
    vector[allele_age_data_t] allele_ages_details( const freqTraj & trajectories,
						   const double minfreq, const unsigned minsojourn ) except +

    freqTraj merge_trajectories_details( const freqTraj & traj1, const freqTraj & traj2 )

ctypedef unsigned uint
cdef extern from "evolve_regions_sampler.hpp" namespace "fwdpy" nogil:
    void evolve_regions_sampler_cpp( GSLrng_t * rng,
				     vector[shared_ptr[singlepop_t]] & pops,
				     vector[unique_ptr[sampler_base]] & samplers,
				     const unsigned * Nvector,
				     const size_t Nvector_length,
				     const double mu_neutral,
				     const double mu_selected,
				     const double littler,
				     const double f,
				     const int sample,
				     const region_manager * rm,
				     const singlepop_fitness & fitness) except +


cdef extern from "sampling_wrappers.hpp" namespace "fwdpy" nogil:
    sample_t sample_single[POPTYPE](gsl_rng * r,const POPTYPE & p, const unsigned nsam, const bool removeFixed)  except +
    sep_sample_t sample_sep_single[POPTYPE](gsl_rng * r,const POPTYPE & p, const unsigned nsam, const bool removeFixed)  except +
    vector[sample_t] sample_single_mloc[POPTYPE](gsl_rng * r,const POPTYPE & p, const unsigned nsam, const bool
            removeFixed, const vector[pair[double,double]] &)  except +
    vector[sep_sample_t] sample_sep_single_mloc[POPTYPE](gsl_rng * r,const POPTYPE & p, const unsigned nsam, const bool
            removeFixed, const vector[pair[double,double]] &)  except +

cdef extern from "fwdpy_add_mutations.hpp" namespace "fwdpy" nogil:
    size_t add_mutation_cpp(singlepop_t * pop,
                            const vector[size_t] & indlist,
                            const vector[short] & clist,
			    const double pos,
			    const double s,
			    const double h) except +
