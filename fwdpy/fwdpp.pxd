"""
This file exposes the fwdpp library to Cython.

It is currently not possible (and not necessarily even desirable) to expose all of fwdpp.

Cython currently doesn't seem to support all of the possible type of C++ template declaration necessary,
making this file a work in progress.

It is possible that some of these issues can be addresses via further evolution of the fwdpp 'sugar' layer,
simplifying the API considerably over the "raw" fwdpp API.
"""

from libcpp cimport bool
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector
from cython_gsl cimport gsl_rng,gsl_ran_discrete_t
from libc.stdint cimport uint16_t,uint32_t

##We will expose some low-level types from fwdpp:
cdef extern from "fwdpp/forward_types.hpp" namespace "KTfwd" nogil:
    cdef cppclass mutation_base:
        double pos
        uint16_t xtra
        bool neutral

    cdef cppclass gamete_base[T=*]:
        unsigned n
        vector[uint32_t] mutations
        vector[uint32_t] smutations

cdef extern from "fwdpp/sugar/popgenmut.hpp" namespace "KTfwd" nogil:
    cdef cppclass popgenmut(mutation_base):
        unsigned g
        double s
        double h

#Wrappers for fwdpp fitness models.  Call operators for custom diploids are exposed
cdef extern from "fwdpp/fitness_models.hpp" namespace "KTfwd" nogil:
    cdef cppclass site_dependent_fitness:
        #Hint: The final double is a starting fitness value.
        double operator()[DIPLOID,GAMETE_CONTAINER,
                          MUTATION_CONTAINER,AA,Aa](const DIPLOID &,
                                                    const GAMETE_CONTAINER &,
                                                    const MUTATION_CONTAINER &,
                                                    const AA &, const Aa &,
                                                    const double &) const

    cdef cppclass additive_diploid:
        #Wrapper around site_depenent fitness for simple additive model w = 1+sum(effects over loc).
        #The final double is the "scaling" term (see fwdpp docs)
        double operator()[DIPLOID,GAMETE_CONTAINER,
                          MUTATION_CONTAINER](const DIPLOID &,
                                              const GAMETE_CONTAINER &,
                                              const MUTATION_CONTAINER &,
                                              const double &) const

    cdef cppclass multiplicative_diploid:
        #Wrapper around site_depenent fitness for simple multiplicative model.
        #The final double is the "scaling" term (see fwdpp docs)
        double operator()[DIPLOID,GAMETE_CONTAINER,
                          MUTATION_CONTAINER](const DIPLOID &,
                                              const GAMETE_CONTAINER &,
                                              const MUTATION_CONTAINER &,
                                              const double &) const
#For hashing/lookup tables:
cdef extern from "fwdpp/fwd_functional.hpp" namespace "KTfwd" nogil:
    cdef cppclass equal_eps:
        bool operator()[T](const T & lhs, const T & rhs) const

#fwdpp's debug functions
cdef extern from "fwdpp/debug.hpp" namespace "KTfwd" nogil:
    bool check_sum[GAMETE_CONTAINER](const GAMETE_CONTAINER & gc, const unsigned twoN)
    bool gamete_is_sorted_n[GAMETE,MUTATION_CONTAINER](const GAMETE & g, const MUTATION_CONTAINER & mc)
    bool gamete_is_sorted_s[GAMETE,MUTATION_CONTAINER](const GAMETE & g, const MUTATION_CONTAINER & mc)
    bool gamete_data_sane[GAMETE_CONTAINER,MUTATION_CONTAINER](const GAMETE_CONTAINER & gc, const MUTATION_CONTAINER & mc, const vector[unsigned] & mcounts)
    bool popdata_sane[DIPLOID_CONTAINER,GAMETE_CONTAINER,MUTATION_CONTAINER](const DIPLOID_CONTAINER & dc, const GAMETE_CONTAINER & gc, const MUTATION_CONTAINER & mc, const vector[unsigned] & mcounts)

#fwdpp's "syntactic sugar" layer makes it easier to develop simulations
#by wrapping lower-level functions from main library
cdef extern from "fwdpp/sugar/generalmut.hpp" namespace "KTfwd" nogil:
    cdef cppclass generalmut_vec(mutation_base):
        vector[double] s
        vector[double] h
        unsigned g        

##See wrappers in fwdpy.pyx to functions in fwdpy's sampling_wrappers.hpp for thin wrappers to single-dem samplers
cdef extern from "fwdpp/sugar/sampling.hpp" namespace "KTfwd" nogil:
    ctypedef vector[pair[double,string]] sample_t
    ctypedef pair[sample_t,sample_t] sep_sample_t

    #Take random samples from a population
    sample_t sample[POPTYPE](const gsl_rng *,const POPTYPE &, 
                            const unsigned nsam, const bint removeFixed)
    sample_t sample_single"KTfwd::sample"[POPTYPE](const gsl_rng *,
                            const POPTYPE &, const unsigned nsam, const bint removeFixed)
    sep_sample_t sample_separate_single"KTfwd::sample_separate"[POPTYPE](const gsl_rng *,
                            const POPTYPE &, const unsigned nsam, const bint removeFixed)
    vector[sample_t] sample_mlocus"KTfwd::sample"[POPTYPE](const gsl_rng *,
                            const POPTYPE &, const unsigned nsam, const bint
                            removeFixed,const vector[pair[double,double]] &) except+
    vector[sep_sample_t] sample_separate_mlocus"KTfwd::sample_separate"[POPTYPE](const gsl_rng *,
                            const POPTYPE &, const unsigned
                            nsam, const bint removeFixed,const vector[pair[double,double]] &) except+
    sample_t sample_meta"KTfwd::sample"[POPTYPE](const gsl_rng *,
                            const POPTYPE &,const unsigned deme,
                            const unsigned nsam , const bool removeFixed) except+
    sep_sample_t sample_separate_meta"KTfwd::sample_separate"[POPTYPE](const gsl_rng *,
                            const POPTYPE &,const unsigned
                            deme , const unsigned nsam , const bool removeFixed) except+

    #Take specific samples from the population
    sample_t sample_single_ind"KTfwd::sample"[POPTYPE](const POPTYPE &,
                            const vector[size_t] &,
                            const bint removeFixed)
    sep_sample_t sample_separate_single_ind"KTfwd::sample_separate"[POPTYPE](const POPTYPE &,
                            const vector[size_t] &,
                            const bint removeFixed)

cdef extern from "fwdpp/sugar/change_neutral.hpp" namespace "KTfwd" nogil:
    void change_neutral[POPTYPE](POPTYPE & p, const size_t mut_index)

cdef extern from "fwdpp/sugar/matrix.hpp" namespace "KTfwd" nogil:
    cdef struct data_matrix:
        vector[char] neutral,selected
        vector[double] neutral_positions,selected_positions,neutral_popfreq,selected_popfreq
        size_t nrow

    pair[vector[pair[size_t,uint]],vector[pair[size_t,uint]]] mutation_keys[POPTYPE](const POPTYPE & pop,
            const vector[size_t] & inividuals,
            const bint include_neutral, 
            const bint include_selected, const size_t deme) except +

    data_matrix fwdpp_genotype_matrix"KTfwd::genotype_matrix"[POPTYPE](const POPTYPE & pop,
            const vector[size_t] & inividuals,
            const vector[pair[size_t,uint]] & neutral_keys,
            const vector[pair[size_t,uint]] & selected_keys, const size_t deme) except +

    data_matrix fwdpp_haplotype_matrix"KTfwd::haplotype_matrix"[POPTYPE](const POPTYPE & pop,
            const vector[size_t] & inividuals,
            const vector[pair[size_t,uint]] & neutral_keys,
            const vector[pair[size_t,uint]] & selected_keys, const size_t deme) except +

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

from libcpp.memory cimport unique_ptr
#These are "internal" fwdpp types

cdef extern from "fwdpp/internal/gsl_discrete.hpp" namespace "fwdpp::internal" nogil:
    #This is a wrapper around gsl_ran_discrete_free.
    #It is used as the deleter type for a unique_ptr
    cdef cppclass gsl_ran_discrete_t_deleter:
        pass

    #Smart-ptr wrapper around gsl_ran_discrete_t that uses custom deleter
    ctypedef unique_ptr[gsl_ran_discrete_t,gsl_ran_discrete_t_deleter] gsl_ran_discrete_t_ptr 
