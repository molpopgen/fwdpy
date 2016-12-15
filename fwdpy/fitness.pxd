#Note: Cython cannot cimport typedefs, which is a bummer
from fwdpy.fwdpy cimport singlepop_t,multilocus_t,diploid_t
from fwdpy.fwdpp cimport popgenmut,gamete_base
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr

cdef extern from "<algorithm>" namespace "std" nogil:
    T max[T](T,T)

ctypedef gamete_base[void] gamete_t
ctypedef vector[gamete_t] gcont_t
ctypedef vector[popgenmut] mcont_t

cdef extern from "fwdpy_fitness.hpp" namespace "fwdpy" nogil:
    ctypedef void(*genotype_fitness_updater)(double &, const popgenmut &)
    ctypedef double(*fitness_function_finalizer)(double)
    ctypedef double(*haplotype_fitness_fxn_finalizer)(double,double)
    ctypedef double(*haplotype_fitness_fxn)(const gamete_t &, const mcont_t &)
    ctypedef double(*single_region_fitness_callback)(const diploid_t &, const gcont_t &, const mcont_t &)
    
    cdef cppclass site_dependent_fitness_wrapper:
        double operator()[DIPLOID,GAMETE_CONTAINER,
                          MUTATION_CONTAINER,
                          AA,Aa,FINAL](const DIPLOID &,
                                       const GAMETE_CONTAINER &,
                                       const MUTATION_CONTAINER &,
                                       const AA &, const Aa &,
                                       const FINAL &,
                                       const double &) const
    #stateless fitness object:
    cdef cppclass singlepop_fitness:
        singlepop_fitness()
        singlepop_fitness(single_region_fitness_callback)
        singlepop_fitness(genotype_fitness_updater Aa,
		          genotype_fitness_updater aa,
		          fitness_function_finalizer wfinal,
		          double starting_fitness)
        singlepop_fitness(haplotype_fitness_fxn h,
		          haplotype_fitness_fxn_finalizer f)
        void update(const singlepop_t *)

    cdef cppclass multilocus_fitness:
        multilocus_fitness()
        void update(const multilocus_t *)

    multilocus_fitness make_mloc_additive_fitness(double scaling)
    multilocus_fitness make_mloc_multiplicative_fitness(double scaling)
    multilocus_fitness make_mloc_additive_trait(double scaling)
    multilocus_fitness make_mloc_multiplicative_trait(double scaling)
    multilocus_fitness make_mloc_gbr_trait()
    multilocus_fitness make_mloc_power_mean_trait(const double SLd,const double MLd,
						  const vector[double] SLp,
						  const vector[double] MLp)
    ctypedef double(*mlocus_fitness_fxn)(const vector[diploid_t] &, const gcont_t &, const mcont_t &)
    multilocus_fitness make_mloc_custom_fitness(mlocus_fitness_fxn f)
    
#Helper functions for making custom fitness functions
cdef inline double return_w(double w) nogil:
    return max[double](0.0,w)

cdef inline double return_w_plus1(double w) nogil:
    return max[double](0.0,1.0+w)

cdef inline double return_trait_value(double w) nogil:
    return w

cdef inline double return_trait_value_minus1(double w) nogil:
    return w-1.0

cdef inline void het_additive_update(double & w, const popgenmut & m) nogil:
    (&w)[0] += m.s*m.h

cdef inline void hom_additive_update_2(double & w, const popgenmut & m) nogil:
    (&w)[0] += 2.0*m.s

cdef inline void hom_additive_update_1(double & w, const popgenmut & m) nogil:
    (&w)[0] += m.s

cdef inline genotype_fitness_updater choose_additive_hom_updater(int scaling) nogil:
    #Defaults to using a scaling of 2
    if scaling==1:
        return <genotype_fitness_updater>hom_additive_update_1
    return <genotype_fitness_updater>hom_additive_update_2
    
cdef inline void het_mult_update(double & w, const popgenmut & m) nogil:
    (&w)[0] *= (1.0+m.s*m.h)

cdef inline void hom_mult_update_2(double & w, const popgenmut & m) nogil:
    (&w)[0] *= (1.0+2.0*m.s)

cdef inline void hom_mult_update_1(double & w, const popgenmut & m) nogil:
    (&w)[0] *= (1.0+m.s)

cdef inline genotype_fitness_updater choose_mult_hom_updater(int scaling) nogil:
    #Defaults to using a scaling of 2
    if scaling==1:
        return <genotype_fitness_updater>hom_additive_update_1
    return <genotype_fitness_updater>hom_additive_update_2

cdef inline double sum_haplotype_effects(const gamete_t & g, const mcont_t & m) nogil:
    cdef size_t i=0,n=g.smutations.size()
    cdef double rv = 0.0
    while i<n:
        rv+=m[g.smutations[i]].s
        i+=1
    return rv
    
cdef class SpopFitness(object):
    """
    Base object for single-deme fitness functions
    """
    cdef unique_ptr[singlepop_fitness] wfxn
    
cdef class SpopAdditive(SpopFitness):
    pass

cdef class SpopMult(SpopFitness):
    pass

cdef class MlocusFitness(object):
    """
    Base object for multi-locus fitness functions
    """
    cdef multilocus_fitness wfxn

cdef class MlocusAdditive(MlocusFitness):
    pass

cdef class MlocusMult(MlocusFitness):
    pass

