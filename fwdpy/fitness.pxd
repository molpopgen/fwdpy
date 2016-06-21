#Note: Cython cannot cimport typedefs, which is a bummer
from fwdpy.fwdpy cimport singlepop_t,multilocus_t#,gamete_t,mcont_t
from fwdpy.fwdpp cimport popgenmut,gamete_base
from libcpp.vector cimport vector

ctypedef gamete_base[void] gamete_t
ctypedef vector[gamete_t] gcont_t
ctypedef vector[popgenmut] mcont_t

cdef extern from "fwdpy_fitness.hpp" namespace "fwdpy" nogil:
    cdef cppclass singlepop_fitness:
        singlepop_fitness()
        void update(const singlepop_t *)

    cdef cppclass multilocus_fitness:
        multilocus_fitness()
        void update(const multilocus_t *)

    ctypedef void(*genotype_fitness_updater)(double &, const popgenmut &)
    ctypedef double(*fitness_function_finalizer)(double)
    ctypedef double(*haplotype_fitness_fxn_finalizer)(double,double)
    ctypedef double(*haplotype_fitness_fxn)(const gamete_t &, const mcont_t &)
    
    singlepop_fitness make_additive_fitness(double scaling)
    singlepop_fitness make_multiplicative_fitness(double scaling)
    singlepop_fitness make_gbr_fitness()
    singlepop_fitness make_custom_fitness(genotype_fitness_updater Aa,
					  genotype_fitness_updater aa,
					  fitness_function_finalizer wfinal,
					  double starting_fitness)
    singlepop_fitness make_custom_haplotype_fitness(haplotype_fitness_fxn h,
                                                    haplotype_fitness_fxn_finalizer f)
                                                    
    
    multilocus_fitness make_mloc_additive_fitness(double scaling)
    multilocus_fitness make_mloc_multiplicative_fitness(double scaling)
    multilocus_fitness make_mloc_gbr_fitness()

#Helper functions for making custom fitness functions
cdef inline double return_w(const double w):
    return max(0.0,w)

cdef inline double return_w_plus1(const double w):
    return max(0.0,1.0+w)

cdef class singlepopFitness(object):
    """
    Base object for single-deme fitness functions
    """
    cdef singlepop_fitness wfxn
    
cdef class singlepopAdditive(singlepopFitness):
    pass

cdef class singlepopGBR(singlepopFitness):
    pass

cdef class singlepopMult(singlepopFitness):
    pass

cdef class multilocusFitness(object):
    """
    Base object for multi-locus fitness functions
    """
    cdef multilocus_fitness wfxn

cdef class multilocusAdditive(multilocusFitness):
    pass

cdef class multilocusGBR(multilocusFitness):
    pass

cdef class multilocusMult(multilocusFitness):
    pass
