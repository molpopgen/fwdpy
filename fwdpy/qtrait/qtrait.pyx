# distutils: language = c++
# distutils: sources = fwdpy/qtrait/qtrait_impl.cc fwdpy/qtrait/ew2010.cc

from cython.operator cimport dereference as deref,preincrement as inc
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.utility cimport pair
from fwdpy.fwdpy cimport *
from fwdpy.fitness cimport *
from libc.math cimport sqrt
import fwdpy.internal as internal

cdef inline double geomean(double a,double b) nogil:
    return sqrt(a*b)

cdef class SpopGBRTrait(SpopFitness):
    """
    The "gene-based recessive" (GBR) model of Thornton et al (2013) PLoS Genetics.

    The genetic value is the geometric mean of haplotype effect sizes.

    .. note:: Be really careful with this one!  Fitnesses are undefined if the sum of effect sizes on a haplotype is :math:`< 0:`.  The intended use case is to calculate a trait value under models with effect sizes :math:`>0`.
    """ 
    def __cinit__(self):
        self.wfxn = unique_ptr[singlepop_fitness](new singlepop_fitness(sum_haplotype_effects,
                                                                      	geomean))

cdef class SpopAdditiveTrait(SpopFitness):
    def __cinit__(self,int scaling = 2):
        self.wfxn = unique_ptr[singlepop_fitness](new singlepop_fitness(het_additive_update,
                                                                        choose_additive_hom_updater(scaling),
                                                                        return_trait_value,
                                                                        0.0))

cdef class SpopMultTrait(SpopFitness):
    def __cinit__(self,int scaling = 2):
        self.wfxn = unique_ptr[singlepop_fitness](new singlepop_fitness(het_mult_update,
                                                                        choose_mult_hom_updater(scaling),
                                                                        return_trait_value_minus1,
                                                                        1.0))

cdef extern from "qtrait_evolve_rules.hpp" namespace "fwdpy::qtrait" nogil:
    cdef cppclass qtrait_model_rules:
        qtrait_model_rules(const double & sigmaE, const double & optimum, const double & VS, const unsigned maxN)

cdef extern from "qtrait_pleiotropic.hpp" namespace "fwdpy::qtrait" nogil:    
    cdef struct ew_mut_details:
        double s
        double e
        double p
        
    map[double,ew_mut_details] ew2010_assign_effects(GSLrng_t * rng, const singlepop_t * pop, const double tau, const double sigma) except +
    vector[double] ew2010_traits_cpp(const singlepop_t * pop, const map[double,ew_mut_details] & effects) except +

include "evolve_qtraits.pyx"
include "ew2010.pyx"
include "misc.pyx"



