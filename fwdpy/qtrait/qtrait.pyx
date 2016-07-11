# distutils: language = c++
# distutils: sources = fwdpy/qtrait/qtrait_impl.cc fwdpy/qtrait/ew2010.cc

from cython.operator cimport dereference as deref,preincrement as inc
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.utility cimport pair
from libcpp.string cimport string
from fwdpy.fwdpy cimport *
from fwdpy.internal.internal cimport shwrappervec
from fwdpy.fitness cimport *
from libc.math cimport sqrt
import fwdpy.internal as internal

cdef inline double geomean(double a,double b):
    return sqrt(a*b)

cdef class SpopGBRTrait(SpopFitness):
    """
    The "gene-based recessive" (GBR) model of Thornton et al (2013) PLoS Genetics.

    The genetic value is the geometric mean of haplotype effect sizes.

    .. note:: Be really careful with this one!  Fitnesses are undefined if the sum of effect sizes on a haplotype is :math:`< 0:`.  The intended use case is to calculate a trait value under models with effect sizes :math:`>0`.
    """ 
    def __cinit__(self):
        self.wfxn = singlepop_fitness(<haplotype_fitness_fxn>sum_haplotype_effects,
                                      <haplotype_fitness_fxn_finalizer>geomean)

cdef class SpopAdditiveTrait(SpopFitness):
    def __cinit__(self,int scaling = 2):
        self.wfxn = singlepop_fitness(<genotype_fitness_updater>het_additive_update,
                                      choose_additive_hom_updater(scaling),
                                      <fitness_function_finalizer>return_trait_value,
                                      0.0)

cdef class SpopMultTrait(SpopFitness):
    def __cinit__(self,int scaling = 2):
        self.wfxn = singlepop_fitness(<genotype_fitness_updater>het_mult_update,
                                      choose_mult_hom_updater(scaling),
                                      <fitness_function_finalizer>return_trait_value_minus1,
                                      1.0)

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

cdef extern from "qtrait_evolve.hpp" namespace "fwdpy::qtrait" nogil:
    void evolve_regions_qtrait_cpp(GSLrng_t * rng,
				   vector[shared_ptr[singlepop_t]] * pops,
				   vector[unique_ptr[sampler_base]] & samplers,
				   const unsigned * Nvector,
				   const size_t Nvector_length,
				   const double neutral,
				   const double selected,
				   const double recrate,
				   const double f,
				   const double sigmaE,
				   const double optimum,
				   const double VS,
				   const int interval,
				   const region_manager * rm,
				   const singlepop_fitness & fitness) except +

include "evolve_qtraits.pyx"
include "ew2010.pyx"
include "misc.pyx"



