# distutils: language = c++
# distutils: sources = fwdpy/qtrait_mloc/qtrait_mloc_impl.cc

from cython.operator cimport dereference as deref,preincrement as inc
from libcpp.vector cimport vector
from fwdpy.fwdpp cimport shmodel
from fwdpy.fwdpy cimport *
from fwdpy.internal.internal cimport shwrappervec,shmodel
from fwdpy.fitness cimport multilocus_fitness
from fwdpy.fitness cimport make_mloc_gbr_trait,make_mloc_power_mean_trait,make_mloc_additive_trait,make_mloc_multiplicative_trait

cdef class MlocusGBRTrait(MlocusFitness):
    """
    The "gene-based recessive" (GBR) trait model of Thornton et al (2013) PLoS Genetics.

    The GBR model is applied to each region, and then the final value is the sum over regions.

    The fitness (or genetic value as it was used in that paper) is the geometric mean 
    of haplotype effect sizes.

    .. note:: Be really careful with this one!  Fitnesses are undefined if the sum of effect sizes on a haplotype is :math:`< 0:`.  The intended use case is to calculate a trait value under models with effect sizes :math:`>0`.
    """ 
    def __cinit__(self):
        self.wfxn=make_mloc_gbr_trait()

cdef class MlocusPowerMeanTrait(MlocusFitness):
    def __cinit__(self,double slp,double mlp,vector[double] & sld, vector[double] & mld):
        if sld.size() != 2:
            raise RuntimeError("sld must be of length 2")
        self.SLp = slp
        self.MLp = mlp
        self.SLd = sld
        self.MLd = mld
        self.wfxn = make_mloc_power_mean_trait(self.SLp,self.MLp,self.SLd,self.MLd)

cdef class MlocusAdditiveTrait(MlocusFitness):
    def __cinit__(self,double scaling=2.0):
        self.wfxn=make_mloc_additive_trait(scaling)

cdef class MlocusMultTrait(MlocusFitness):
    def __cinit__(self,double scaling = 2.0):
        self.wfxn=make_mloc_multiplicative_trait(scaling)

cdef extern from "qtrait_evolve_mlocus.hpp" namespace "fwdpy::qtrait" nogil:
    void evolve_qtrait_mloc_cpp( GSLrng_t * rng,
			         vector[shared_ptr[multilocus_t]] * pops,
			         vector[unique_ptr[sampler_base]] & samplers,
			         const unsigned * Nvector,
			         const size_t Nvector_length,
			         const vector[double] & neutral_mutation_rates,
			         const vector[double] & selected_mutation_rates,
                                 const vector[shmodel] & shmodels,
			         const vector[double] & within_region_rec_rates,
			         const vector[double] & between_region_rec_rates,
			         const double f,
			         const double sigmaE,
			         const double optimum,
			         const double VS,
                                 const int sample,
			         const multilocus_fitness & fitness ) except +

    void evolve_qtrait_mloc_regions_cpp(GSLrng_t *rng,
            vector[shared_ptr[multilocus_t]] *pops,
            vector[unique_ptr[sampler_base]] &samplers,
            const unsigned *Nvector, const size_t Nvector_length,
            const region_manager * rm,
            const vector[double] &between_region_rec_rates,
            const double f, const double sigmaE, const double optimum,
            const double VS, const int interval,
            const multilocus_fitness &fitness) except +
    
include "evolve_qtraits_mloc.pyx"
