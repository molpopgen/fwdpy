from fwdpy.fwdpp cimport popgenmut,gamete_base
from fwdpy.fitness cimport SpopFitness
from fwdpy.fwdpy cimport singlepop_t,sampler_base,singlepop_fitness,GSLrng_t
from fwdpy.internal.internal cimport shwrappervec,region_manager
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr,unique_ptr
cdef class SpopGBRTrait(SpopFitness):
    pass

cdef class SpopAdditiveTrait(SpopFitness):
    pass

cdef class SpopMultTrait(SpopFitness):
    pass

cdef extern from "qtrait_evolve.hpp" namespace "fwdpy::qtrait" nogil:
    void evolve_regions_qtrait_cpp(GSLrng_t * rng,
				   vector[shared_ptr[singlepop_t]] & pops,
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
