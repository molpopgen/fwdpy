from fwdpy.fwdpp cimport popgenmut,gamete_base
from fwdpy.fitness cimport SpopFitness

cdef class SpopGBRTrait(SpopFitness):
    pass

cdef class SpopAdditiveTrait(SpopFitness):
    pass

cdef class SpopMultTrait(SpopFitness):
    pass
