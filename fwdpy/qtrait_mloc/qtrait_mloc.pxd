from fwdpy.fwdpp cimport gamete_base,popgenmut
from fwdpy.fwdpy cimport diploid_t,gamete_t,gcont_t,mcont_t
from fwdpy.fitness cimport MlocusFitness
from libcpp.vector cimport vector

cdef class MlocusAdditiveTrait(MlocusFitness):
    pass

cdef class MlocusGBRTrait(MlocusFitness):
    pass

cdef class MlocusMultTrait(MlocusFitness):
    pass

cdef class MlocusPowerMeanTrait(MlocusFitness):
    cdef vector[double] SLd,MLd
    cdef double SLp,MLp

