from fwdpy.fwdpp cimport gamete_base,popgenmut
from fwdpy.fwdpy cimport diploid_t
from fwdpy.fitness cimport MlocusFitness
from libcpp.vector cimport vector

ctypedef gamete_base[void] gamete_t
ctypedef vector[gamete_t] gcont_t
ctypedef vector[popgenmut] mcont_t

cdef class MlocusAdditiveTrait(MlocusFitness):
    pass

cdef class MlocusGBRTrait(MlocusFitness):
    pass

cdef class MlocusMultTrait(MlocusFitness):
    pass

cdef class MlocusPowerMeanTrait(MlocusFitness):
    cdef vector[double] SLd,MLd
    cdef double SLp,MLp
