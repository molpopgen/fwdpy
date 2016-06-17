from fwdpy.fwdpy cimport singlepop_t

cdef extern from "fwdpy_fitness.hpp" namespace "fwdpy" nogil:
    cdef cppclass singlepop_fitness:
        singlepop_fitness()
        void update(const singlepop_t *)

    singlepop_fitness make_additive_fitness(double scaling)
    singlepop_fitness make_multiplicative_fitness(double scaling)
    
cdef class singlepopFitness(object):
    """
    Base object for single-deme fitness functions
    """
    cdef singlepop_fitness wfxn
    
cdef class singlepopAdditive(singlepopFitness):
    pass

cdef class singlepopMult(singlepopFitness):
    pass
