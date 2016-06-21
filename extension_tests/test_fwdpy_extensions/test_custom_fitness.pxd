from fwdpy.fwdpp cimport popgenmut
from fwdpy.fitness cimport singlepopFitness

cdef class additiveFitnessTesting(singlepopFitness):
    pass

cdef class Aa_only_testing(singlepopFitness):
    pass

cdef inline void additive_het_testing(double & w, const popgenmut & m):
    (&w)[0]+=m.s*m.h

cdef inline void additive_hom_testing(double & w, const popgenmut & m):
    (&w)[0]+=m.s*2.0
    
cdef inline void Aa_only_het_testing(double & w, const popgenmut & m):
    (&w)[0] += m.s*m.h

cdef inline void Aa_only_hom_testing(double & w, const popgenmut & m):
    return
