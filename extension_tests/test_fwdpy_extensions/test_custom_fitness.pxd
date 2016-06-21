from fwdpy.fwdpp cimport popgenmut,gamete_base
from fwdpy.fitness cimport singlepopFitness
from libcpp.vector cimport vector
from libc.math cimport sqrt
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

ctypedef gamete_base[void] gamete_t
ctypedef vector[popgenmut] mcont_t;

cdef inline double addEsizes(const gamete_t & g, const mcont_t & m):
    cdef size_t i=0,n=g.smutations.size()
    cdef double sum = 0.0
    while i<n:
        sum+=m[g.smutations[i]].s
        i+=1
    return sum

cdef inline double geomean(double e1, double e2):
    return sqrt(e1*e2)

cdef class GBR(singlepopFitness):
    pass
    
