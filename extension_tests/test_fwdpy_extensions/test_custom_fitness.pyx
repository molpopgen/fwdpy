#Import the types needed from fwdpp
from fwdpy.fwdpp cimport popgenmut,gamete_base
#Bring all of fwdpy.fitness into scope for convenience:
from fwdpy.fitness cimport *
#We'll need C's sqrt function for the 'GBR' model:
from libc.math cimport sqrt


cdef inline void additive_het_testing(double & w, const popgenmut & m):
    (&w)[0]+=m.s*m.h

cdef inline void additive_hom_testing(double & w, const popgenmut & m):
    (&w)[0]+=m.s*2.0
    
cdef inline void Aa_only_het_testing(double & w, const popgenmut & m):
    (&w)[0] += m.s*m.h

cdef inline void Aa_only_hom_testing(double & w, const popgenmut & m):
    return

cdef inline double addEsizes(const gamete_t & g, const mcont_t & m):
    cdef size_t i=0,n=g.smutations.size()
    cdef double sum = 0.0
    while i<n:
        sum+=m[g.smutations[i]].s
        i+=1
    return sum

cdef inline double geomean(double e1, double e2):
    return sqrt(e1*e2)

cdef class AdditiveFitnessTesting(SpopFitness):
    def __cinit__(self):
        self.wfxn = singlepop_fitness(<genotype_fitness_updater>additive_het_testing,
                                      <genotype_fitness_updater>additive_hom_testing,
                                      <fitness_function_finalizer>return_w_plus1,
                                      0.0)

cdef class AaOnlyTesting(SpopFitness):
    def __cinit__(self):
        self.wfxn = singlepop_fitness(<genotype_fitness_updater>Aa_only_het_testing,
                                      <genotype_fitness_updater>Aa_only_hom_testing,
                                      <fitness_function_finalizer>return_w_plus1,
                                      0.0)

cdef class GBRFitness(SpopFitness):
    def __cinit__(self):
        self.wfxn = singlepop_fitness(<haplotype_fitness_fxn>addEsizes,
                                      <haplotype_fitness_fxn_finalizer>geomean)
