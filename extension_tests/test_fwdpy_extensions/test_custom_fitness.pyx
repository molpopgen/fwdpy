from fwdpy.fitness cimport singlepopFitness,genotype_fitness_updater,fitness_function_finalizer,make_custom_fitness,return_w_plus1
from fwdpy.fwdpp cimport popgenmut

cdef void additive_het_testing(double & w, const popgenmut & m):
    (&w)[0]+=m.s*m.h

cdef void additive_hom_testing(double & w, const popgenmut & m):
    (&w)[0]+=m.s*2.0

cdef class additiveFitnessTesting(singlepopFitness):
    def __cinit__(self):
        self.wfxn = make_custom_fitness(<genotype_fitness_updater>additive_het_testing,
                                        <genotype_fitness_updater>additive_hom_testing,
                                        <fitness_function_finalizer>return_w_plus1,
                                        0.0)

cdef void Aa_only_het_testing(double & w, const popgenmut & m):
    (&w)[0] += m.s*m.h

cdef void Aa_only_hom_testing(double & w, const popgenmut & m):
    return

cdef class AA_only_testing(singlepopFitness):
    def __cinit__(self):
        self.wfxn = make_custom_fitness(<genotype_fitness_updater>Aa_only_het_testing,
                                        <genotype_fitness_updater>Aa_only_hom_testing,
                                        <fitness_function_finalizer>return_w_plus1,
                                        0.0)
