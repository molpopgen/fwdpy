from fwdpy.fitness cimport genotype_fitness_updater,fitness_function_finalizer,make_custom_fitness,return_w_plus1

cdef class additiveFitnessTesting(singlepopFitness):
    def __cinit__(self):
        self.wfxn = make_custom_fitness(<genotype_fitness_updater>additive_het_testing,
                                        <genotype_fitness_updater>additive_hom_testing,
                                        <fitness_function_finalizer>return_w_plus1,
                                        0.0)

cdef class Aa_only_testing(singlepopFitness):
    def __cinit__(self):
        self.wfxn = make_custom_fitness(<genotype_fitness_updater>Aa_only_het_testing,
                                        <genotype_fitness_updater>Aa_only_hom_testing,
                                        <fitness_function_finalizer>return_w_plus1,
                                        0.0)
