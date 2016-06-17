cdef class singlepopAdditive(singlepopFitness):
    def __cinit__(self,double scaling = 2.0):
        self.wfxn = make_additive_fitness(scaling)

cdef class singlepopMult(singlepopFitness):
    def __cinit__(self,double scaling = 2.0):
        self.wfxn = make_multiplicative_fitness(scaling)

cdef class singlepopGBR(singlepopFitness):
    def __cinit__(self):
        self.wfxn = make_gbr_fitness()
