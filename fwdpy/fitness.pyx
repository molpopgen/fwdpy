cdef class SpopAdditive(SpopFitness):
    """
    Additive fitness model for a single deme.
    """
    def __cinit__(self,double scaling = 2.0):
        """
        Constructor

        :param scaling: For a single mutation, fitness is calculated as 1, 1+sh, and 1+scaling*s for genotypes AA, Aa, and aa, respectively
        """
        self.wfxn = make_additive_fitness(scaling)

cdef class SpopMult(SpopFitness):
    """
    Multiplicative fitness model for a single deme.
    """
    def __cinit__(self,double scaling = 2.0):
        """
        Constructor

        :param scaling: For a single mutation, fitness is calculated as 1, 1+sh, and 1+scaling*s for genotypes AA, Aa, and aa, respectively
        """
        self.wfxn = make_multiplicative_fitness(scaling)
cdef class MlocusAdditive(MlocusFitness):
    """
    Additive fitness for a multi-region model.

    Fitness is additive with dominance within regions,
    then strictly additive across regions.
    """
    def __cinit__(self,double scaling = 2.0):
        """
        Constructor

        :param scaling: For a single mutation, fitness is calculated as 1, 1+sh, and 1+scaling*s for genotypes AA, Aa, and aa, respectively
        """
        self.wfxn=make_mloc_additive_fitness(scaling)

cdef class MlocusMult(MlocusFitness):
    """
    Multiplicative fitness for a multi-region model.

    Fitness is multiplicative with dominance within regions,
    then strictly additive across regions.
    """
    def __cinit__(self,double scaling = 2.0):
        """
        Constructor

        :param scaling: For a single mutation, fitness is calculated as 1, 1+sh, and 1+scaling*s for genotypes AA, Aa, and aa, respectively
        """
        self.wfxn=make_mloc_multiplicative_fitness(scaling)

