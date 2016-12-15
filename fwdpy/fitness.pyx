cdef class SpopAdditive(SpopFitness):
    """
    Additive fitness model for a single deme.
    """
    def __cinit__(self,int scaling = 2):
        """
        Constructor

        :param scaling: For a single mutation, fitness is calculated as 1, 1+sh, and 1+scaling*s for genotypes AA, Aa, and aa, respectively
        """
        self.wfxn = unique_ptr[singlepop_fitness](new singlepop_fitness(het_additive_update,
                                                                        choose_additive_hom_updater(scaling),
                                                                        return_w_plus1,
                                                                        0.0))
        
cdef class SpopMult(SpopFitness):
    """
    Multiplicative fitness model for a single deme.
    """
    def __cinit__(self,int scaling = 1):
        """
        Constructor

        :param scaling: For a single mutation, fitness is calculated as 1, 1+sh, and 1+scaling*s for genotypes AA, Aa, and aa, respectively
        """
        self.wfxn = unique_ptr[singlepop_fitness](new singlepop_fitness(het_mult_update,
                                                                        choose_mult_hom_updater(scaling),
                                                                        return_w,
                                                                        1.0))
        
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
