cdef class singlepopAdditive(singlepopFitness):
    """
    Additive fitness model for a single deme.
    """
    def __cinit__(self,double scaling = 2.0):
        """
        Constructor

        :param scaling: For a single mutation, fitness is calculated as 1, 1+sh, and 1+scaling*s for genotypes AA, Aa, and aa, respectively
        """
        self.wfxn = make_additive_fitness(scaling)

cdef class singlepopMult(singlepopFitness):
    """
    Multiplicative fitness model for a single deme.
    """
    def __cinit__(self,double scaling = 2.0):
        """
        Constructor

        :param scaling: For a single mutation, fitness is calculated as 1, 1+sh, and 1+scaling*s for genotypes AA, Aa, and aa, respectively
        """
        self.wfxn = make_multiplicative_fitness(scaling)

cdef class singlepopGBR(singlepopFitness):
    """
    The "gene-based recessive" (GBR) model of Thornton et al (2013) PLoS Genetics.

    The fitness (or genetic value as it was used in that paper) is the geometric mean 
    of haplotype effect sizes.

    .. note:: Be really careful with this one!  Fitnesses are undefined if the sum
    of effect sizes on a haplotype is :math:`< 0:`.  The intended use case is to calculate
    a trait value under models with effect sizes :math:`>0`.
    """ 
    def __cinit__(self):
        self.wfxn = make_gbr_fitness()
