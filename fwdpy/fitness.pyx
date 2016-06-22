cdef class SinglepopAdditive(SinglepopFitness):
    """
    Additive fitness model for a single deme.
    """
    def __cinit__(self,double scaling = 2.0):
        """
        Constructor

        :param scaling: For a single mutation, fitness is calculated as 1, 1+sh, and 1+scaling*s for genotypes AA, Aa, and aa, respectively
        """
        self.wfxn = make_additive_fitness(scaling)

cdef class SinglepopMult(SinglepopFitness):
    """
    Multiplicative fitness model for a single deme.
    """
    def __cinit__(self,double scaling = 2.0):
        """
        Constructor

        :param scaling: For a single mutation, fitness is calculated as 1, 1+sh, and 1+scaling*s for genotypes AA, Aa, and aa, respectively
        """
        self.wfxn = make_multiplicative_fitness(scaling)

cdef class SinglepopGBR(SinglepopFitness):
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

cdef class MultilocusAdditive(MultilocusFitness):
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

cdef class MultilocusGBR(MultilocusFitness):
    """
    The "gene-based recessive" (GBR) model of Thornton et al (2013) PLoS Genetics.

    The GBR model is applied to each region, and then the final value is the sum over regions.

    The fitness (or genetic value as it was used in that paper) is the geometric mean 
    of haplotype effect sizes.

    .. note:: Be really careful with this one!  Fitnesses are undefined if the sum
    of effect sizes on a haplotype is :math:`< 0:`.  The intended use case is to calculate
    a trait value under models with effect sizes :math:`>0`.
    """ 
    def __cinit__(self):
        self.wfxn=make_mloc_gbr_fitness()

cdef class MultilocusMult(MultilocusFitness):
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
