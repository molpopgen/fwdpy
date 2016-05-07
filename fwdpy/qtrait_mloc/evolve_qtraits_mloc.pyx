cimport cython

def evolve_qtraits_mloc_sample(GSLrng rng_evolve,
                               GSLrng rng_sample,
                               popvec_mloc pops,
                               unsigned[:] nlist,
                               const vector[double] & mu_neutral,
                               const vector[double] & mu_selected,
                               const vector[double] & sigmas,
                               const vector[double] & recrates_within,
                               const vector[double] & recrates_between,
                               int sample,
                               int nsam,
                               double optimum = 0.0,
                               double sigmaE = 0.0,
                               double f = 0.0,
                               double VS = 1.0):
    """
    Evolve a set of partially-linked regions contributing to variation
    in a trait under real Gaussian stabilizing selection. 

    This version takes a sample from the population at regular intervals.

    :param rng_evolve: a :class:`fwdpy.fwdpy.GSLrng` used for evolving populations
    :param rng_sample: a :class:`fwdpy.fwdpy.GSLrng` used for sampling from populations
    :param pops: a :class:`fwdpy.fwdpy.popvec_mloc`
    :param nlist: Population sizes over time (array of 32 bit unsigned integers)
    :param mu_neutral: Mutation rates to neutral variants at each locus. Per gamete, per generation
    :param mu_selected: Mutation rates to non-neutral variants at each locus. Per gamete, per generation
    :param sigmas: Effect sizes are N(0,sigma[i]) at the i-th locus.
    :param recrates_within: Recombination rates within each locus.  Per diploid, per generation
    :param recrates_betwen: Genetic distance, in centiMorgans, between each pair of loci.
    :param sample: Take a sample of size nsam every "sample" generations
    :param nsam: Sample size (no. chromosomes)
    :param optimum: Optimum trait value
    :param sigmaE: Std. dev. of Gaussian noise to add to trait values.
    :param f: Selfing probability
    :param VS: Intensity of selection against extreme phenotypes.  Smaller = more selection.

    This is a pure additive effects model.

    When evolving to an equilibrium around an optimum of 0, the total :math:`VG` for the trait would be 
    approximately :math:`4VS\\sum_i\\mu_i`.

    If :math:`P` is a diploid's phenotype, fitness is :math:`w=e^-\\frac{(P-Opt)^2}{2VS}`.
    """
    #Note: exceptions are thrown from the C++ side...
    cdef size_t listlen = len(nlist)
    return evolve_qtrait_mloc_sample_async(rng_evolve.thisptr,
                                           rng_sample.thisptr,
                                           &pops.pops,
                                           &nlist[0],listlen,
                                           mu_neutral,mu_selected,sigmas,
                                           recrates_within,recrates_between,
                                           f,sigmaE,optimum,VS,sample,nsam)


def evolve_qtraits_mloc_track(GSLrng rng_evolve,
                              popvec_mloc pops,
                              unsigned[:] nlist,
                              const vector[double] & mu_neutral,
                              const vector[double] & mu_selected,
                              const vector[double] & sigmas,
                              const vector[double] & recrates_within,
                              const vector[double] & recrates_between,
                              int sample,
                              double optimum = 0.0,
                              double sigmaE = 0.0,
                              double f = 0.0,
                              double VS = 1.0):
    """
    Evolve a set of partially-linked regions contributing to variation
    in a trait under real Gaussian stabilizing selection. 

    This version records the frequencies of all mutations affecting trait/fitness
    at regular intervals.

    :param rng_evolve: a :class:`fwdpy.fwdpy.GSLrng` used for evolving populations
    :param pops: a :class:`fwdpy.fwdpy.popvec_mloc`
    :param nlist: Population sizes over time (array of 32 bit unsigned integers)
    :param mu_neutral: Mutation rates to neutral variants at each locus. Per gamete, per generation
    :param mu_selected: Mutation rates to non-neutral variants at each locus. Per gamete, per generation
    :param sigmas: Effect sizes are N(0,sigma[i]) at the i-th locus.
    :param recrates_within: Recombination rates within each locus.  Per diploid, per generation
    :param recrates_betwen: Genetic distance, in centiMorgans, between each pair of loci.
    :param nsam: Sample size (no. chromosomes)
    :param optimum: Optimum trait value
    :param sigmaE: Std. dev. of Gaussian noise to add to trait values.
    :param f: Selfing probability
    :param VS: Intensity of selection against extreme phenotypes.  Smaller = more selection.

    This is a pure additive effects model.

    When evolving to an equilibrium around an optimum of 0, the total :math:`VG` for the trait would be 
    approximately :math:`4VS\\sum_i\\mu_i`.

    If :math:`P` is a diploid's phenotype, fitness is :math:`w=e^-\\frac{(P-Opt)^2}{2VS}`.
    """
    cdef size_t listlen = len(nlist)
    return evolve_qtrait_mloc_track_async(rng_evolve.thisptr,
                                          &pops.pops,
                                          &nlist[0],listlen,
                                          mu_neutral,mu_selected,sigmas,
                                          recrates_within,recrates_between,
                                          f,sigmaE,optimum,VS,sample)

def evolve_qtraits_mloc_VA(GSLrng rng_evolve,
                           popvec_mloc pops,
                           unsigned[:] nlist,
                           const vector[double] & mu_neutral,
                           const vector[double] & mu_selected,
                           const vector[double] & sigmas,
                           const vector[double] & recrates_within,
                           const vector[double] & recrates_between,
                           int sample,
                           double optimum = 0.0,
                           double sigmaE = 0.0,
                           double f = 0.0,
                           double VS = 1.0):
    """
    Evolve a set of partially-linked regions contributing to variation
    in a trait under real Gaussian stabilizing selection. 

    This version tracks the contribution of mutations to additive genetic variance.

    :param rng_evolve: a :class:`fwdpy.fwdpy.GSLrng` used for evolving populations
    :param pops: a :class:`fwdpy.fwdpy.popvec_mloc`
    :param nlist: Population sizes over time (array of 32 bit unsigned integers)
    :param mu_neutral: Mutation rates to neutral variants at each locus. Per gamete, per generation
    :param mu_selected: Mutation rates to non-neutral variants at each locus. Per gamete, per generation
    :param sigmas: Effect sizes are N(0,sigma[i]) at the i-th locus.
    :param recrates_within: Recombination rates within each locus.  Per diploid, per generation
    :param recrates_betwen: Genetic distance, in centiMorgans, between each pair of loci.
    :param nsam: Sample size (no. chromosomes)
    :param optimum: Optimum trait value
    :param sigmaE: Std. dev. of Gaussian noise to add to trait values.
    :param f: Selfing probability
    :param VS: Intensity of selection against extreme phenotypes.  Smaller = more selection.

    This is a pure additive effects model.

    When evolving to an equilibrium around an optimum of 0, the total :math:`VG` for the trait would be 
    approximately :math:`4VS\\sum_i\\mu_i`.

    If :math:`P` is a diploid's phenotype, fitness is :math:`w=e^-\\frac{(P-Opt)^2}{2VS}`.
    """
    cdef size_t listlen = len(nlist)
    return evolve_qtrait_mloc_VA_async(rng_evolve.thisptr,
                                       &pops.pops,
                                       &nlist[0],listlen,
                                       mu_neutral,mu_selected,sigmas,
                                       recrates_within,recrates_between,
                                       f,sigmaE,optimum,VS,sample)

def evolve_qtraits_mloc_popstats(GSLrng rng_evolve,
                                 popvec_mloc pops,
                                 unsigned[:] nlist,
                                 const vector[double] & mu_neutral,
                                 const vector[double] & mu_selected,
                                 const vector[double] & sigmas,
                                 const vector[double] & recrates_within,
                                 const vector[double] & recrates_between,
                                 int sample,
                                 double optimum = 0.0,
                                 double sigmaE = 0.0,
                                 double f = 0.0,
                                 double VS = 1.0):
    """
    Evolve a set of partially-linked regions contributing to variation
    in a trait under real Gaussian stabilizing selection. 

    This version records various "quantitative gentics" statistics from the population
    at regular intervals.

    :param rng_evolve: a :class:`fwdpy.fwdpy.GSLrng` used for evolving populations
    :param pops: a :class:`fwdpy.fwdpy.popvec_mloc`
    :param nlist: Population sizes over time (array of 32 bit unsigned integers)
    :param mu_neutral: Mutation rates to neutral variants at each locus. Per gamete, per generation
    :param mu_selected: Mutation rates to non-neutral variants at each locus. Per gamete, per generation
    :param sigmas: Effect sizes are N(0,sigma[i]) at the i-th locus.
    :param recrates_within: Recombination rates within each locus.  Per diploid, per generation
    :param recrates_betwen: Genetic distance, in centiMorgans, between each pair of loci.
    :param nsam: Sample size (no. chromosomes)
    :param optimum: Optimum trait value
    :param sigmaE: Std. dev. of Gaussian noise to add to trait values.
    :param f: Selfing probability
    :param VS: Intensity of selection against extreme phenotypes.  Smaller = more selection.

    This is a pure additive effects model.

    When evolving to an equilibrium around an optimum of 0, the total :math:`VG` for the trait would be 
    approximately :math:`4VS\\sum_i\\mu_i`.

    If :math:`P` is a diploid's phenotype, fitness is :math:`w=e^-\\frac{(P-Opt)^2}{2VS}`.
    """
    cdef size_t listlen = len(nlist)
    return evolve_qtrait_mloc_popstats_async(rng_evolve.thisptr,
                                             &pops.pops,
                                             &nlist[0],listlen,
                                             mu_neutral,mu_selected,sigmas,
                                             recrates_within,recrates_between,
                                             f,sigmaE,optimum,VS,sample)

@cython.boundscheck(False)
def evolve_qtraits_mloc(GSLrng rng_evolve,
                        unsigned npops,
                        unsigned N,
                        unsigned nloci,
                        unsigned[:] nlist,
                        const vector[double] & mu_neutral,
                        const vector[double] & mu_selected,
                        const vector[double] & sigmas,
                        const vector[double] & recrates_within,
                        const vector[double] & recrates_between,
                        double optimum = 0.0,
                        double sigmaE = 0.0,
                        double f = 0.0,
                        double VS = 1.0):
    """
    Evolve a set of partially-linked regions contributing to variation
    in a trait under real Gaussian stabilizing selection. 

    This version records various "quantitative gentics" statistics from the population
    at regular intervals.

    :param rng_evolve: a :class:`fwdpy.fwdpy.GSLrng` used for evolving populations
    :param pops: a :class:`fwdpy.fwdpy.popvec_mloc`
    :param nlist: Population sizes over time (array of 32 bit unsigned integers)
    :param mu_neutral: Mutation rates to neutral variants at each locus. Per gamete, per generation
    :param mu_selected: Mutation rates to non-neutral variants at each locus. Per gamete, per generation
    :param sigmas: Effect sizes are N(0,sigma[i]) at the i-th locus.
    :param recrates_within: Recombination rates within each locus.  Per diploid, per generation
    :param recrates_betwen: Genetic distance, in centiMorgans, between each pair of loci.
    :param nsam: Sample size (no. chromosomes)
    :param optimum: Optimum trait value
    :param sigmaE: Std. dev. of Gaussian noise to add to trait values.
    :param f: Selfing probability
    :param VS: Intensity of selection against extreme phenotypes.  Smaller = more selection.

    This is a pure additive effects model.

    When evolving to an equilibrium around an optimum of 0, the total :math:`VG` for the trait would be 
    approximately :math:`4VS\\sum_i\\mu_i`.

    If :math:`P` is a diploid's phenotype, fitness is :math:`w=e^-\\frac{(P-Opt)^2}{2VS}`.
    """
    cdef size_t listlen = len(nlist)
    pops = popvec_mloc(npops,N,nloci)
    with nogil:
        evolve_qtrait_mloc_no_sampling_async(rng_evolve.thisptr,
                                             &pops.pops,
                                             &nlist[0],listlen,
                                             mu_neutral,mu_selected,sigmas,
                                             recrates_within,recrates_between,
                                             f,sigmaE,optimum,VS)
    return pops

#POWER MEAN FUNCTIONS

def evolve_qtraits_mloc_pm_sample(GSLrng rng_evolve,
                                  GSLrng rng_sample,
                                  popvec_mloc pops,
                                  unsigned[:] nlist,
                                  const vector[double] & mu_neutral,
                                  const vector[double] & mu_selected,
                                  const vector[double] & sigmas,
                                  const vector[double] & recrates_within,
                                  const vector[double] & recrates_between,
                                  int sample,
                                  int nsam,
                                  vector[double] & SLd,
                                  double & SLp,
                                  vector[double] & MLd,
                                  double & MLp,
                                  double optimum = 0.0,
                                  double sigmaE = 0.0,
                                  double f = 0.0,
                                  double VS = 1.0):
    """
    Evolve a set of partially-linked regions contributing to variation
    in a trait under real Gaussian stabilizing selection. Let trait values be
    determined by a nested weighted power mean function. 

    This version takes a sample from the population at regular intervals.

    :param rng_evolve: a :class:`fwdpy.fwdpy.GSLrng` used for evolving populations
    :param rng_sample: a :class:`fwdpy.fwdpy.GSLrng` used for sampling from populations
    :param pops: a :class:`fwdpy.fwdpy.popvec_mloc`
    :param nlist: Population sizes over time (array of 32 bit unsigned integers)
    :param mu_neutral: Mutation rates to neutral variants at each locus. Per gamete, per generation
    :param mu_selected: Mutation rates to non-neutral variants at each locus. Per gamete, per generation
    :param sigmas: Effect sizes are N(0,sigma[i]) at the i-th locus.
    :param recrates_within: Recombination rates within each locus.  Per diploid, per generation
    :param recrates_betwen: Genetic distance, in centiMorgans, between each pair of loci.
    :param sample: Take a sample of size nsam every "sample" generations
    :param nsam: Sample size (no. chromosomes)
    :param SLd: Single locus power mean weights
    :param SLp: Single locus power mean power
    :param MLd: Multilocus power mean weights
    :param MLp: Multilocus power mean power
    :param optimum: Optimum trait value
    :param sigmaE: Std. dev. of Gaussian noise to add to trait values.
    :param f: Selfing probability
    :param VS: Intensity of selection against extreme phenotypes.  Smaller = more selection.

    

    When evolving to an equilibrium around an optimum of 0, if all weights and powers were equal to 1, the total :math:`VG` for the trait would be 
    approximately :math:`4VS\\sum_i\\mu_i`.

    If :math:`P` is a diploid's phenotype, fitness is :math:`w=e^-\\frac{(P-Opt)^2}{2VS}`.
    """
    #Note: exceptions are thrown from the C++ side...
    cdef size_t listlen = len(nlist)
    return evolve_qtrait_mloc_pm_sample_async(rng_evolve.thisptr,
                                              rng_sample.thisptr,
                                              &pops.pops,
                                              &nlist[0],listlen,
                                              mu_neutral,mu_selected,sigmas,
                                              recrates_within,recrates_between,
                                              f,sigmaE,optimum,VS,SLd,SLp,MLd,MLp,sample,nsam)


def evolve_qtraits_mloc_pm_track(GSLrng rng_evolve,
                                 popvec_mloc pops,
                                 unsigned[:] nlist,
                                 const vector[double] & mu_neutral,
                                 const vector[double] & mu_selected,
                                 const vector[double] & sigmas,
                                 const vector[double] & recrates_within,
                                 const vector[double] & recrates_between,
                                 int sample,
                                 vector[double] & SLd,
                                 double & SLp,
                                 vector[double] & MLd,
                                 double & MLp,
                                 double optimum = 0.0,
                                 double sigmaE = 0.0,
                                 double f = 0.0,
                                 double VS = 1.0):
    """
    Evolve a set of partially-linked regions contributing to variation
    in a trait under real Gaussian stabilizing selection. Let trait values be
    determined by a nested weighted power mean function. 


    This version records the frequencies of all mutations affecting trait/fitness
    at regular intervals.

    :param rng_evolve: a :class:`fwdpy.fwdpy.GSLrng` used for evolving populations
    :param pops: a :class:`fwdpy.fwdpy.popvec_mloc`
    :param nlist: Population sizes over time (array of 32 bit unsigned integers)
    :param mu_neutral: Mutation rates to neutral variants at each locus. Per gamete, per generation
    :param mu_selected: Mutation rates to non-neutral variants at each locus. Per gamete, per generation
    :param sigmas: Effect sizes are N(0,sigma[i]) at the i-th locus.
    :param recrates_within: Recombination rates within each locus.  Per diploid, per generation
    :param recrates_betwen: Genetic distance, in centiMorgans, between each pair of loci.
    :param nsam: Sample size (no. chromosomes)
    :param SLd: Single locus power mean weights
    :param SLp: Single locus power mean power
    :param MLd: Multilocus power mean weights
    :param MLp: Multilocus power mean power
    :param optimum: Optimum trait value
    :param sigmaE: Std. dev. of Gaussian noise to add to trait values.
    :param f: Selfing probability
    :param VS: Intensity of selection against extreme phenotypes.  Smaller = more selection.


    When evolving to an equilibrium around an optimum of 0, if all weights and powers were equal to 1,the total :math:`VG` for the trait would be 
    approximately :math:`4VS\\sum_i\\mu_i`. 

    If :math:`P` is a diploid's phenotype, fitness is :math:`w=e^-\\frac{(P-Opt)^2}{2VS}`.
    """
    cdef size_t listlen = len(nlist)
    return evolve_qtrait_mloc_pm_track_async(rng_evolve.thisptr,
                                          &pops.pops,
                                          &nlist[0],listlen,
                                          mu_neutral,mu_selected,sigmas,
                                          recrates_within,recrates_between,
                                          f,sigmaE,optimum,VS,SLd,SLp,MLd,MLp,sample)

def evolve_qtraits_mloc_pm_VA(GSLrng rng_evolve,
                           popvec_mloc pops,
                           unsigned[:] nlist,
                           const vector[double] & mu_neutral,
                           const vector[double] & mu_selected,
                           const vector[double] & sigmas,
                           const vector[double] & recrates_within,
                           const vector[double] & recrates_between,
                           int sample,
                           vector[double] & SLd,
                           double & SLp,
                           vector[double] & MLd,
                           double & MLp,
                           double optimum = 0.0,
                           double sigmaE = 0.0,
                           double f = 0.0,
                           double VS = 1.0):
    """
    Evolve a set of partially-linked regions contributing to variation
    in a trait under real Gaussian stabilizing selection. Let trait values
    be determined by a nested weighted power mean function.

    This version tracks the contribution of mutations to additive genetic variance.

    :param rng_evolve: a :class:`fwdpy.fwdpy.GSLrng` used for evolving populations
    :param pops: a :class:`fwdpy.fwdpy.popvec_mloc`
    :param nlist: Population sizes over time (array of 32 bit unsigned integers)
    :param mu_neutral: Mutation rates to neutral variants at each locus. Per gamete, per generation
    :param mu_selected: Mutation rates to non-neutral variants at each locus. Per gamete, per generation
    :param sigmas: Effect sizes are N(0,sigma[i]) at the i-th locus.
    :param recrates_within: Recombination rates within each locus.  Per diploid, per generation
    :param recrates_betwen: Genetic distance, in centiMorgans, between each pair of loci.
    :param nsam: Sample size (no. chromosomes)
    :param SLd: Single locus power mean weights
    :param SLp: Single locus power mean power
    :param MLd: Multilocus power mean weights
    :param MLp: Multilocus power mean power
    :param optimum: Optimum trait value
    :param sigmaE: Std. dev. of Gaussian noise to add to trait values.
    :param f: Selfing probability
    :param VS: Intensity of selection against extreme phenotypes.  Smaller = more selection.


    When evolving to an equilibrium around an optimum of 0,if all weights and powers were equal to 1, the total :math:`VG` for the trait would be 
    approximately :math:`4VS\\sum_i\\mu_i`.

    If :math:`P` is a diploid's phenotype, fitness is :math:`w=e^-\\frac{(P-Opt)^2}{2VS}`.
    """
    cdef size_t listlen = len(nlist)
    return evolve_qtrait_mloc_pm_VA_async(rng_evolve.thisptr,
                                       &pops.pops,
                                       &nlist[0],listlen,
                                       mu_neutral,mu_selected,sigmas,
                                       recrates_within,recrates_between,
                                       f,sigmaE,optimum,VS,SLd,SLp,MLd,MLp,sample)

def evolve_qtraits_mloc_pm_popstats(GSLrng rng_evolve,
                                 popvec_mloc pops,
                                 unsigned[:] nlist,
                                 const vector[double] & mu_neutral,
                                 const vector[double] & mu_selected,
                                 const vector[double] & sigmas,
                                 const vector[double] & recrates_within,
                                 const vector[double] & recrates_between,
                                 int sample,
                                 vector[double] & SLd,
                                 double & SLp,
                                 vector[double] & MLd,
                                 double & MLp,
                                 double optimum = 0.0,
                                 double sigmaE = 0.0,
                                 double f = 0.0,
                                 double VS = 1.0):
    """
    Evolve a set of partially-linked regions contributing to variation
    in a trait under real Gaussian stabilizing selection. Let trait values
    be determined by a nested weighted power mean function.

    This version records various "quantitative gentics" statistics from the population
    at regular intervals.

    :param rng_evolve: a :class:`fwdpy.fwdpy.GSLrng` used for evolving populations
    :param pops: a :class:`fwdpy.fwdpy.popvec_mloc`
    :param nlist: Population sizes over time (array of 32 bit unsigned integers)
    :param mu_neutral: Mutation rates to neutral variants at each locus. Per gamete, per generation
    :param mu_selected: Mutation rates to non-neutral variants at each locus. Per gamete, per generation
    :param sigmas: Effect sizes are N(0,sigma[i]) at the i-th locus.
    :param recrates_within: Recombination rates within each locus.  Per diploid, per generation
    :param recrates_betwen: Genetic distance, in centiMorgans, between each pair of loci.
    :param nsam: Sample size (no. chromosomes)
    :param SLd: Single locus power mean weights
    :param SLp: Single locus power mean power
    :param MLd: Multilocus power mean weights
    :param MLp: Multilocus power mean power
    :param optimum: Optimum trait value
    :param sigmaE: Std. dev. of Gaussian noise to add to trait values.
    :param f: Selfing probability
    :param VS: Intensity of selection against extreme phenotypes.  Smaller = more selection.

    When evolving to an equilibrium around an optimum of 0, if all weights and powers were equal to 1, the total :math:`VG` for the trait would be 
    approximately :math:`4VS\\sum_i\\mu_i`.

    If :math:`P` is a diploid's phenotype, fitness is :math:`w=e^-\\frac{(P-Opt)^2}{2VS}`.
    """
    cdef size_t listlen = len(nlist)
    return evolve_qtrait_mloc_pm_popstats_async(rng_evolve.thisptr,
                                             &pops.pops,
                                             &nlist[0],listlen,
                                             mu_neutral,mu_selected,sigmas,
                                             recrates_within,recrates_between,
                                             f,sigmaE,optimum,VS,SLd,SLp,MLd,MLp,sample)

@cython.boundscheck(False)
def evolve_qtraits_mloc_pm(GSLrng rng_evolve,
                        unsigned npops,
                        unsigned N,
                        unsigned nloci,
                        unsigned[:] nlist,
                        const vector[double] & mu_neutral,
                        const vector[double] & mu_selected,
                        const vector[double] & sigmas,
                        const vector[double] & recrates_within,
                        const vector[double] & recrates_between,
                        vector[double] & SLd,
                        double & SLp,
                        vector[double] & MLd,
                        double & MLp,
                        double optimum = 0.0,
                        double sigmaE = 0.0,
                        double f = 0.0,
                        double VS = 1.0):
    """
    Evolve a set of partially-linked regions contributing to variation
    in a trait under real Gaussian stabilizing selection. Let trait values
    be determined by a nested weighted power mean function.

    This version records various "quantitative gentics" statistics from the population
    at regular intervals.

    :param rng_evolve: a :class:`fwdpy.fwdpy.GSLrng` used for evolving populations
    :param pops: a :class:`fwdpy.fwdpy.popvec_mloc`
    :param nlist: Population sizes over time (array of 32 bit unsigned integers)
    :param mu_neutral: Mutation rates to neutral variants at each locus. Per gamete, per generation
    :param mu_selected: Mutation rates to non-neutral variants at each locus. Per gamete, per generation
    :param sigmas: Effect sizes are N(0,sigma[i]) at the i-th locus.
    :param recrates_within: Recombination rates within each locus.  Per diploid, per generation
    :param recrates_betwen: Genetic distance, in centiMorgans, between each pair of loci.
    :param nsam: Sample size (no. chromosomes)
    :param SLd: Single locus power mean weights
    :param SLp: Single locus power mean power
    :param MLd: Multilocus power mean weights
    :param MLp: Multilocus power mean power
    :param optimum: Optimum trait value
    :param sigmaE: Std. dev. of Gaussian noise to add to trait values.
    :param f: Selfing probability
    :param VS: Intensity of selection against extreme phenotypes.  Smaller = more selection.


    When evolving to an equilibrium around an optimum of 0, if all weights and powers were equal to 1,the total :math:`VG` for the trait would be 
    approximately :math:`4VS\\sum_i\\mu_i`.

    If :math:`P` is a diploid's phenotype, fitness is :math:`w=e^-\\frac{(P-Opt)^2}{2VS}`.
    """
    cdef size_t listlen = len(nlist)
    pops = popvec_mloc(npops,N,nloci)
    with nogil:
        evolve_qtrait_mloc_pm_no_sampling_async(rng_evolve.thisptr,
                                             &pops.pops,
                                             &nlist[0],listlen,
                                             mu_neutral,mu_selected,sigmas,
                                             recrates_within,recrates_between,
                                             f,sigmaE,optimum,VS,SLd,SLp,MLd,MLp)
    return pops
