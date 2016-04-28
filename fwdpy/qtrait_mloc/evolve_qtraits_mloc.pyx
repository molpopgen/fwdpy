def evolve_qtraits_mloc_sample(GSLrng rng_evolve,
                               GSLrng rng_sample,
                               popvec_mloc pops,
                               unsigned[:] nlist,
                               list mu_neutral,
                               list mu_selected,
                               list sigmas,
                               list recrates_within,
                               list recrates_between,
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

