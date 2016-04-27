def evolve_qtraits_mloc_sample(GSLrng rng_evolve,
                               GSLrng rng_sample,
                               popvec_mloc pops,
                               unsigned[:] nlist,
                               list mu_neutral,
                               list mu_selected,
                               list sigmas,
                               list recrates_within,
                               list recrates_between,
                               double f,
                               double sigmaE,
                               double optimum,
                               double VS,
                               int sample,
                               int nsam):
    """
    Evolve a set of partially-linked regions contributing to variation
    in a trait under real stabilizing selection.
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

