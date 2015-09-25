def evolve_qtrait(GSLrng rng,
                    int npops,
                    int N,
                    unsigned[:] nlist,
                    double mu_neutral,
                    double mu_selected,
                    double recrate,
                    list nregions,
                    list sregions,
                    list recregions,
                    double sigmaE,
                    double optimum = 0.,
                    bint track = False,
                    double f = 0.):
    """
    Evolve a quantitative trait with variable mutation, fitness effects, and recombination rates.

    :param rng: a :class:`GSLrng`
    :param npops: The number of populations to simulate.  This is equal to the number of threads that will be used!
    :param N: The diploid population size to simulate
    :param nlist: An array view of a NumPy array.  This represents the population sizes over time.  The length of this view is the length of the simulation in generations. The view must be of an array of 32 bit, unsigned integers (see example).
    :param mu_neutral: The mutation rate to variants not affecting fitness ("neutral" mutations).  The unit is per gamete, per generation.
    :param mu_selected: The mutation rate to variants affecting fitness ("selected" mutations).  The unit is per gamete, per generation.
    :param recrate: The recombination rate in the regions (per diploid, per generation)
    :param nregions: A list specifying where neutral mutations occur
    :param sregions: A list specifying where selected mutations occur
    :param recregions: A list specifying how the genetic map varies along the region
    :param sigmaE: The standard deviation in random variation to add to trait value
    :param optimum: The optimum trait value.
    :param track: whether or not to record the frequency trajectories of mutations.  True = simulations are much slower!
    :param f: The selfing probabilty
    """
    pops = popvec(npops,N)
    nreg = internal.process_regions(nregions)
    sreg = internal.process_regions(sregions)
    recreg = internal.process_regions(recregions)
    v = shwrappervec()
    internal.process_sregion_callbacks(v,sregions)
    evolve_qtraits_t(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,optimum,track,
                    nreg['beg'].tolist(),nreg['end'].tolist(),nreg['weight'].tolist(),
                    sreg['beg'].tolist(),sreg['end'].tolist(),sreg['weight'].tolist(),&v.vec,
                    recreg['beg'].tolist(),recreg['end'].tolist(),recreg['weight'].tolist())
    return pops

def evolve_qtrait_more(GSLrng rng,
                    popvec pops,
                    unsigned[:] nlist,
                    double mu_neutral,
                    double mu_selected,
                    double recrate,
                    list nregions,
                    list sregions,
                    list recregions,
                    double sigmaE,
                    double optimum = 0.,
                    bint track = False,
                    double f = 0.):
    """
    Continue to evolve a quantitative trait with variable mutation, fitness effects, and recombination rates.

    :param rng: a :class:`GSLrng`
    :param pops: A list of populations simulated by :func:`fwdpy.qtrait.qtrait.evolve_qtrait`
    :param N: The diploid population size to simulate
    :param nlist: An array view of a NumPy array.  This represents the population sizes over time.  The length of this view is the length of the simulation in generations. The view must be of an array of 32 bit, unsigned integers (see example).
    :param mu_neutral: The mutation rate to variants not affecting fitness ("neutral" mutations).  The unit is per gamete, per generation.
    :param mu_selected: The mutation rate to variants affecting fitness ("selected" mutations).  The unit is per gamete, per generation.
    :param recrate: The recombination rate in the regions (per diploid, per generation)
    :param nregions: A list specifying where neutral mutations occur
    :param sregions: A list specifying where selected mutations occur
    :param recregions: A list specifying how the genetic map varies along the region
    :param sigmaE: The standard deviation in random variation to add to trait value
    :param optimum: The optimum trait value.
    :param track: whether or not to record the frequency trajectories of mutations.  True = simulations are much slower!
    :param f: The selfing probabilty
    """
    nreg = internal.process_regions(nregions)
    sreg = internal.process_regions(sregions)
    recreg = internal.process_regions(recregions)
    v = shwrappervec()
    internal.process_sregion_callbacks(v,sregions)
    evolve_qtraits_t(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,optimum,track,
                    nreg['beg'].tolist(),nreg['end'].tolist(),nreg['weight'].tolist(),
                    sreg['beg'].tolist(),sreg['end'].tolist(),sreg['weight'].tolist(),&v.vec,
                    recreg['beg'].tolist(),recreg['end'].tolist(),recreg['weight'].tolist())
