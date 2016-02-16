import warnings
from cython.view cimport array as cvarray
from cpython cimport array
cimport cython

@cython.boundscheck(False)
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
                  double f = 0.,
                  double VS=1,):
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
    :param f: The selfing probabilty
    :param VS: The total variance in selection intensity

    :raises: RuntimeError if parameters do not pass checks
    """
    if mu_neutral < 0:
        raise RuntimeError("mutation rate to neutral variants must be >= 0.")
    if mu_selected < 0:
        raise RuntimeError("mutation rate to selected variants must be >= 0.")
    if recrate < 0:
        raise RuntimeError("recombination rate must be >= 0.")
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    if sigmaE < 0.:
        raise RuntimeError("sigmaE must be >= 0.")
    if VS < 0.:
        raise RuntimeError("VS must be >= 0.")

    pops = popvec(npops,N)
    rmgr = region_manager_wrapper();
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    cdef unsigned listlen = len(nlist)
    with nogil:
        evolve_qtraits_t(rng.thisptr,&pops.pops,&nlist[0],listlen,mu_neutral,mu_selected,recrate,f,sigmaE,optimum,VS,
                         rmgr.thisptr)
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
                    double f = 0.,
                    double VS = 1,):
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
    :param track: whether or not to record the frequency trajectories of mutations.  If value is x > 0, values are recorded every x generations.  Values < 0 result in a RuntimeError being raised.
    :param trackStats: whether or not to track VG, etc.  If value is x > 0, stats are recorded every x generations.  Values < 0 result in a RuntimeError being raised.
    :param f: The selfing probabilty
    :param VS: The variance in the Gaussian fitness function.  Under certaing strong assumtions, :math:`V(G) \approx 4\times\mu\timesV(S)`, where :math:`\mu` is mu_selected.

    :raises: RuntimeError if parameters do not pass checks
    """
    if mu_neutral < 0:
        raise RuntimeError("mutation rate to neutral variants must be >= 0.")
    if mu_selected < 0:
        raise RuntimeError("mutation rate to selected variants must be >= 0.")
    if recrate < 0:
        raise RuntimeError("recombination rate must be >= 0.")
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    if sigmaE < 0.:
        raise RuntimeError("sigmaE must be >= 0.")
    if VS < 0.:
        raise RuntimeError("VS must be >= 0.")
    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    evolve_qtraits_t(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,optimum,VS,
                     rmgr.thisptr)

def evolve_qtrait_sample(GSLrng rng,
                         popvec pops,
                         unsigned[:] nlist,
                         double mu_neutral,
                         double mu_selected,
                         double recrate,
                         list nregions,
                         list sregions,
                         list recregions,
                         double sigmaE,
                         int trackSamples,
                         unsigned nsam,
                         double optimum = 0.,
                         double f = 0.,
                         double VS = 1):
    if mu_neutral < 0:
        raise RuntimeError("mutation rate to neutral variants must be >= 0.")
    if mu_selected < 0:
        raise RuntimeError("mutation rate to selected variants must be >= 0.")
    if recrate < 0:
        raise RuntimeError("recombination rate must be >= 0.")
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    if sigmaE < 0.:
        raise RuntimeError("sigmaE must be >= 0.")
    if VS < 0.:
        raise RuntimeError("VS must be >= 0.")
    if trackSamples < 0:
        raise RuntimeError("trackSamples must be >= 0.")
    if nsam == 0:
        raise RuntimeError("Sample size (nsam) must be > 0")
    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    return evolve_qtrait_sample_async(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,optimum,VS,trackSamples,nsam,
                                      rmgr.thisptr)

def evolve_qtrait_popstats(GSLrng rng,
                           popvec pops,
                           unsigned[:] nlist,
                           double mu_neutral,
                           double mu_selected,
                           double recrate,
                           list nregions,
                           list sregions,
                           list recregions,
                           double sigmaE,
                           int trackStats,
                           double optimum = 0.,
                           double f = 0.,
                           double VS = 1):
    if mu_neutral < 0:
        raise RuntimeError("mutation rate to neutral variants must be >= 0.")
    if mu_selected < 0:
        raise RuntimeError("mutation rate to selected variants must be >= 0.")
    if recrate < 0:
        raise RuntimeError("recombination rate must be >= 0.")
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    if sigmaE < 0.:
        raise RuntimeError("sigmaE must be >= 0.")
    if VS < 0.:
        raise RuntimeError("VS must be >= 0.")
    if trackStats < 0:
        raise RuntimeError("trackSamples must be >= 0.")

    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    return evolve_qtrait_popstats_async(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,optimum,VS,trackStats,
                                      rmgr.thisptr)
