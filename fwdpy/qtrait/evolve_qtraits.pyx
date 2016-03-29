import warnings,fwdpy
from cython.view cimport array as cvarray
from cpython cimport array
cimport cython

def check_input_params(double sigmaE, double VS):
    if sigmaE < 0.:
        raise RuntimeError("sigmaE must be >= 0.")
    if VS < 0.:
        raise RuntimeError("VS must be >= 0.")

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
    fwdpy.check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    check_input_params(sigmaE,VS)
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0

    pops = popvec(npops,N)
    rmgr = region_manager_wrapper();
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    cdef size_t listlen = len(nlist)
    with nogil:
        evolve_qtrait_no_sampling_async(rng.thisptr,&pops.pops,&nlist[0],listlen,mu_neutral,mu_selected,recrate,f,sigmaE,optimum,VS,
                                        rmgr.thisptr)
    return pops

@cython.boundscheck(False)
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
    :param f: The selfing probabilty
    :param VS: The variance in the Gaussian fitness function.  Under certaing strong assumtions, :math:`V(G) \approx 4\times\mu\timesV(S)`, where :math:`\mu` is mu_selected.

    :raises: RuntimeError if parameters do not pass checks
    """
    fwdpy.check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    check_input_params(sigmaE,VS)
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
        rmgr = region_manager_wrapper()
        internal.make_region_manager(rmgr,nregions,sregions,recregions)
    cdef size_t listlen = len(nlist)
    with nogil:
        evolve_qtrait_no_sampling_async(rng.thisptr,&pops.pops,&nlist[0],listlen,mu_neutral,mu_selected,recrate,f,sigmaE,optimum,VS,
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
    fwdpy.check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    check_input_params(sigmaE,VS)
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
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
    fwdpy.check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    check_input_params(sigmaE,VS)
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    if trackStats < 0:
        raise RuntimeError("trackSamples must be >= 0.")

    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    return evolve_qtrait_popstats_async(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,optimum,VS,trackStats,
                                        rmgr.thisptr)

def evolve_qtrait_track(GSLrng rng,
                        popvec pops,
                        unsigned[:] nlist,
                        double mu_neutral,
                        double mu_selected,
                        double recrate,
                        list nregions,
                        list sregions,
                        list recregions,
                        double sigmaE,
                        int track,
                        double optimum = 0.,
                        double f = 0.,
                        double VS = 1):
    fwdpy.check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    check_input_params(sigmaE,VS)
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    if track < 0:
        raise RuntimeError("trackSamples must be >= 0.")

    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    return evolve_qtrait_track_async(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,optimum,VS,track,
                                     rmgr.thisptr)

#Below are functions related to the 'gene-based' recessive models of doi:10.1371/journal.pgen.1003258

def check_gbr_sdist(sregions):
    for i in sregions:
        if isinstance(i,fwdpy.GaussianS):
            raise RuntimeError("Gaussian effects not allowed for this model")
        elif (isinstance(i,fwdpy.GammaS) and i.mean < 0) or (isinstance(i,fwdpy.ExpS) and i.mean<0):
            raise RuntimeError("mean effect size must be >= 0")
        elif isinstance(i,fwdpy.ConstantS) and i.s < 0:
            raise RuntimeError("effect size must be >= 0")
        elif isinstance(i,fwdpy.UniformS) and (i.lo<0 or i.hi<0):
            raise RuntimeError("min and max effect size must be >= 0")

        
@cython.boundscheck(False)
def evolve_gbr(GSLrng rng,
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
    fwdpy.check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    check_gbr_sdist(sregions)
    check_input_params(sigmaE,VS)
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0

    pops = popvec(npops,N)
    rmgr = region_manager_wrapper();
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    cdef size_t listlen = len(nlist)
    with nogil:
        evolve_gbr_no_sampling_async(rng.thisptr,&pops.pops,&nlist[0],listlen,mu_neutral,mu_selected,recrate,f,sigmaE,optimum,VS,
                                     rmgr.thisptr)
    return pops

@cython.boundscheck(False)
def evolve_gbr_sample(GSLrng rng,
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
    fwdpy.check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    check_gbr_sdist(sregions)
    check_input_params(sigmaE,VS)
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    if trackSamples < 0:
        raise RuntimeError("trackSamples must be >= 0.")
    if nsam == 0:
        raise RuntimeError("Sample size (nsam) must be > 0")
    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    
    return evolve_qtrait_sample_async(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,optimum,VS,trackSamples,nsam,
                                      rmgr.thisptr)

@cython.boundscheck(False)
def evolve_gbr_popstats(GSLrng rng,
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
    fwdpy.check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    check_gbr_sdist(sregions)
    check_input_params(sigmaE,VS)
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    if trackStats < 0:
        raise RuntimeError("trackSamples must be >= 0.")

    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    return evolve_gbr_popstats_async(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,optimum,VS,trackStats,
                                     rmgr.thisptr)

@cython.boundscheck(False)
def evolve_gbr_track(GSLrng rng,
                     popvec pops,
                     unsigned[:] nlist,
                     double mu_neutral,
                     double mu_selected,
                     double recrate,
                     list nregions,
                     list sregions,
                     list recregions,
                     double sigmaE,
                     int track,
                     double optimum = 0.,
                     double f = 0.,
                     double VS = 1):
    fwdpy.check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    check_gbr_sdist(sregions)
    check_input_params(sigmaE,VS)
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    if track < 0:
        raise RuntimeError("trackSamples must be >= 0.")

    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    return evolve_gbr_track_async(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,optimum,VS,track,
                                  rmgr.thisptr)
