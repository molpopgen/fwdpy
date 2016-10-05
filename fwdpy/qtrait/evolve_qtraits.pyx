import warnings,fwdpy
from cython.view cimport array as cvarray
from cpython cimport array
from cython.operator cimport dereference as deref
cimport cython
from fwdpy.fitness cimport SpopFitness
from fwdpy.fitness cimport SpopAdditive

def check_input_params(double sigmaE, double VS):
    if sigmaE < 0.:
        raise RuntimeError("sigmaE must be >= 0.")
    if VS < 0.:
        raise RuntimeError("VS must be >= 0.")

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
def evolve_regions_qtrait(GSLrng rng,
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
    :param optimum: The optimum trait value. **Default = 0.0**
    :param f: The selfing probabilty. **Default = 0.0**
    :param VS: The total variance in selection intensity. **Default = 1.0**

    :raises: RuntimeError if parameters do not pass checks
    """
    pops = SpopVec(npops,N)
    donothing = NothingSampler(npops)
    fitness = SpopAdditiveTrait()
    evolve_regions_qtrait_sampler_fitness(rng,pops,donothing,fitness,nlist,
                                          mu_neutral,mu_selected,recrate,
                                          nregions,sregions,recregions,
                                          len(nlist),sigmaE,optimum,f,VS)
                                          
    return pops

def evolve_regions_qtrait_more(GSLrng rng,
                               SpopVec pops,
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
    donothing = NothingSampler(len(pops))
    fitness = SpopAdditiveTrait()
    evolve_regions_qtrait_sampler_fitness(rng,pops,donothing,fitness,nlist,
                                          mu_neutral,mu_selected,recrate,
                                          nregions,sregions,recregions,
                                          len(nlist),sigmaE,optimum,f,VS)

@cython.boundscheck(False)
def evolve_regions_qtrait_sampler(GSLrng rng,
                                  SpopVec pops,
                                  TemporalSampler slist,
                                  unsigned[:] nlist,
                                  double mu_neutral,
                                  double mu_selected,
                                  double recrate,
                                  list nregions,
                                  list sregions,
                                  list recregions,
                                  int sample,
                                  double sigmaE,
                                  double optimum = 0.0,
                                  double f = 0,
                                  double VS = 1.0):
    fitness = SpopAdditiveTrait()
    evolve_regions_qtrait_sampler_fitness(rng,pops,slist,fitness,nlist,
                                          mu_neutral,mu_selected,recrate,
                                          nregions,sregions,recregions,
                                          sample,sigmaE,optimum,f,VS)
    
@cython.boundscheck(False)
def evolve_regions_qtrait_sampler_fitness(GSLrng rng,
                                          SpopVec pops,
                                          TemporalSampler slist,
                                          SpopFitness fitness_function,
                                          unsigned[:] nlist,
                                          double mu_neutral,
                                          double mu_selected,
                                          double recrate,
                                          list nregions,
                                          list sregions,
                                          list recregions,
                                          int sample,
                                          double sigmaE,
                                          double optimum = 0.0,
                                          double f = 0,
                                          double VS = 1.0):
    fwdpy.check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    if isinstance(fitness_function,SpopGBRTrait):
        check_gbr_sdist(sregions)
    if sample < 0:
        raise RuntimeError("sample must be >= 0")
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    cdef size_t listlen = len(nlist)
    evolve_regions_qtrait_cpp(rng.thisptr,pops.pops,
                              slist.vec,&nlist[0],listlen,mu_neutral,mu_selected,recrate,f,sigmaE,optimum,VS,sample,rmgr.thisptr,deref(fitness_function.wfxn.get()))
