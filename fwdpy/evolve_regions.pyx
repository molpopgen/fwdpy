#See http://docs.cython.org/src/userguide/memoryviews.html
from cython.view cimport array as cvarray
from cpython cimport array
import warnings,math
cimport cython
import internal

def check_input_params(double mu_neutral, double mu_selected, double recrate,
                       list nregions, list sregions, list recregions) :
    for i in [mu_neutral,mu_selected,recrate]:
        if math.isnan(i) or math.isinf(i):
            raise RuntimeError("non-finite value encountered for mutation and/or recombination rates")
    if mu_neutral < 0:
        raise RuntimeError("mutation rate to neutral variants must be >= 0.")
    if mu_selected < 0:
        raise RuntimeError("mutation rate to selected variants must be >= 0.")
    if recrate < 0:
        raise RuntimeError("recombination rate must be >= 0.")
    if mu_neutral > 0 and len(nregions)==0:
        raise RuntimeError("Neutral regions must be defined when mutation rate > 0")
    if mu_selected > 0 and len(sregions)==0:
        raise RuntimeError("Selected regions must be defined when mutation rate > 0")
    if recrate > 0 and len(recregions)==0:
        raise RuntimeError("Recombination regions must be defined when recombination rate > 0")
                       
@cython.boundscheck(False)
def evolve_regions(GSLrng rng,
                    int npops,
                    int N,
                    unsigned[:] nlist,
                    double mu_neutral,
                    double mu_selected,
                    double recrate,
                    list nregions,
                    list sregions,
                    list recregions,
                    double f = 0,
                    const char * fitness = "multiplicative"):
    """
    Evolve a region with variable mutation, fitness effects, and recombination rates.

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
    :param f: The selfing probabilty
    :param fitness: The fitness model.  Must be either "multiplicative" or "additive".

    :raises: RuntimeError if parameters do not pass checks

    Example:

    >>> import fwdpy
    >>> import numpy as np
    >>> #Our "neutral" regions will be from positions [0,1) and [2,3).
    >>> #The regions will have equal weights, and thus will each get
    >>> #1/2 of newly-arising, neutral mutations
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> #Our "selected" mutations will have positions on the continuous
    >>> #interval [1,2).  There will be two classes of such mutations,
    >>> #each with exponentially-distrubted selection coefficients.
    >>> #The first class will have a mean of s = -0.1 (deleterious), and the second
    >>> #will have a mean of 1e-3 (adaptive).  The former will be 100x more common than
    >>> #the latter, as the weights are 1 and 0.01, respectively.
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> #Recombination will be uniform along the interval [0,3).
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> #popsizes are NumPy arrays of 32bit unsigned integers
    >>> #Initial pop size will be N = 1,000, and
    >>> #the type must be uint32 (32 bit unsigned integer)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> #The simulation will be for 10*N generations,
    >>> #so we replicate that value in the array
    >>> popsizes=np.tile(popsizes,100)
    >>> #Simulate 1 deme under this model.
    >>> #The total neutral mutation rate is 1e-3,
    >>> #which is also the recombination rate.
    >>> #The total mutation rate to selected variants is 0.1*(the neutral mutation rate).
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    """
    check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    pops = popvec(npops,N)
    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    cdef size_t listlen = len(nlist)
    with nogil:
        evolve_regions_no_sampling_async(rng.thisptr,&pops.pops,&nlist[0],listlen,mu_neutral,mu_selected,recrate,f,rmgr.thisptr,fitness)
    return pops

@cython.boundscheck(False)
def evolve_regions_more(GSLrng rng,
                        popvec pops,
                        unsigned[:] nlist,
                        double mu_neutral,
                        double mu_selected,
                        double recrate,
                        list nregions,
                        list sregions,
                        list recregions,
                        double f = 0,
                        const char * fitness = "multiplicative"):
    """
    Continue to evolve a region with variable mutation, fitness effects, and recombination rates.

    :param rng: a :class:`GSLrng`
    :param pops: A list of populations simulated by :func:`evolve_regions`
    :param N: The diploid population size to simulate
    :param nlist: An array view of a NumPy array.  This represents the population sizes over time.  The length of this view is the length of the simulation in generations. The view must be of an array of 32 bit, unsigned integers (see example).
    :param mu_neutral: The mutation rate to variants not affecting fitness ("neutral" mutations).  The unit is per gamete, per generation.
    :param mu_selected: The mutation rate to variants affecting fitness ("selected" mutations).  The unit is per gamete, per generation.
    :param recrate: The recombination rate in the regions (per diploid, per generation)
    :param nregions: A list specifying where neutral mutations occur
    :param sregions: A list specifying where selected mutations occur
    :param recregions: A list specifying how the genetic map varies along the region
    :param f: The selfing probabilty
    :param fitness: The fitness model.  Must be either "multiplicative" or "additive".

    :raises: RuntimeError if parameters do not pass checks

    Example:

    >>> # See docstring for fwdpy.evolve_regions for the gory details
    >>> import fwdpy
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> # Evolve for 5N generations initially
    >>> popsizes=np.tile(popsizes,5000)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> # Evolve for another 5N generations
    >>> fwdpy.evolve_regions_more(rng,pops,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    """
    check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    cdef size_t listlen = len(nlist)
    with nogil:
        evolve_regions_no_sampling_async(rng.thisptr,&pops.pops,&nlist[0],listlen,mu_neutral,mu_selected,recrate,f,rmgr.thisptr,fitness)

@cython.boundscheck(False)
def evolve_regions_sample(GSLrng rng,
                            popvec pops,
                            unsigned[:] nlist,
                            double mu_neutral,
                            double mu_selected,
                            double recrate,
                            list nregions,
                            list sregions,
                            list recregions,
                            int sample,
                            unsigned nsam,
                            double f = 0,
                            const char * fitness = "multiplicative"):
    """
    Evolve a set of populations, taking a random sample of size "nsam" every "sample" generations.

    :param rng: a :class:`GSLrng`
    :param pops: A list of populations simulated by :func:`evolve_regions`
    :param N: The diploid population size to simulate
    :param nlist: An array view of a NumPy array.  This represents the population sizes over time.  The length of this view is the length of the simulation in generations. The view must be of an array of 32 bit, unsigned integers (see example).
    :param mu_neutral: The mutation rate to variants not affecting fitness ("neutral" mutations).  The unit is per gamete, per generation.
    :param mu_selected: The mutation rate to variants affecting fitness ("selected" mutations).  The unit is per gamete, per generation.
    :param recrate: The recombination rate in the regions (per diploid, per generation)
    :param nregions: A list specifying where neutral mutations occur
    :param sregions: A list specifying where selected mutations occur
    :param recregions: A list specifying how the genetic map varies along the region
    :param sample: How often to record data about non-neutral mutations
    :param nsam: The sample size
    :param f: The selfing probabilty
    :param fitness: The fitness model.  Must be either "multiplicative" or "additive".

    :return: One list per population in "pops".  Each list contains a sample of genotypes + the generation from which it was taken.

    :raises: RuntimeError if parameters do not pass checks
    """
    check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    if sample <= 0:
        raise RuntimeError("sample must be > 0")
    if nsam == 0:
        raise RuntimeError("nsam must be > 0")
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    cdef size_t listlen = len(nlist)
    return evolve_regions_sample_async(rng.thisptr,&pops.pops,&nlist[0],listlen,mu_neutral,mu_selected,recrate,f,sample,nsam,rmgr.thisptr,fitness)

@cython.boundscheck(False)
def evolve_regions_track(GSLrng rng,
                            popvec pops,
                            unsigned[:] nlist,
                            double mu_neutral,
                            double mu_selected,
                            double recrate,
                            list nregions,
                            list sregions,
                            list recregions,
                            int sample,
                            double f = 0,
                            const char * fitness = "multiplicative"):
    """
    Evolve a set of populations, recording data on non-neutral every "sample" generations.  

    :param rng: a :class:`GSLrng`
    :param pops: A list of populations simulated by :func:`evolve_regions`
    :param N: The diploid population size to simulate
    :param nlist: An array view of a NumPy array.  This represents the population sizes over time.  The length of this view is the length of the simulation in generations. The view must be of an array of 32 bit, unsigned integers (see example).
    :param mu_neutral: The mutation rate to variants not affecting fitness ("neutral" mutations).  The unit is per gamete, per generation.
    :param mu_selected: The mutation rate to variants affecting fitness ("selected" mutations).  The unit is per gamete, per generation.
    :param recrate: The recombination rate in the regions (per diploid, per generation)
    :param nregions: A list specifying where neutral mutations occur
    :param sregions: A list specifying where selected mutations occur
    :param recregions: A list specifying how the genetic map varies along the region
    :param sample: How often to record data about non-neutral mutations
    :param f: The selfing probabilty
    :param fitness: The fitness model.  Must be either "multiplicative" or "additive".

    :return: One list per population in "pops".  Each list contains information on mutation frequency, effect size, etc., over time.

    :raises: RuntimeError if parameters do not pass checks
    """
    check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    if sample <= 0:
        raise RuntimeError("sample must be > 0")
    if f < 0.:
        warnings.warn("f < 0 will be treated as 0")
        f=0
    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    cdef size_t listlen = len(nlist)
    return evolve_regions_track_async(rng.thisptr,&pops.pops,&nlist[0],listlen,mu_neutral,mu_selected,recrate,f,sample,rmgr.thisptr,fitness)
    
@cython.boundscheck(False)
def evolve_regions_split(GSLrng rng,
                            popvec pops,
                            unsigned[:] nlist1,
                            unsigned[:] nlist2,
                            double mu_neutral,
                            double mu_selected,
                            double recrate,
                            list nregions,
                            list sregions,
                            list recregions,
                            vector[double] fs,
                            const char * fitness = "multiplicative"):
    """
    Take the output of a single-deme simulation, split into two demes, and evolve.

    :raises: RuntimeError if parameters do not pass checks

    Example:

    >>> #The first part is the same as the example for fwdpy.evolve_regions
    >>> import fwdpy
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> # Evolve for 5N generations initially
    >>> popsizes=np.tile(popsizes,100)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> #Now, "bud" off a daughter population of same size, and evolve both for another 100 generations
    >>> mpops = fwdpy.evolve_regions_split(rng,pops,popsizes[0:100],popsizes[0:100],0.001,0.0001,0.001,nregions,sregions,rregions,[0.,0.])

    """
    check_input_params(mu_neutral,mu_selected,recrate,nregions,sregions,recregions)
    cdef size_t i
    for i in range(fs.size()):
        if fs[i] < 0.:
            warnings.warn("f[i] < 0 will be treated as 0")
            fs[i]=0.0
    if len(nlist1) != len(nlist2):
        raise RuntimeError("len(nlist1) must equal len(nlist2)")

    mpv = mpopvec(0,[0])
    for i in range(len(pops)):
        #Step 1: Make the ith single pop the first deme in each metapop
        mpv.mpops.push_back(shared_ptr[metapop_t](new metapop_t([0])))
        re_init_mpop(mpv.mpops[i].get(),pops.pops[0].get())
        #Step 2: copy the first deme into a second deme
        copy_deme(mpv.mpops[i].get(),0)
        #Step 2a: make sure all is cool under the hood.
        if not mpv.mpops[i].get().sane():
            raise RuntimeError("fwdpy.evolve_regions_split: metapopulation data structures out of sync")
        #Step 3: make the python list have another shared_ptr to the metapop
        tt = metapop()
        tt.mpop = mpv.mpops[i]
        mpv.pympops.append(tt)

    rmgr = region_manager_wrapper()
    internal.make_region_manager(rmgr,nregions,sregions,recregions)
    cdef size_t listlen1 = len(nlist1)
    cdef size_t listlen2 = len(nlist2)
    with nogil:
        split_and_evolve_t(rng.thisptr,&mpv.mpops,&nlist1[0],listlen1,&nlist2[0],listlen2,mu_neutral,mu_selected,recrate,fs,rmgr.thisptr,fitness)
    return mpv

