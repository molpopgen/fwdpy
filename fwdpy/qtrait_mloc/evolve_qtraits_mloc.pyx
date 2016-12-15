from fwdpy.fitness cimport MlocusFitness
from fwdpy.internal import process_sregion_callbacks,make_region_manager

def evolve_qtraits_mloc_sample_fitness(GSLrng rng,
                                       MlocusPopVec pops,
                                       TemporalSampler slist,
                                       MlocusFitness fitness_function,
                                       unsigned[:] nlist,
                                       const vector[double] & mu_neutral,
                                       const vector[double] & mu_selected,
                                       list sregions,
                                       const vector[double] & recrates_within,
                                       const vector[double] & recrates_between,
                                       int sample,
                                       double optimum = 0.0,
                                       double sigmaE = 0.0,
                                       double f = 0.0,
                                       double VS = 1.0):
    if sample<0:
        raise RuntimeError("sample must be >= 0")
    cdef size_t nlen=len(nlist)
    sh = shwrappervec()
    process_sregion_callbacks(sh,sregions)
    evolve_qtrait_mloc_cpp(rng.thisptr,&pops.pops,slist.vec,
                           &nlist[0],nlen,mu_neutral,mu_selected,
                           sh.vec,
                           recrates_within,
                           recrates_between,f,sigmaE,optimum,VS,sample,
                           fitness_function.wfxn)

def evolve_qtraits_mloc_regions_sample_fitness(GSLrng rng,
                                       MlocusPopVec pops,
                                       TemporalSampler slist,
                                       MlocusFitness fitness_function,
                                       unsigned[:] nlist,
                                       list nregions,
                                       list sregions,
                                       list recregions,
                                       const vector[double] & recrates_between,
                                       int sample,
                                       double optimum = 0.0,
                                       double sigmaE = 0.0,
                                       double f = 0.0,
                                       double VS = 1.0):
    if sample<0:
        raise RuntimeError("sample must be >= 0")
    if recrates_between.size() != len(nregions)-1:
        raise RuntimeError("There must be i-1 between-locus crossover rates for an i-locus simulation")

    cdef size_t nlen=len(nlist)
    rmgr = region_manager_wrapper()
    make_region_manager(rmgr,nregions,sregions,recregions)
    evolve_qtrait_mloc_regions_cpp(rng.thisptr,&pops.pops,slist.vec,
                           &nlist[0],nlen,rmgr.thisptr,
                           recrates_between,f,sigmaE,optimum,VS,sample,
                           fitness_function.wfxn)
