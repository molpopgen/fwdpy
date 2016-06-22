from fwdpy.fitness cimport MlocusFitness
from fwdpy.internal import process_sregion_callbacks

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
    cdef size_t nlen=len(nlist)
    sh = shwrappervec()
    process_sregion_callbacks(sh,sregions)
    evolve_qtrait_mloc_cpp(rng.thisptr,&pops.pops,slist.vec,
                           &nlist[0],nlen,mu_neutral,mu_selected,
                           sh.vec,
                           recrates_within,
                           recrates_between,f,sigmaE,optimum,VS,sample,
                           fitness_function.wfxn)
