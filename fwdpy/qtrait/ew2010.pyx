# def evolve_ewvw(GSLrng rng,
#                 int npops,
#                 int N,
#                 unsigned[:] nlist,
#                 double mu_neutral,
#                 double mu_selected,
#                 double recrate,
#                 list nregions,
#                 list sregions,
#                 list recregions,
#                 double sigmaE,
#                 double VS_total,
#                 double optimum = 0.,
#                 bint track = False,
#                 double f = 0.):
#     """
#     Eyre-Walker's 2010 model, modified so that the trait contributes to a portion of total variance in fitness.

#     :param rng: a :class:`GSLrng`
#     :param npops: The number of populations to simulate.  This is equal to the number of threads that will be used!
#     :param N: The diploid population size to simulate
#     :param nlist: An array view of a NumPy array.  This represents the population sizes over time.  The length of this view is the length of the simulation in generations. The view must be of an array of 32 bit, unsigned integers (see example).
#     :param mu_neutral: The mutation rate to variants not affecting fitness ("neutral" mutations).  The unit is per gamete, per generation.
#     :param mu_selected: The mutation rate to variants affecting fitness ("selected" mutations).  The unit is per gamete, per generation.
#     :param recrate: The recombination rate in the regions (per diploid, per generation)
#     :param nregions: A list specifying where neutral mutations occur
#     :param sregions: A list specifying where selected mutations occur
#     :param recregions: A list specifying how the genetic map varies along the region
#     :param sigmaE: Gaussian noise to add to fitness due to variation in "the trait".
#     :param VS_total: The total variance in fitness. This is the variance in fitness due to "the trait", plus all remaining variance. This parameter must be >= 1.
#     :param optimum: The optimum trait value.
#     :param track: whether or not to record the frequency trajectories of mutations.  True = simulations are much slower!
#     :param f: The selfing probabilty
#     """
#     if mu_neutral < 0.0:
#         raise RuntimeError("neutral mutation rate must be >= 0")
#     if mu_selected < 0.0:
#         raise RuntimeError("selected mutation rate must be >= 0")
#     if recrate < 0.0:
#         raise RuntimeError("recombination rate must be >= 0")
#     if sigmaE < 0.:
#         raise RuntimeError("sigmaE must be >= 0")
#     if VS_total < 1:
#         raise RuntimeError("VS_total must be >= 1")
#     pops = popvec(npops,N)
#     nreg = internal.process_regions(nregions)
#     sreg = internal.process_regions(sregions)
#     recreg = internal.process_regions(recregions)
#     v = shwrappervec()
#     internal.process_sregion_callbacks(v,sregions)
#     evolve_ewvw_t(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,VS_total,optimum,track,
#                     nreg['beg'].tolist(),nreg['end'].tolist(),nreg['weight'].tolist(),
#                     sreg['beg'].tolist(),sreg['end'].tolist(),sreg['weight'].tolist(),&v.vec,
#                     recreg['beg'].tolist(),recreg['end'].tolist(),recreg['weight'].tolist())
#     return pops

def ew2010_effects(GSLrng rng, singlepop pop, double tau, double sigma):
    """
    :param pop: A :class:`fwdpy.fwdpy.singlepop` simulated using :func:`fwdpy.qtrait.qtrait.evolve_qtrait` and/or :func:`fwdp.qtrait.qtrait.evolve_qtrait_more`
    :param tau: The coupling of trait value to fitness effect of mutation
    :param sigma: The standard deviation for Gaussian noise applied to trait value.  Generates the :math:`\epsilon` term in E-W's paper
    """
    if tau < 0.:
        raise RuntimeError("tau cannot be < 0.")
    if sigma < 0.:
        raise RuntimeError("sigma cannot be < 0.")
    cdef map[double,ew_mut_details] x = ew2010_assign_effects(rng.thisptr,pop.pop.get(),tau,sigma)
    cdef map[double,ew_mut_details].iterator itr = x.begin()
    rv = list()
    while itr != x.end():
        temp = deref(itr)
        rv.append( {temp.first : {'s':temp.second.s,'e':temp.second.e,'p':temp.second.p}} )
        inc(itr)
    return rv
    
def ew2010_traits(singlepop pop,list effects):
    """
    Implement model of Eyre-Walker 2010
    
    .. note:: The citation is www.pnas.org/cgi/doi/10.1073/pnas.0906182107.
       We implement the simple additive case here.
    """
    cdef map[double,ew_mut_details] temp
    cdef ew_mut_details tstruct
    for i in effects:
        for j in i:
            tstruct.s = i[j]['s']
            tstruct.e = i[j]['e']
            tstruct.p = i[j]['p']
            temp.insert( pair[double,ew_mut_details](j,tstruct) )
    return ew2010_traits_cpp(pop.pop.get(),temp)
