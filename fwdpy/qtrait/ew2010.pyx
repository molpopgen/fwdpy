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
