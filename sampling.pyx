def ms_sample(GSLrng rng, popvec pops, int nsam):
    """
    Take a sample from a set of simulated populations.

    :param rng: a GSLrng
    :param pops: a vector of populations ("popvec")
    :param nsam: the sample size (no. chromosomes) to sample

    Example:
    
    >>> import fwdpy
    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*1000,50,50)
    >>> s = fwdpy.ms_sample(rng,pop,10)
    """
    return take_sample_from_pop(rng.thisptr,pops.pops,nsam)

def get_sample_details(list ms_samples, popvec pops):
    """
    Get additional details for population samples

    :param ms_samples: A list returned by fwdpy.ms_sample
    :param pops: A a popvec containing simulated populations

    :return: A pandas.DataFrame containing the selection coefficient (s), dominance (h), populations frequency (p), and age (a) for each mutation.

    :rtype: pandas.DataFrame
    """
    rv = list()
    cdef vector[double] h
    cdef vector[double] s
    cdef vector[double] p
    cdef vector[double] a
    for i in range(len(ms_samples)):
        s.clear()
        h.clear()
        p.clear()
        a.clear()
        get_sh(ms_samples,pops.pops,i,&s,&h,&p,&a)
        rv.append( pandas.DataFrame({'s':s,'h':h,'p':p,'a':a}) )
    return rv

def TajimasD( vector[pair[double,string]] data ):
    """
    Calculate Tajima's D statistic from a sample

    :param data: a sample from a population

    Example:
    
    >>> import fwdpy
    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*1000,50,50)
    >>> s = fwdpy.ms_sample(rng,pop,10)
    >>> d = [fwdpy.TajimasD(si) for si in s]
    """
    return tajd(data)
