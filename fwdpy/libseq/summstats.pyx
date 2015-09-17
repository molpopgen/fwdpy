def TajimasD( vector[pair[double,string]] data ):
    """
    Calculate Tajima's D statistic from a sample

    :param data: a sample from a population.  The return value of :func:`ms_sample`.

    Example:
    
    >>> import fwdpy
    >>> import fwdpy.libseq as lseq
    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*1000,50,50)
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    >>> d = [lseq.TajimasD(si) for si in s]
    """
    return tajd(data)
