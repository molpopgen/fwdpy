def summstats(vector[pair[double,string]] data ):
    """
    Calculate common summary statistics from a sample

    :param data: a sample from a population.  For example, the return value of :func:`fwdpy.fwdpy.ms_sample` or :func:`fwdpy.fwdpy.get_samples`

    :raises: RuntimeError if data are not sorted in increasing order according to position
    
    Example:
    
    >>> import fwdpy
    >>> import fwdpy.libseq as lseq
    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*1000,50,50)
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    >>> d = [lseq.summstats(si) for si in s]
    """
    return libseq_basic_stats(data)

def hapstats(vector[pair[double,string]] data,
             double minfreq,
             double binsize,
             gmap = None):
    """
    Calculate haplotype-based statistics from a sample

    :param data: a sample from a population.  For example, the return value of :func:`fwdpy.fwdpy.ms_sample` or :func:`fwdpy.fwdpy.get_samples`
    :param minfreq: exclude sites with minor allele frequency < minfreq
    :param binsize: The frequency window size for binning
    :param gmap: A list of values containing the position of each variant in the data's location on the genetic map.
    """
    if gmap is None:
        return libseq_extra_ld_stats(data,minfreq,binsize)
    else:
        return libseq_extra_ld_stats(data,minfreq,binsize,gmap)
    
