def summstats(vector[pair[double,string]] data ):
    """
    Calculate common summary statistics from a sample

    :param data: a sample from a population.  For example, the return value of :func:`fwdpy.fwdpy.ms_sample` or :func:`fwdpy.fwdpy.get_samples`
    
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

    :raises: RuntimeError if data are not sorted in increasing order according to position
    
    Example:
    
    >>> import fwdpy
    >>> import fwdpy.libseq as lseq
    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*1000,50,50)
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    >>> d = [lseq.hapstats(si,0.05,0.05) for si in s]
    """
    if gmap is None:
        return libseq_extra_ld_stats(data,minfreq,binsize)
    else:
        return libseq_extra_ld_stats(data,minfreq,binsize,gmap)

def lHAF( vector[pair[double,string]] data, double l ):
    """
    Calculate :math:`l-HAF` from doi:10.1371/journal.pgen.1005527.g001

    :param data: A sample from a population.  See :func:`fwdpy.fwdpy.get_samples` for details
    :param l: the exponent i :math:`l-HAF`

    :rtype: A list of values, which is the :math:`l-HAF` score for each haplotype 
    """
    return lhaf(data,l)
    
