##Undocumented fxns are wrappers to enable run-time polymorphism within the Py environs.
##These fxns make calls to the C++ layer

def ms_sample_single_deme(GSLrng rng, singlepop pop, int nsam, bint removeFixed):
    return take_sample_from_pop(rng.thisptr,pop.pop.get(),nsam, int(removeFixed))

def ms_sample_single_deme_sep(GSLrng rng, singlepop pop, int nsam, bint removeFixed):
    return take_sample_from_pop_sep(rng.thisptr,pop.pop.get(),nsam, int(removeFixed))

def ms_sample_metapop_sep(GSLrng rng, metapop pop, int nsam, bint removeFixed,int deme):
    if deme >= len(pop):
        raise RuntimeError("value for deme out of range. len(pop) = "+str(len(pop))+", but deme = "+str(deme))
    return take_sample_from_metapop_sep(rng.thisptr,pop.mpop.get(),nsam, int(removeFixed), deme)

def diploid_view_singlepop(singlepop pop, int ind, bint removeFixed):
    return diploid_view_cpp(pop.pop.get(),ind,removeFixed)

def diploid_view_metapop(metapop mpop, int ind, bint removeFixed, int deme):
    return diploid_view_cpp(mpop.mpop.get(),ind,removeFixed,deme)

cdef get_sh_single(const vector[pair[double,string] ] & ms_sample,
                    singlepop pop,
                    vector[double] * s,
                    vector[double] * h,
                    vector[double] * p,
                    vector[double] * a):
    get_sh(ms_sample,pop.pop.get(),s,h,p,a)

cdef get_sh_metapop(const vector[pair[double,string] ] & ms_sample,
                    metapop pop,
                    vector[double] * s,
                    vector[double] * h,
                    vector[double] * p,
                    vector[double] * a):
    get_sh(ms_sample,pop.mpop.get(),s,h,p,a)  

def ms_sample(GSLrng rng, poptype pop, int nsam, bint removeFixed = True):
    """
    Take a sample from a set of simulated populations.

    :param rng: a :class:`GSLrng`
    :param pops: An object inheriting from :class:`poptype`
    :param nsam: List of sample sizes (no. chromosomes) to sample.
    :param removeFixed: if True, only polymorphic sites are retained

    Note: nsam will likely be changed to a list soon, to accomodate multi-deme simulations
    
    Example:
    
    >>> import fwdpy
    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*1000,50,50)
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    """
    if isinstance(pop,singlepop):
        return ms_sample_single_deme(rng,pop,nsam,removeFixed)
    else:
        raise ValueError("ms_sample: unsupported type of popcontainer")

def get_samples(GSLrng rng, poptype pop, int nsam, bint removeFixed = True, deme = None):
    """
    Take a sample from a set of simulated populations.

    :param rng: a :class:`GSLrng`
    :param pop: An object inheriting from :class:`poptype`
    :param nsam: The sample size to take.
    :param removeFixed: if True, only polymorphic sites are retained
    :param deme: Optional.  If 'pop' is a :class:`metapop`, deme is required and represents the sub-population to sample.
    
    :return: A list. Element 0 is neutral mutations, and element 1 is selected mutations.  Within each list is a tuple of size 2.  The first element is the mutation position.  The second element is the genotype for each of the 'nsam' chromosomes.  Genotypes are coded as 0 = the ancestral state and 1 = the derived state.  For each site, each pair of genotypes constitutes a single diploid.  In other words, for nsam = 50, the data will represent the complete haplotypes of 25 diploids.

    Please note that if you desire an odd 'nsam', you should input nsam+2 and randomly remove one haplotype to obtain your desired sample size.  This is due to an issue with how we are sampling chromosomes from the population.

    Example:
    
    >>> import fwdpy
    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*1000,50,50)
    >>> s = [fwdpy.get_samples(rng,i,10) for i in pop]
    """
    if isinstance(pop,singlepop):
        return ms_sample_single_deme_sep(rng,pop,nsam,removeFixed)
    elif isinstance(pop,metapop):
        if deme is None:
            raise RuntimeError("deme may not be set to None when sampling from a meta-population")
        return ms_sample_metapop_sep(rng,pop,nsam,removeFixed,deme)
    else:
        raise ValueError("ms_sample: unsupported type of popcontainer")

def get_sample_details( vector[pair[double,string]] ms_sample, poptype pop ):
    """
    Get additional details for population samples

    :param ms_samples: A list returned by :func:`ms_sample`
    :param pops: A :class:`poptype`

    :return: A pandas.DataFrame containing the selection coefficient (s), dominance (h), populations frequency (p), and age (a) for each mutation.

    :rtype: pandas.DataFrame

    Example:
    
    >>> import fwdpy
    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*1000,50,50)
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    >>> details = [fwdpy.get_sample_details(i,j) for i,j in zip(s,pop)]
    """
    cdef vector[double] h
    cdef vector[double] s
    cdef vector[double] p
    cdef vector[double] a
    if isinstance(pop,singlepop):
        get_sh_single(ms_sample,pop,&s,&h,&p,&a)
    elif isinstance(pop,metapop):
        get_sh_metapop(ms_sample,pop,&s,&h,&p,&a)
    return pandas.DataFrame({'s':s,'h':h,'p':p,'a':a})

def diploid_view( poptype pop, list indlist, bint removeFixed = False, deme = None ):
    """
    Get detailed information about a list of diploids.

    :param pop: A :class:`poptype`
    :param indlist: A list of *indexes* of individuals to sample. (Start counting from 0.)
    :param removeFixed: If non-zero, fixations will be excluded.
    :param deme: If pop is of type :class:`metapop`, deme is the index of the sub-population from which to get the individuals

    :return: A pandas.DataFrame containing information about each mutation for each individual in indlist.

    :rtype: pandas.DataFrame

    :raises: IndexError if any item in indlist is out of range, or if deme is out of range.

    .. note:: This return value of this function does not allow the calculation of fixation times.
       In order to do that, a change must be made to fwdpp, which may or may not happen
       soon.

    Example:

    >>> import fwdpy as fp
    >>> import numpy as np
    >>> rng = fp.GSLrng(100)
    >>> nregions = [fp.Region(0,1,1),fp.Region(2,3,1)]
    >>> sregions = [fp.ExpS(1,2,1,-0.1),fp.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fp.Region(0,3,1)]
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> # Evolve for 5N generations initially
    >>> popsizes=np.tile(popsizes,10000)
    >>> pops = fp.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> #Take a "view" of the first 5 diploids:
    >>> view = fp.diploid_view(pops[0],[0,1,2,3,4])
    """
    if isinstance(pop,singlepop):
        rv = pandas.DataFrame()
        for i in indlist:
            temp = diploid_view_singlepop(pop,i,removeFixed)
            rv = pandas.concat([rv,pandas.DataFrame.from_dict(temp)])
        return rv
    elif isinstance(pop,metapop):
        if deme is None:
            raise RuntimeError("deme may not be set to None when taking a view from a meta-population")
        rv = pandas.DataFrame()
        for i in indlist:
            temp = diploid_view_metapop(pop,i,removeFixed)
            rv = pandas.concat([rv,pandas.DataFrame.from_dict(temp)])
        return rv
    else:
        raise ValueError("diploid_view: type of pop is not supported")

###### Functions for manipulating samples.
def nderived_site(tuple site):
    """
    Get the number of derived mutations at a site.

    :param site: A tuple.  See example

    .. note:: In general, it will be more convenient to call :func:`fwdpy.fwdpy.nderived` on a list of tuples.
    
    Example:

    >>> import fwdpy
    >>> #Create a site at position 0.1 with the
    >>> #genotypes as the second element. 0/1 = ancestral/derived
    >>> x = (0.1,'01100111')
    >>> fwdpy.nderived_site(x)
    5
    >>> #Process simulation results:
    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve_pops_t(rng,1,1000,[1000]*1000,50,50)
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    >>> for i in s[0]: ndi = fwdpy.nderived_site(i)
    """
    return site[1].count('1')

def nderived( list sample ):
    """
    Convenience wrapper around :func:`fwdpy.fwdpy.nderived`

    :param sample: a sample from a population.  For example, the return value of :func:`fwdpy.fwdpy.ms_sample` or :func:`fwdpy.fwdpy.get_samples`

    Example:

    >>> import fwdpy

    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*1000,50,50)
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    >>> nd = [fwdpy.nderived(i) for i in s]
    """
    return [nderived_site(i) for i in sample]

def getfreq(tuple site,bint derived = True):
    """
    Get mutation frequencies

    :param site: A tuple. See example
    :param derived:  If True, report derived allele frequency (DAF).  If False, return minor allele freqency (MAF)

    .. note:: Do **not** use this function to calculate :math:`\pi` (a.k.a. :math:`\\hat\\theta_\pi`, a.k.a. "sum of site heterozygosity").
       :math:`\pi` for a **sample** is not :math:`2\sum_ip_iq_i`. because the sample is *finite*.  Please use :func:`fwdpy.libseq.libseq.summstats` instead.
       In general, it will be more convenient to call :func:`fwdpy.fwdpy.getfreqs` on a list of tuples.
       
    Example:

    >>> import fwdpy
    >>> #DAF = 7/10
    >>> #MAF = 3/10
    >>> x = (0.1,'0111111100')
    >>> round(fwdpy.getfreq(x,True),3)
    0.7
    >>> round(fwdpy.getfreq(x,False),3)
    0.3
    """
    o = nderived_site(site)
    dfreq = float(o)/float(len(site[1]))
    if derived is True:
        return dfreq
    return min(dfreq,1.0-dfreq)

def getfreqs(list sample,bint derived = True):
    """
    Convenience wrapper around :func:`fwdpy.fwdpy.getfreq`

    :param sample: a sample from a population.  For example, the return value of :func:`fwdpy.fwdpy.ms_sample` or :func:`fwdpy.fwdpy.get_samples`
    :param derived: If True, report derived allele frequency (DAF).  If False, return minor allele freqency (MAF).
    
    Example:

    >>> import fwdpy

    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*1000,50,50)
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    >>> freqs = [fwdpy.getfreqs(i) for i in s]
    """
    return [getfreq(i,derived) for i in sample]

def freqfilter( list sample,
                float minfreq,
                bint derived = True ):
    """
    Remove low-frequency variants from a sample.

    :param sample: a sample from a population.  For example, the return value of :func:`fwdpy.fwdpy.ms_sample` or :func:`fwdpy.fwdpy.get_samples`
    :param minfreq: Remove all sites with frequency < minfreq
    :param derived: if True, filter on derived allele frequency.  If False, filter on minor allele frequency.

    .. note:: Do **not** use this function to calculate :math:`\pi` (a.k.a. :math:`\\hat\\theta_\\pi`, a.k.a. "sum of site heterozygosity").
       :math:`\pi` for a **sample** is not :math:`2\sum_ip_iq_i`. because the sample is *finite*.  Please use :func:`fwdpy.libseq.libseq.summstats` instead.
       
    Example:
    
    >>> import fwdpy
    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*1000,50,50)
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    >>> s2 = [fwdpy.freqfilter(i,0.2) for i in s]
    """
    rv=list()
    for i in sample:
        if type(i) is not tuple:
            raise RuntimeError("values is sample must be tuples, not"+str(type(i)))
        if getfreq(i,derived)>=minfreq:
            rv.append(i)
    return rv
