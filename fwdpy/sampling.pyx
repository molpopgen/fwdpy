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

def windows(vector[pair[double,string]] ms_sample, double windowSize,
            double stepSize, double startPos, double endPos):
    """
    Split a sample up into "windows" based on physical distance.

    :param ms_sample: The return value of a function like :func:`get_samples`
    :param windowSize: The length of each window, in same units as your simulation's regions.
    :param stepSize: The step size between windows, in same units as your simulation's regions.
    :param startPos: The minimum position possible (see Example below)
    :param endPos: The maximum position possible (see Example below)

    :return: A list of samples for each window.  Each element in the list is a list of tuples (the same type as the input data ms_sample).

    :rtype: A list
    
    .. note:: Empty elements in the return value represent windows with no variation.
    
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
    >>> #Get sample of size n = 10
    >>> s = fp.get_samples(rng,pops[0],10)
    >>> #Split the neutral variants in the sample up into non-overlapping windows of size 0.1
    >>> #The minimum position in all of nregions, sregions, and rregions is 0,
    >>> #and so 0 must be passed as 'startPos'.  Likewise, endPos must be 3.
    >>> #(If you were to input a value < 0., you'd get a bunch of empty windows in the return value.)
    >>> windows = fp.windows(s[0],0.1,0.1,0,3)
    """
    if windowSize <= 0.:
        raise RuntimeError("fwdpy.windows: windowSize must be > 0.")
    if stepSize <= 0.:
        raise RuntimeError("fwdpy.windows: stepSize must be > 0.")
    return sliding_windows_cpp(ms_sample,windowSize,stepSize,startPos,endPos)
