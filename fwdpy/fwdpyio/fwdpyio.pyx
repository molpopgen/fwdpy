# distutils: language = c++
# distutils: sources = fwdpy/fwdpyio/serialize.cc
    
##Undocumented fxns are implementation details
def serialize_single(singlepop pop):
    return serialize_singlepop(pop.pop.get())

def serialize_meta(metapop mpop):
    return serialize_metapop(mpop.mpop.get())

def serialize(poptype pop):
    """
    Return a binary representation of an evolved population

    :param pop: A list of :class:`fwdpy.fwdpy.poptype`
    
    Example:

    >>> import fwdpy
    >>> import fwdpy.fwdpyio as fpio
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> popsizes=np.tile(popsizes,100)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> strings = [fpio.serialize(i) for i in pops]
    """
    if isinstance(pop,singlepop):
        return serialize_single(pop)
    elif isinstance(pop,metapop):
        return serialize_meta(pop)
    else:
        raise RuntimeError("fwdpyio.serialize: unsupported poptype "+str(type(pop)))

def deserialize_singlepops(list strings):
    """
    Convert binary representation back to a :class:`fwdpy.fwdpy.popvec`

    :param strings: A list of populations in binary format.  This should be the value returned by :func:`fwdpy.fwdpyio.fwdpyio.serialize`

    :returns: :func:`fwdpy.fwdpy.popvec`

    .. note:: len(strings) determines the length of the return value, and therefore the number of threads to use if the population is evolved further.
        
    Example:
    
    >>> import fwdpy
    >>> import fwdpy.fwdpyio as fpio
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> popsizes=np.tile(popsizes,100)
    >>> pops = fwdpy.evolve_regions(rng,4,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> strings = [fpio.serialize(i) for i in pops]
    >>> len(strings)
    4
    >>> pops2 = fpio.deserialize_singlepops(strings)
    """
    cdef vector[shared_ptr[singlepop_t]] temp = deserialize_singlepop(strings)
    pops=popvec(0,0)
    pops.reset(temp)
    return pops

def deserialize_metapops(list strings):
    """
    Convert binary representation of populations back to a :class:`fwdpy.fwdpy.mpopvec`

    :param strings: A list of populations in binary format.  This should be the value returned by :func:`fwdpy.fwdpyio.fwdpyio.serialize`

    :returns: :func:`fwdpy.fwdpy.mpopvec`

    .. note:: len(strings) determines the length of the return value, and therefore the number of threads to use if the population is evolved further.

    Example:

    >>> #The first part is the same as the example for :func:`fwdpy.fwdpy.evolve_regions_split`
    >>> import fwdpy
    >>> import fwdpy.fwdpyio as fpio
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> # Evolve for 5N generations initially
    >>> popsizes=np.tile(popsizes,100)
    >>> pops = fwdpy.evolve_regions(rng,4,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> #Now, "bud" off a daughter population of same size, and evolve both for another 100 generations
    >>> mpops = fwdpy.evolve_regions_split(rng,pops,popsizes[0:100],popsizes[0:100],0.001,0.0001,0.001,nregions,sregions,rregions,[0.,0.])
    >>> #Serialize
    >>> bstrings = [fpio.serialize(i) for i in mpops]
    >>> len(bstrings)
    4
    >>> #Deserialize
    >>> mpops2 = fpio.deserialize_metapops(bstrings)
    """
    cdef vector[shared_ptr[metapop_t]] temp = deserialize_metapop(strings)
    mpops = mpopvec(0,[0]*1)
    mpops.reset(temp)
    return mpops
