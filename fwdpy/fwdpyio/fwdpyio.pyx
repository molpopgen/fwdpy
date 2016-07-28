# distutils: language = c++
# distutils: sources = fwdpy/fwdpyio/serialize.cc
    
##Undocumented fxns are implementation details
def serialize_single(Spop pop):
    return serialize_singlepop(pop.pop.get())

def serialize_meta(MetaPop mpop):
    return serialize_metapop(mpop.mpop.get())

def serialize_mlocus(MlocusPop pop):
    return serialize_multilocus(pop.pop.get())

def serialize(PopType pop):
    """
    Return a binary representation of an evolved population

    :param pop: A list of :class:`fwdpy.fwdpy.PopType`
    
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
    if isinstance(pop,Spop):
        return serialize_single(pop)
    elif isinstance(pop,MetaPop):
        return serialize_meta(pop)
    elif isinstance(pop,MlocusPop):
        return serialize_mlocus(pop)
    else:
        raise RuntimeError("fwdpyio.serialize: unsupported PopType "+str(type(pop)))

def deserialize_singlepops(list strings):
    """
    Convert binary representation back to a :class:`fwdpy.fwdpy.PopVec`

    :param strings: A list of populations in binary format.  This should be the value returned by :func:`fwdpy.fwdpyio.fwdpyio.serialize`

    :returns: :func:`fwdpy.fwdpy.PopVec`

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
    pops=SpopVec(0,0)
    pops.reset(temp)
    return pops

def deserialize_metapops(list strings):
    """
    Convert binary representation of populations back to a :class:`fwdpy.fwdpy.mpopvec`

    :param strings: A list of populations in binary format.  This should be the value returned by :func:`fwdpy.fwdpyio.fwdpyio.serialize`

    :returns: :func:`fwdpy.fwdpy.mpopvec`

    .. note:: len(strings) determines the length of the return value, and therefore the number of threads to use if the population is evolved further.

    Example:

    TODO
    """
    cdef vector[shared_ptr[metapop_t]] temp = deserialize_metapop(strings)
    mpops = MetaPopVec(0,[0]*1)
    mpops.reset(temp)
    return mpops
