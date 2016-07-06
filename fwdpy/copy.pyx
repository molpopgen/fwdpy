import fwdpyio as fpio

def copypop(PopType pop):
    """
    Copy a population

    :param pop: the population to copy.

    :rtype: A :class:`fwdpy.fwdpy.PopType` of the same (derived) type as pop

    .. note:: The function is implemented in terms of the serialization routines in the fwdpy.fwdpyio module.  The return value can be evolved and not affect the input value.
    """
    s = fpio.serialize(pop)
    if isinstance(pop,Spop):
        return fpio.deserialize_singlepops(list(s))
    elif isinstance(pop,MetaPop):
        return fpio.deserialize_metapops(list(s))
    else:
        raise RuntimeError("fwdpy.copypop: PopType "+str(type(pop))+" is not supported")

def copypops(PopVec pops):
    """
    Copy a population

    :param pops: the list of population to copy.

    :rtype: A :class:`fwdpy.fwdpy.PopVec` of the same (derived) type as pop

    .. note:: The function is implemented in terms of the serialization routines in the fwdpy.fwdpyio module.  The return value can be evolved and not affect the input value.
    """
    s = [fpio.serialize(pop) for pop in pops]
    if isinstance(pops,SpopVec):
        return fpio.deserialize_singlepops(s)
    elif isinstance(pops,MetaPopVec):
        return fpio.deserialize_metapops(s)
    else:
        raise RuntimeError("fwdpy.copypopvec: popvec type "+str(type(pops))+" is not supported")
