def getmuts(poptype pop, bint all = True, bint fixations = False):
    """
    Get summary of mutations in population

    :param pop: A :class:`fwdpy.fwdpy.poptype`
    :param all: if True, include all mutations.  Otherwise, only include non-neutral variants
    :param fixations: if True, include fixations

    :rtype: pandas.DataFrame
    """
    if isinstance(pop,singlepop):
        return internal.getmuts_singlepop(pop,all,fixations)
    elif isinstance(pop,metapop):
        return internal.getmuts_metapop(pop,all,fixations)
    else:
        raise RuntimeError("fwdpy.getmuts: unsupported poptype "+str(type(pop)))
