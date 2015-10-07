def popstats( singlepop pop ):
    """
    Get statistics like VG, etc. from a population

    :return: A pandas.DataFrame.

    :rtype: pandas.DataFrame
    """
    return pandas.DataFrame(qtrait_pop_props(pop.pop.get()).items(),columns=['stat','value'])

def esize_freq(singlepop pop):
    """
    Returns effect size vs frequency

    :param pop: A :class:`fwdpy.fwdpy.singlepop` simulated using :func:`fwdpy.qtrait.qtrait.evolve_qtrait` and/or :func:`fwdp.qtrait.qtrait.evolve_qtrait_more`

    :rtype: pandas.DataFrame
    """
    return pandas.DataFrame.from_dict(qtrait_esize_freq(pop.pop.get()))
