def popstats( singlepop pop ):
    """
    Get statistics like VG, etc. from a population

    :return: A pandas.DataFrame.

    :rtype: pandas.DataFrame
    """
    return pandas.DataFrame(qtrait_pop_props(pop.pop.get()).items(),columns=['stat','value'])

def trajectories( singlepop pop, int minsojourn = 0, double minfreq = 0.):
    """
    Get the frequency trajactories of mutations

    :param pop: A :class:`fwdpy.fwdpy.singlepop` simulated using :func:`fwdpy.qtrait.qtrait.evolve_qtrait` and/or :func:`fwdp.qtrait.qtrait.evolve_qtrait_more`
    :param minsojounrn: Exclude all mutations that did not exist in the population for :math:`\geq\ \mathrm{minsojourn}` generations
    :param minfreq: Exclude all mutations that did not read a minimum frequency of :math:`\geq\ \mathrm{minfreq}`

    :return: A pandas.DataFrame

    :rtype: pandas.DataFrame
    """
    return pandas.DataFrame.from_dict(get_qtrait_traj(pop.pop.get(),minsojourn,minfreq))

def esize_freq(singlepop pop):
    """
    Returns effect size vs frequency

    :param pop: A :class:`fwdpy.fwdpy.singlepop` simulated using :func:`fwdpy.qtrait.qtrait.evolve_qtrait` and/or :func:`fwdp.qtrait.qtrait.evolve_qtrait_more`

    :rtype: pandas.DataFrame
    """
    return pandas.DataFrame.from_dict(qtrait_esize_freq(pop.pop.get()))
