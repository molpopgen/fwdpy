def trajectories( singlepop pop, int minsojourn = 0, double minfreq = 0.):
    """
    Get the frequency trajactories of mutations

    :param pop: A :class:`fwdpy.fwdpy.singlepop` simulated using :func:`fwdpy.fwdpy.evolve_regions` and/or :func:`fwdpy.fwdpy.evolve_regions_more`
    :param minsojounrn: Exclude all mutations that did not exist in the population for :math:`\geq\ \mathrm{minsojourn}` generations
    :param minfreq: Exclude all mutations that did not read a minimum frequency of :math:`\geq\ \mathrm{minfreq}`

    :return: A pandas.DataFrame

    :rtype: pandas.DataFrame
    """
    return pandas.DataFrame.from_dict(get_singlepop_traj(pop.pop.get(),minsojourn,minfreq))
