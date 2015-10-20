def trajectories( singlepop pop, int minsojourn = 0, double minfreq = 0.):
    """
    Get the frequency trajectories of (non-neutral) mutations.

    :param pop: A :class:`fwdpy.fwdpy.singlepop` simulated using :func:`fwdpy.fwdpy.evolve_regions` and/or :func:`fwdpy.fwdpy.evolve_regions_more`
    :param minsojounrn: Exclude all mutations that did not exist in the population for :math:`\geq\ \mathrm{minsojourn}` generations
    :param minfreq: Exclude all mutations that did not read a minimum frequency of :math:`\geq\ \mathrm{minfreq}`

    :return: A pandas.DataFrame

    :rtype: pandas.DataFrame

    Example:

    >>> #Same example as fwdpy.evolve_regions
    >>> import fwdpy
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> popsizes=np.tile(popsizes,100)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions,track=True)
    >>> traj = [fwdpy.trajectories(i) for i in pops]
    """
    return pandas.DataFrame.from_dict(get_singlepop_traj(pop.pop.get(),minsojourn,minfreq))
