
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
    >>> import fwdpy.libseq as lseq
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
    >>> windows = lseq.windows(s[0],0.1,0.1,0,3)
    """
    if windowSize <= 0.:
        raise RuntimeError("fwdpy.windows: windowSize must be > 0.")
    if stepSize <= 0.:
        raise RuntimeError("fwdpy.windows: stepSize must be > 0.")
    return sliding_windows_cpp(ms_sample,windowSize,stepSize,startPos,endPos)
