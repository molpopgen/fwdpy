def evolve_pops_t(GSLrng rng,int npops, int N, list nlist, double theta, double rho):
    """
    Simulate populations under idealized Wright-Fisher model

    :param rng: a GSLrng
    :param npops: The number of populations to simulate.  This is equal to the number of threads that will be used!
    :param N: The diploid population size to simulate
    :param nlist: A list of population sizes over time.  This allows the simulation of arbitrary changes in N
    :param theta: :math:`\\theta = 4N_e\\mu` is the scaled mutation rate to variants not affecting fitness ("neutral mutations")
    :param rho: :math:`\\rho = 4N_er` is the scaled recombination rate

    Example:
    
    >>> import fwdpy
    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*1000,50,50)
    """
    p=popvec(npops,N)
    #call the C++ fxn
    evolve_pop(rng.thisptr,&p.pops,nlist,theta,rho)
    return p

def evolve_pops_more_t(GSLrng rng, popvec pops, list nlist, double theta, double rho):
    """
    Continue evolving a set of populations

    :param rng: a GSLrng
    :param pops: A list of populations simulated using evolve_pops_t
    :param nlist: A list of population sizes over time.  This allows the simulation of arbitrary changes in N
    :param theta: :math:`\\theta = 4N_e\\mu` is the scaled mutation rate to variants not affecting fitness ("neutral mutations")
    :param rho: :math:`\\rho = 4N_er` is the scaled recombination rate
    """
    evolve_pop(rng.thisptr,&pops.pops,nlist,theta,rho);
