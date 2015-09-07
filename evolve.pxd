from libcpp.vector cimport vector

include "fwdpy.pxi"
include "fwdpy.pyx"
        
##Now, wrap the functions
cdef extern from "neutral.hpp" namespace "fwdpy":
    void evolve_pop(GSLrng_t * rng, vector[shared_ptr[singlepop_t]] * pops, const vector[unsigned] nlist, const double & theta, const double & rho)

def evolve_pops_more_t(GSLrng rng, popvec pops, list nlist, double theta, double rho):
    """
    Continue evolving a set of populations

    :param rng: a GSLrng
    :param pops: A list of populations simulated using evolve_pops_t
    :param nlist: A list of population sizes over time.  This allows the simulation of arbitrary changes in N
    :param theta: :math:`\\theta = 4N_e\\mu` is the scaled mutation rate to variants not affecting fitness ("neutral mutations")
    :param rho: :math:`\\rho = 4N_er` is the scaled recombination rate
    """
    evolve_pop(rng.thisptr,pops.pops,nlist,theta,rho);
