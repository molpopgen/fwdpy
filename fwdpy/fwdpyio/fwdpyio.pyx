# distutils: language = c++
# distutils: sources = fwdpy/fwdpyio/serialize.cc

from libcpp.string cimport string
from fwdpy.fwdpy cimport *

cdef extern from "serialize.hpp" namespace "fwdpy::serialize":
    string serialize_singlepop(const singlepop_t * pop)
    vector[shared_ptr[singlepop_t]] deserialize_singlepop(const vector[string] & strings)
    
##Undocumented fxns are implementation details
def serialize_single(singlepop pop):
    return serialize_singlepop(pop.pop.get())

def serialize(poptype pop):
    """
    Return a binary representation of an evolved population

    Example:

    >>> import fwdpy
    >>> import fwdpy.io as fpio
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> popsizes=np.tile(popsizes,10000)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> strings = [fpio.serialize(i) for i in pops]
    """
    if isinstance(pop,singlepop):
        return serialize_single(pop)


def deserialize_singlepops(list strings):
    cdef vector[shared_ptr[singlepop_t]] test = deserialize_singlepop(strings)
    pops=popvec(0,0)
    pops.reset(test)
    return pops
