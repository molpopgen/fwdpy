# distutils: language = c++
# distutils: sources = src/sample.cpp src/neutral.cpp

from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

import threading

##Create hooks to C++ types

#Wrap the classes:
cdef extern from "types.hpp" namespace "fwdpy":
    cdef cppclass popvector:
        popvector(unsigned,unsigned)
    cdef cppclass singlepop_t:
        singlepop_t(unsigned)
    cdef cppclass GSLrng_t:
        GSLrng_t(unsigned)

##Now, wrap the functions
cdef extern from "neutral.hpp" namespace "fwdpy":
    void evolve_pop(GSLrng_t * rng, popvector * pops, const unsigned & ngens, const double & theta, const double & rho)


cdef extern from "sample.hpp" namespace "fwdpy":
    vector[vector[pair[double,string]]] take_sample_from_pop(GSLrng_t * rng,const popvector * pop,const unsigned & nsam)
    double tajd( const vector[pair[double,string]] & __data )
  
##Creat the python classes
cdef class popvec:
    """
    Vector of single-deme objects

    The constructor takes two objects: the number of pops, and the initial population size (which is the same for each pop)
    """
    cdef popvector *thisptr
    def __cinit__(self,unsigned npops,unsigned N):
        self.thisptr = new popvector(npops,N)
    def __dealloc__(self):
        del self.thisptr
    
cdef class GSLrng:
    """
    A wrapper around a random number generator (rng) 
    from the GNU Scientific Library (GSL).

    The constructor takes a seed (int) as an argument.

    Example:

    >>> import fwdpy
    >>> rng = fwdpy.GSLrng(100)
    """
    cdef GSLrng_t * thisptr
    def __cinit__(self, int seed):
        self.thisptr = new GSLrng_t(seed)
    def __dealloc__(self):
        del self.thisptr

def evolve_pops_t(GSLrng rng,int npops, int N, int ngens, double theta, double rho):
    p=popvec(npops,N)
    #call the C++ fxn
    evolve_pop(rng.thisptr,p.thisptr,ngens,theta,rho)
    return p

def ms_sample(GSLrng rng, popvec pops, int nsam):
    return take_sample_from_pop(rng.thisptr,pops.thisptr,nsam)
    
def TajimasD( vector[pair[double,string]] data ):
    """
    Calculate Tajima's D statistic from a sample

    :param data: a sample from a population

    Example:
    
    >>> import fwdpy
    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve(rng,1000,1000,50,50)
    >>> s = fwdpy.ms_sample(rng,pop,10)
    >>> d = fwdpy.TajimasD(s)
    """
    return tajd(data)



    
