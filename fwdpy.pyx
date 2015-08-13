# distutils: language = c++
# distutils: sources = src/sample.cpp src/neutral.cpp

from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

import threading

##Create hooks to C++ types

#Wrap the classes:
cdef extern from "types.hpp" namespace "fwdpy":
    cdef cppclass singlepop_t:
        singlepop_t(unsigned)
    cdef cppclass GSLrng_t:
        GSLrng_t(unsigned)

##Now, wrap the functions
cdef extern from "neutral.hpp" namespace "fwdpy":

  void evolve_pop(GSLrng_t * rng, singlepop_t * pop, const unsigned & ngens, const double & theta, const double & rho)
#  vector[int] sfs_from_sample(GSLrng_t * rng,const singlepop_t * pop,const unsigned & nsam)

cdef extern from "sample.hpp" namespace "fwdpy":
    vector[pair[double,string]] take_sample_from_pop(GSLrng_t * rng,const singlepop_t * pop,const unsigned & nsam)
    double tajd( const vector[pair[double,string]] & __data )
  
##Creat the python classes
cdef class Singlepop:
    """
    Single-deme object

    The constructor takes a single argument, N, which is the initial population number.
    """
    cdef singlepop_t *thisptr
    def __cinit__(self,unsigned N):
        self.thisptr = new singlepop_t(N)
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


def ms_sample(GSLrng rng, Singlepop pop, int nsam):
    """
    Return a sample of size nsam from a population.

    :param rng: a random-number generator of type GSLrng
    :param pop: an object of type Singlepop
    :param nsam: the desired sample size (should be << than the size of pop)   
    """
    return take_sample_from_pop(rng.thisptr,pop.thisptr,nsam)

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
   
def evolve(GSLrng rng,int N,int ngens,double theta, double rho):
    """
    Evolve a single population

    :param rng: a GSLrng from this module
    :param ngens: number of generations to simulate
    :param theta: :math:`\\theta = 4N_e\\mu` is the scaled mutation rate to variants not affecting fitness ("neutral mutations")
    :param rho: :math:`\\rho = 4N_er` is the scaled recombination rate

    Example:

    >>> import fwdpy
    >>> rng = fwdpy.GSLrng(100)
    >>> pop = fwdpy.evolve(rng,1000,1000,50,50)

    """
    pop = Singlepop(N)
    evolve_pop(rng.thisptr,pop.thisptr,ngens,theta,rho)
    return pop

def evolve_t_details(GSLrng rng,Singlepop p,int ngens,double theta, double rho):
    evolve_pop(rng.thisptr,p.thisptr,ngens,theta,rho)

def evolve_t(GSLrng rng,int nthreads,int N,int ngens,double theta, double rho):
    plist = list()
    threads = []
    for i in range(0,nthreads):
        plist.append(Singlepop(N))
        threads.append(threading.Thread(target=evolve_t_details,args=(rng,plist[i],ngens,theta,rho)))
    for i in range(0,nthreads):
        threads[i].start()
    return plist
    
