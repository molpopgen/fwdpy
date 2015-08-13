# distutils: language = c++
# distutils: sources = src/sample.cpp

from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

##Create hooks to C++ types

#Wrap the classes:
cdef extern from "types.hpp" namespace "fwdpy":
    cdef cppclass singlepop_t:
        singlepop_t(unsigned)
    cdef cppclass GSLrng_t:
        GSLrng_t(unsigned)

##Now, wrap the functions
#cdef extern from "neutral.hpp" namespace "fwdpy":
#  void evolve_pop(GSLrng_t * rng, singlepop_t * pop, const unsigned & ngens, const double & theta, const double & rho)
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

    >>> rng = GSLrng(100)
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
   return tajd(data)
