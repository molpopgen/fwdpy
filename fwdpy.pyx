# distutils: language = c++
# distutils: sources = src/sample.cpp

from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.string cimport string

##Create hooks to C++ types

#Wrap the classes:
cdef extern from "types.hpp" namespace "fwdpy":
    cdef cppclass singlepop_t:
        singlepop_t(unsigned,unsigned)
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
    cdef singlepop_t *thisptr
    def __cinit__(self,int N, int reserve_size = 100):
        self.thisptr = new singlepop_t(N,reserve_size)
    def __dealloc__(self):
        del self.thisptr

cdef class GSLrng:
    cdef GSLrng_t * thisptr
    def __cinit__(self, int seed):
        self.thisptr = new GSLrng_t(seed)
    def __dealloc__(self):
        del self.thisptr


##OK--this works!

#def sfs_sample(GSLrng rng, Singlepop pop, int nsam):
#    return sfs_from_sample(rng.thisptr,pop.thisptr,nsam)

# def evolve(GSLrng rng,int N,int ngens,double theta, double rho):
#     pop = Singlepop(N)
#     evolve_pop(rng.thisptr,pop.thisptr,ngens,theta,rho)
#     return pop

#def ms_sample(GSLrng rng, Singlepop pop, int nsam):
#    return take_sample_from_pop(rng.thisptr,pop.thisptr,nsam)

#def TajimasD( vector[pair[double,string]] data ):
#    return tajd(data)
