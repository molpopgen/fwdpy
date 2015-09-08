##Create the python classes

cdef class singlepop:
    """
    Object representing data structures for single-deme simulations.

    Users are not expected to construct these on their own.  Rather,
    they should be working with "popvec".  This type exists as
    the output of iterating through a "popvec".
    """
    cdef shared_ptr[singlepop_t] pop
    def __del__(self):
       self.pop.reset()
    def popsize(self):
       cdef singlepop_t * pp = self.pop.get()
       return pp.popsize()
    def gen(self):
        """
        Returns the generation that the population is currently evolved to
        """
        return self.pop.get().generation
    def popsize(self):
        """
        Returns the size of the population
        """
        return self.pop.get().N
    def sane(self):
        """
        Makes sure that the population is in a sane state.

        Internally, this checks that pop.N == pop.diploids.size(),
        which it should be if the C++ code behind this all is properly updating
        the data structures!

        """
        return self.pop.get().sane()
    
cdef class popvec:
    """
    Vector of single-deme objects

    The constructor takes two objects: the number of pops, and the initial population size (which is the same for each pop).

    Internally, the class contains both a C++ vector of populations and a list of populations.  These two containers
    have pointers to the same objects.  This organization adds little overhead and makes a popvec iterable in the "usual"
    Python way.
    """
    cdef vector[shared_ptr[singlepop_t]] pops
    pypops = list()
    def __cinit__(self,unsigned npops,unsigned N):
        for i in range(npops):
            self.pops.push_back(shared_ptr[singlepop_t](new singlepop_t(N)))
            pi = singlepop()
            pi.pop = self.pops[i]
            self.pypops.append(pi)
    def __iter__(self):
        return iter(self.pypops)
    def __next__(self):
        return next(self.pypops)
    def __getitem__(self, int i):
        return self.pypops[i]
    def __len__(self):
        return self.pops.size()
    def size(self):
        """
        Returns number of populations (size of underlying C++ vector)
        """
        return self.pops.size()

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

##These are the callback wrappers from fwdpp
cdef class shwrapper:
    cdef shmodel thisptr
    def __cinit__(self):
        self.thisptr = shmodel();

# cdef class constsh:
#     cdef constant * thisptr
#     def __cinit__(self,double val):
#         self.thisptr = new constant(val)
#     def __dealloc__(self):
#         del self.thisptr

# cdef class expsh:
#     cdef exponential * thisptr
#     def __cinit__(self,double mean):
#         self.thisptr = new exponential(mean)
#     def __dealloc__(self):
#         del self.thisptr

# cdef class betash:
#     cdef beta * thisptr
#     def __cinit__(self,double a, double b, double factor = 1):
#         self.thisptr = new beta(a,b,factor)
#     def __dealloc__(self):
#         del self.thisptr

# cdef class gaussiansh:
#     cdef gaussian * thisptr
#     def __cinit__(self,double sd):
#         self.thisptr = new gaussian(sd)
#     def __dealloc__(self):
#         del self.thisptr

# cdef class gammash:
#     cdef gamma * thisptr
#     def __cinit__(self,double mean,double shape):
#         self.thisptr = new gamma(mean,shape)
#     def __dealloc__(self):
#         del self.thisptr
