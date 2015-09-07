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
    def get(self,unsigned i):
        return self.pypops[i]
    def __len__(self):
        return self.thisptr.size()
    def size(self):
        """
        Returns number of populations (size of underlying C++ vector)
        """
        return self.pops.size()
    def generation(self,unsigned i):
        """
        Returns the generation that population 'i' is currently evolved to

        :param i: index of the population for which to return the generation
        """
        cdef const singlepop_t * pp = self.pops[i].get()
        return pp.gen()
    def popsize(self,unsigned i):
        """
        Returns the size of population 'i'

        :param i: index of the population for which to return the population size
        """
        cdef const singlepop_t * pp = self.pops[i].get()
        return pp.popsize()
    def sane(self,unsigned i):
        """
        Makes sure that population 'i' is in a sane state.

        Internally, this checks that pop[i]->N == pop[i]->diploids.size(),
        which it should be if the C++ code behind this all is properly updating
        the data structures!

        :param i: index of the population to check
        """
        cdef const singlepop_t * pp = self.pops[i].get()
        return pp.sane()

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
