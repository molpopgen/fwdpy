##Create the python classes

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
    def __len__(self):
        return self.thisptr.size()
    def size(self):
        """
        Returns number of populations (size of underlying C++ vector)
        """
        return self.thisptr.size()
    def generation(self,unsigned i):
        """
        Returns the generation that population 'i' is currently evolved to

        :param i: index of the population for which to return the generation
        """
        return self.thisptr.generation(i)
    def popsize(self,unsigned i):
        """
        Returns the size of population 'i'

        :param i: index of the population for which to return the population size
        """
        return self.thisptr.popsize(i)
    def sane(self,unsigned i):
        """
        Makes sure that population 'i' is in a sane state.

        Internally, this checks that pop[i]->N == pop[i]->diploids.size(),
        which it should be if the C++ code behind this all is properly updating
        the data structures!

        :param i: index of the population to check
        """
        return self.thisptr.sane(i)
    
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
