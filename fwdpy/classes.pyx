##Create the python classes
from cython.operator import dereference as deref

cdef class Spop(PopType):
    """
    Object representing data structures for single-deme simulations
    based on a mutation type having a single 's' and 'h' term.

    Users are not expected to construct these on their own.  Rather,
    they should be working with :class:`SpopVec`.  This type exists as
    the output of iterating through a :class:`SpopVec`.
    """
    def __dealloc__(self):
       self.pop.reset()
    cpdef gen(self):
        """
        Returns the generation that the population is currently evolved to
        """
        return self.pop.get().generation
    cpdef popsize(self):
        """
        Returns the size of the population
        """
        return self.pop.get().N
    cpdef sane(self):
        """
        Makes sure that the population is in a sane state.

        Internally, this checks that pop.N == pop.diploids.size(),
        which it should be if the C++ code behind this all is properly updating
        the data structures!

        """
        return self.pop.get().sane()

cdef class SpopGenMut(PopType):
    """
    Object representing data structures for single-deme simulations
    based on a mutation type having a vector of properties.

    Users are not expected to construct these on their own.  Rather,
    they should be working with :class:`SpopGenMutVec`.  This type exists as
    the output of iterating through a :class:`SpopGenMut`.

    ..note:: Currently, there are no functions in fwdpy using this type!  See :class:`fwdpy.fwdpy.Spop` instead.
    """
    def __dealloc__(self):
        self.pop.reset()
    cpdef gen(self):
        """
        Returns the generation that the population is currently evolved to
        """
        return self.pop.get().generation
    cpdef popsize(self):
        """
        Returns the size of the population
        """
        return self.pop.get().N
    cpdef sane(self):
        """
        Makes sure that the population is in a sane state.

        Internally, this checks that pop.N == pop.diploids.size(),
        which it should be if the C++ code behind this all is properly updating
        the data structures!

        """
        return self.pop.get().sane()

cdef class MlocusPop(PopType):
    """
    Object representing data structures for single-deme, multi-locus/region simulations
    based on a mutation type having a single 's' and 'h' term.

    Users are not expected to construct these on their own.  Rather,
    they should be working with :class:`MlocusPopVec`.  This type exists as
    the output of iterating through a :class:`MlocusPopVec`.
    """
    def __dealloc__(self):
        self.pop.reset()
    cpdef gen(self):
        """
        Returns the generation that the population is currently evolved to
        """
        return self.pop.get().generation
    cpdef popsize(self):
        """
        Returns the size of the population
        """
        return self.pop.get().N
    cpdef sane(self):
        """
        Makes sure that the population is in a sane state.

        Internally, this checks that pop.N == pop.diploids.size(),
        which it should be if the C++ code behind this all is properly updating
        the data structures!

        """
        return self.pop.get().sane()  

cdef class SpopVec(PopVec):
    """
    Vector of single-deme objects

    Internally, the class contains "smart" pointers to the underlying C++ types.

    The class is iterable, yielding :class:`fwdpy.fwdpy.Spop` objects.

    See :func:`evolve_regions` for use cases.
    """
    def __cinit__(self,unsigned npops,unsigned N):
        """
        Constructor:

        :param npops: The number of populations
        :param N: Initial population number for each population
        """
        for i in range(npops):
            self.pops.push_back(shared_ptr[singlepop_t](new singlepop_t(N)))
    def __dealloc__(self):
        self.pops.clear()
    def __iter__(self):
        for i in range(self.pops.size()):
            x = Spop()
            x.pop=self.pops[i]
            yield x
    def __next__(self):
        return next(self)
    def __getitem__(self, int i):
        x = Spop()
        x.pop=self.pops[i]
        return x
    def __len__(self):
        return self.pops.size()
    cdef reset(self,const vector[shared_ptr[singlepop_t]] & newpops):
        self.pops=newpops
    def __append_details__(self,SpopVec p):
        for i in range(p.pops.size()):
            self.pops.push_back(p.pops[i])
    cpdef append(self,SpopVec p):
        """
        Append 'p' into this object.

        This is done via a serialized copy, meaning that 
        this object and p will not share any pointers
        """
        self.__append_details__(copypops(p))
    cpdef size(self):
        """
        Returns number of populations (size of underlying C++ vector)
        """
        return self.pops.size()
    def clear(self):
        """
        Frees the memory allocated by each simulated population.

        Result is an empty object.
        """
        self.pops.clear()

cdef class MlocusPopVec(PopVec):
    """
    Vector of single-deme objects representing multiple partially-linked regions.

    Internally, the class contains "smart" pointers to the underlying C++ types.

    The class is iterable, yielding :class:`fwdpy.fwdpy.MlocusPop` objects.

    See :func:`evolve_regions` for use cases.
    """
    def __cinit__(self,unsigned npops,unsigned N,unsigned nloci):
        """
        Constructor:

        :param npops: The number of populations
        :param N: Initial population number for each population
        :param nloci: Number of loci/regions
        """
        for i in range(npops):
            self.pops.push_back(shared_ptr[multilocus_t](new multilocus_t(N,nloci)))
    def __dealloc__(self):
        self.pops.clear()
    def __iter__(self):
        for i in range(self.pops.size()):
            pi = MlocusPop()
            pi.pop = self.pops[i]
            yield pi
    def __next__(self):
        return next(self)
    def __getitem__(self, int i):
        pi = MlocusPop()
        pi.pop = self.pops[i]
        return pi
    def __len__(self):
        return self.pops.size()
    cdef reset(self,const vector[shared_ptr[multilocus_t]] & newpops):
        self.pops=newpops
    def __append_details__(self,MlocusPopVec p):
        for i in range(len(p)):
            self.pops.push_back(p.pops[i])
    cpdef append(self,MlocusPopVec p):
        """
        Append 'p' into this object.

        This is done via a serialized copy, meaning that 
        this object and p will not share any pointers
        """
        self.__append_details__(copypops(p))
             
    cpdef size(self):
        """
        Returns number of populations (size of underlying C++ vector)
        """
        return self.pops.size()
    def clear(self):
        """
        Frees the memory allocated by each simulated population.

        Result is an empty object.
        """
        self.pops.clear()
    
cdef class SpopGenMutVec(PopVec):
    def __cinit__(self,unsigned npops,unsigned N):
        """
        Constructor:

        :param npops: The number of populations
        :param N: Initial population number for each population
        """
        for i in range(npops):
            self.pops.push_back(shared_ptr[singlepop_gm_vec_t](new singlepop_gm_vec_t(N)))
    def __dealloc__(self):
        self.pops.clear()
    def __iter__(self):
        for i in range(self.pops.size()):
            pi = SpopGenMut()
            pi.pop = self.pops[i]
            yield pi
    def __next__(self):
        return next(self)
    def __getitem__(self, int i):
        pi = SpopGenMut()
        pi.pop = self.pops[i]
        return pi
    def __len__(self):
        return self.pops.size()
    cdef reset(self,const vector[shared_ptr[singlepop_gm_vec_t]] & newpops):
        self.pops=newpops
    cpdef size(self):
        """
        Returns number of populations (size of underlying C++ vector)

        Result is an empty object.
        """
        return self.pops.size()
    def clear(self):
        """
        Frees the memory allocated by each simulated population.

        Result is an empty object.
        """
        self.pops.clear()
    
cdef class MetaPop(PopType):
    """
    Object representing data structures for single-deme simulations.

    Users are not expected to construct these on their own.  Rather,
    they should be working with :class:`MetaPopVec`.  This type exists as
    the output of iterating through a :class:`MetaPopVec`.
    """
    def __dealloc__(self):
       self.mpop.reset()
    def __len__(self):
        return self.mpop.get().size()
    cpdef gen(self):
        """
        Returns the generation that the population is currently evolved to
        """
        return self.mpop.get().generation
    cpdef popsizes(self):
        """
        Returns the size of the population
        """
        return self.mpop.get().Ns
    cpdef sane(self):
        """
        Makes sure that the population is in a sane state.

        Internally, this checks that pop.N == pop.diploids.size(),
        which it should be if the C++ code behind this all is properly updating
        the data structures!

        """
        return self.mpop.get().sane()
    cpdef from_Spop(self,Spop p):
         self.mpop.reset(new metapop_t(deref(p.pop.get())))
         
cdef class MetaPopVec(PopVec):
    """
    Vector of metapopulation objects.

    Internally, the class contains "smart" pointers to the underlying C++ types.

    The class is iterable, yielding :class:`fwdpy.fwdpy.MetaPop` objects.

    """
    def __cinit__(self,unsigned nmpops,vector[unsigned] Ns):
        """
        Constructor:

        :param nmpops: Number of metapopulations
        :param Ns: A list of population sizes.  The length of this list is the number of demes in each metapopulation
        """
        for i in range(len(Ns)):
            if Ns[i] < 0:
                raise ValueError("MetaPopVec: deme size < 0 encountered")
        for i in range(nmpops):
            self.mpops.push_back(shared_ptr[metapop_t](new metapop_t(Ns)))
    def __dealloc__(self):
        self.mpops.clear()
    def __iter__(self):
        for i in range(self.mpops.size()):
            pi = MetaPop()
            pi.mpop = self.mpops[i]
            yield pi
    def __next__(self):
        return next(self)
    def __getitem__(self, int i):
        pi = MetaPop()
        pi.mpop = self.mpops[i]
        return pi
    def __len__(self):
        return self.mpops.size()
    cpdef size(self):
        """
        Returns number of populations (size of underlying C++ vector)
        """
        return self.mpops.size()
    cdef reset(self,const vector[shared_ptr[metapop_t]] & mpops):
        self.mpops = mpops
    def __append_details__(self,MetaPopVec p):
        for i in range(len(p)):
            self.mpops.push_back(p.mpops[i])
    cpdef append(self,MetaPopVec p):
        """
        Append 'p' into this object.

        This is done via a serialized copy, meaning that 
        this object and p will not share any pointers
        """
        self.__append_details__(copypops(p))
    def clear(self):
        """
        Frees the memory allocated by each simulated population.

        Result is an empty object.
        """
        self.mpops.clear()

cdef class GSLrng:
    """
    A wrapper around a random number generator (rng) 
    from the GNU Scientific Library (GSL).

    The constructor takes a seed (int) as an argument.

    Example:

    >>> import fwdpy
    >>> rng = fwdpy.GSLrng(100)
    """
    def __cinit__(self, int seed):
        """
        Constructor:

        :param seed: The seed for the RNG
        """
        self.thisptr = new GSLrng_t(seed)
    def __dealloc__(self):
        del self.thisptr
