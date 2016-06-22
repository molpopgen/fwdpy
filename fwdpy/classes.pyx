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
    def __del__(self):
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
    def __del__(self):
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
    def __del__(self):
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

    Internally, the class contains both a C++ vector of populations and a list of populations.  These two containers
    have pointers to the same objects.  This organization adds little overhead and makes a SpopVec iterable in the "usual"
    Python way.

    See :func:`evolve_regions` for use cases.
    """
    def __cinit__(self,unsigned npops,unsigned N):
        """
        Constructor:

        :param npops: The number of populations
        :param N: Initial population number for each population
        """
        self.pypops=list()
        for i in range(npops):
            self.pops.push_back(shared_ptr[singlepop_t](new singlepop_t(N)))
            pi = Spop()
            pi.pop = self.pops[i]
            self.pypops.append(pi)
    def __iter__(self):
        return iter(self.pypops)
    def __next__(self):
        return next(self.pypops)
    def __getitem__(self, int i):
        return self.pypops[i]
    def __len__(self):
        cdef size_t size_ = len(self.pypops)
        if self.pops.size() != size_:
            raise RuntimeError("fwdpy.SpopVec internal data structures out of sync")
        return self.pops.size()
    cdef reset(self,const vector[shared_ptr[singlepop_t]] & newpops):
        self.pops=newpops
        self.pypops=list()
        for i in range(self.pops.size()):
            pi = Spop()
            pi.pop=self.pops[i]
            self.pypops.append(pi)
    def __append_details__(self,SpopVec p):
        for i in range(len(p)):
            self.pops.push_back(p.pops[i])
            self.pypops.append(p[i])        
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

cdef class MlocusPopVec(PopVec):
    """
    Vector of single-deme objects representing multiple partially-linked regions.

    Internally, the class contains both a C++ vector of populations and a list of populations.  These two containers
    have pointers to the same objects.  This organization adds little overhead and makes a popvec iterable in the "usual"
    Python way.

    See :func:`evolve_regions` for use cases.
    """
    def __cinit__(self,unsigned npops,unsigned N,unsigned nloci):
        """
        Constructor:

        :param npops: The number of populations
        :param N: Initial population number for each population
        :param nloci: Number of loci/regions
        """
        self.pypops=list()
        for i in range(npops):
            self.pops.push_back(shared_ptr[multilocus_t](new multilocus_t(N,nloci)))
            pi = MlocusPop()
            pi.pop = self.pops[i]
            self.pypops.append(pi)
    def __iter__(self):
        return iter(self.pypops)
    def __next__(self):
        return next(self.pypops)
    def __getitem__(self, int i):
        return self.pypops[i]
    def __len__(self):
        cdef size_t size_ = len(self.pypops)
        if self.pops.size() != size_:
            raise RuntimeError("fwdpy.MlocusPopVec internal data structures out of sync")
        return self.pops.size()
    cdef reset(self,const vector[shared_ptr[multilocus_t]] & newpops):
        self.pops=newpops
        self.pypops=list()
        for i in range(self.pops.size()):
            pi = MlocusPop()
            pi.pop=self.pops[i]
            self.pypops.append(pi)
    def __append_details__(self,MlocusPopVec p):
        for i in range(len(p)):
            self.pops.push_back(p.pops[i])
            self.pypops.append(p[i])        
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
    
cdef class SpopGenMutVec(PopVec):
    def __cinit__(self,unsigned npops,unsigned N):
        """
        Constructor:

        :param npops: The number of populations
        :param N: Initial population number for each population
        """
        self.pypops=list()
        for i in range(npops):
            self.pops.push_back(shared_ptr[singlepop_gm_vec_t](new singlepop_gm_vec_t(N)))
            pi = SpopGenMut()
            pi.pop = self.pops[i]
            self.pypops.append(pi)
    def __iter__(self):
        return iter(self.pypops)
    def __next__(self):
        return next(self.pypops)
    def __getitem__(self, int i):
        return self.pypops[i]
    def __len__(self):
        cdef size_t size_ = len(self.pypops)
        if self.pops.size() != size_:
            raise RuntimeError("fwdpy.popvec internal data structures out of sync")
        return self.pops.size()
    cdef reset(self,const vector[shared_ptr[singlepop_gm_vec_t]] & newpops):
        self.pops=newpops
        self.pypops=list()
        for i in range(self.pops.size()):
            pi = SpopGenMut()
            pi.pop=self.pops[i]
            self.pypops.append(pi)
    cpdef size(self):
        """
        Returns number of populations (size of underlying C++ vector)
        """
        return self.pops.size()
    
cdef class MetaPop(PopType):
    """
    Object representing data structures for single-deme simulations.

    Users are not expected to construct these on their own.  Rather,
    they should be working with :class:`MetaPopVec`.  This type exists as
    the output of iterating through a :class:`MetaPopVec`.
    """
    def __del__(self):
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
    Vector of metapopulation objects
    """
    def __cinit__(self,unsigned nmpops,vector[unsigned] Ns):
        """
        Constructor:

        :param nmpops: Number of metapopulations
        :param Ns: A list of population sizes.  The length of this list is the number of demes in each metapopulation
        """
        self.pympops=[]
        for i in range(len(Ns)):
            if Ns[i] < 0:
                raise ValueError("MetaPopVec: deme size < 0 encountered")
        for i in range(nmpops):
            self.mpops.push_back(shared_ptr[metapop_t](new metapop_t(Ns)))
            pi = MetaPop()
            pi.mpop = self.mpops[i]
            self.pympops.append(pi)
    def __iter__(self):
        return iter(self.pympops)
    def __next__(self):
        return next(self.pympops)
    def __getitem__(self, int i):
        return self.pympops[i]
    def __len__(self):
        cdef size_t size_ = len(self.pympops)
        if self.mpops.size() != size_:
            raise RuntimeError("fwdpy.MetaPopVec internal data structures out of sync")
        return self.mpops.size()
    cpdef size(self):
        """
        Returns number of populations (size of underlying C++ vector)
        """
        return self.mpops.size()
    cdef reset(self,const vector[shared_ptr[metapop_t]] & mpops):
        self.mpops = mpops
        self.pympops = []
        for i in range(self.mpops.size()):
            pi = MetaPop()
            pi.mpop = self.mpops[i]
            self.pympops.append(pi)
    def __append_details__(self,MetaPopVec p):
        for i in range(len(p)):
            self.mpops.push_back(p.mpops[i])
            self.pympops.append(p[i])        
    cpdef append(self,MetaPopVec p):
        """
        Append 'p' into this object.

        This is done via a serialized copy, meaning that 
        this object and p will not share any pointers
        """
        self.__append_details__(copypops(p))
        
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
