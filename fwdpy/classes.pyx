##Create the python classes
from cython.operator import dereference as deref

cdef class singlepop(poptype):
    """
    Object representing data structures for single-deme simulations.

    Users are not expected to construct these on their own.  Rather,
    they should be working with :class:`popvec`.  This type exists as
    the output of iterating through a :class:`popvec`.
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

cdef class singlepop_gm_vec(poptype):
    def __del__(self):
        self.pop.reset()
    cpdef gen(self):
        return self.pop.get().generation
    cpdef popsize(self):
        return self.pop.get().N
    cpdef sane(self):
        return self.pop.get().sane();

cdef class popvec(popcont):
    """
    Vector of single-deme objects

    Internally, the class contains both a C++ vector of populations and a list of populations.  These two containers
    have pointers to the same objects.  This organization adds little overhead and makes a popvec iterable in the "usual"
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
        cdef size_t size_ = len(self.pypops)
        if self.pops.size() != size_:
            raise RuntimeError("fwdpy.popvec internal data structures out of sync")
        return self.pops.size()
    cdef reset(self,const vector[shared_ptr[singlepop_t]] & newpops):
        self.pops=newpops
        self.pypops=list()
        for i in range(self.pops.size()):
            pi = singlepop();
            pi.pop=self.pops[i]
            self.pypops.append(pi)
    def __append_details__(self,popvec p):
        for i in range(len(p)):
            self.pops.push_back(p.pops[i])
            self.pypops.append(p[i])        
    cpdef append(self,popvec p):
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

cdef class popvec_gmv(popcont):
    def __cinit__(self,unsigned npops,unsigned N):
        """
        Constructor:

        :param npops: The number of populations
        :param N: Initial population number for each population
        """
        self.pypops=list()
        for i in range(npops):
            self.pops.push_back(shared_ptr[singlepop_gm_vec_t](new singlepop_gm_vec_t(N)))
            pi = singlepop_gm_vec()
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
            pi = singlepop_gm_vec();
            pi.pop=self.pops[i]
            self.pypops.append(pi)
    cpdef size(self):
        """
        Returns number of populations (size of underlying C++ vector)
        """
        return self.pops.size()
    
cdef class metapop(poptype):
    """
    Object representing data structures for single-deme simulations.

    Users are not expected to construct these on their own.  Rather,
    they should be working with :class:`mpopvec`.  This type exists as
    the output of iterating through a :class:`mpopvec`.
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
    
cdef class mpopvec(popcont):
    """
    Vector of metapopulation objects
    """
    def __cinit__(self,unsigned nmpops,list Ns):
        """
        Constructor:

        :param nmpops: Number of metapopulations
        :param Ns: A list of population sizes.  The length of this list is the number of demes in each metapopulation
        """
        self.pympops=[]
        for i in range(len(Ns)):
            if Ns[i] < 0:
                raise ValueError("mpopvec: deme size < 0 encountered")
        for i in range(nmpops):
            self.mpops.push_back(shared_ptr[metapop_t](new metapop_t(Ns)))
            pi = metapop()
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
            raise RuntimeError("fwdpy.mpopvec internal data structures out of sync")
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
            pi = metapop()
            pi.mpop = self.mpops[i]
            self.pympops.append(pi)
    def __append_details__(self,mpopvec p):
        for i in range(len(p)):
            self.mpops.push_back(p.mpops[i])
            self.pympops.append(p[i])        
    cpdef append(self,mpopvec p):
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
