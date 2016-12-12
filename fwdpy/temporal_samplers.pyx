from libcpp.string cimport string as cppstring
from cython.operator cimport dereference as deref

# distutils: language = c++
cdef class TemporalSampler:
    cpdef size_t size(self):
        return self.vec.size()
    def __dealloc__(self):
        clear_samplers(self.vec)
    def force_clear(self):
        clear_samplers(self.vec)

cdef class NothingSampler(TemporalSampler):
    """
    A :class:`fwdpy.fwdpy.TemporalSampler` that does nothing.

    This is useful during evolution to equilibrium or in other situations
    where you are not interested in recording.
    """
    def __cinit__(self, unsigned n):
        """
        Constructor

        :param n: A length.  Must correspond to number of simulations that will be run simultaneously.
        """
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[no_sampling](new no_sampling()))
    def get(self):
        """
        Retrieve the data from the sampler.

        :return: None
        """
        return None
    
cdef class QtraitStatsSampler(TemporalSampler):
    """
    A :class:`fwdpy.fwdpy.TemporalSampler` that records various statistics about the population.

    This type is a model of an iterable container.  Return values may be either yielded
    or accessed via [i].

    .. note:: This is not useful for the standard fwdpy population.  It only actually records anything meaningful in the qtrait and qtrait_mloc modules.  This will change in a future release.
    """
    def __cinit__(self, unsigned n, double optimum):
        """
        Constructor
        
        :param n: A length.  Must correspond to number of simulations that will be run simultaneously.
        :param optimum: The value of the optimum trait/fitness value.
        """
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[pop_properties](new pop_properties(optimum)))
    def __iter__(self):
        for i in range(self.vec.size()):
            yield (<pop_properties*>(self.vec[i].get())).final()
    def __next__(self):
        return next(self)
    def __len__(self):
        return self.vec.size()
    def __getitem__(self,i):
        if i>=self.vec.size():
            raise IndexError("index out of range")
        return (<pop_properties*>(self.vec[i].get())).final()
    def get(self):
        """
        Retrieve the data from the sampler.

        .. note:: This returns all data as a list.  It is more RAM-friendly to iterate over the object.
        """
        cdef vector[vector[qtrait_stats_cython]] rv
        cdef size_t i=0
        for i in range(self.vec.size()):
            rv.push_back((<pop_properties*>(self.vec[i].get())).final())
        return rv

cdef class PopSampler(TemporalSampler):
    """
    A :class:`fwdpy.fwdpy.TemporalSampler` that takes a sample of size :math:`n \leq N` from the population.

    This type is a model of an iterable container.  Return values may be either yielded
    or accessed via [i].
    """
    def __cinit__(self, unsigned n, unsigned nsam,GSLrng
            rng,removeFixed=True,neutral_file=None,selected_file=None,boundaries=None,append=False,recordSamples=True,recordDetails=True):
        """
        Constructor
        
        :param n: A length.  Must correspond to number of simulations that will be run simultaneously.
        :param nsam: The sample size to take
        :param rng: A :class:`fwdpy.fwdpy.GSLrng`
        :param removeFixed: (False) Whether or not to include fixations in the output.
        :param neutral_file: (None) File name (or file name prefix) where neutral data will be written in "ms" format.
        :param selected_file: (None) File name (or file name prefix) where selected data will be written in "ms" format.
        :param boundaries: (None) For a multi-locus simulation, this must be a list of tuples specifying the positional boundaries of each locus
        :param append: (False) Whether or not to append to output files, or over-write them.

        ..note:: 
        
            When n==1, the output file names will be neutral_file and selected file.  When n > 1,
            the names will be neutral_file.i.gz and selected_file.i.gz for all :math:`0\leq i \le n`

        """
        cdef cppstring sfile,nfile
        cdef vector[pair[double,double]] locus_boundaries
        if boundaries is not None:
            locus_boundaries=boundaries
        if n==1:
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[sample_n](new
                sample_n(nsam,rng.thisptr.get(),neutral_file,selected_file,removeFixed,recordSamples,recordDetails,locus_boundaries,append)))
        else:
            for i in range(n):
                sfile.clear()
                nfile.clear()
                if selected_file is not None:
                    temp=selected_file.encode('utf-8')+b'.'+str(i).encode('utf-8')+b'.gz'
                    sfile=temp
                if neutral_file is not None:
                    temp=neutral_file.encode('utf-8')+b'.'+str(i).encode('utf-8')+b'.gz'
                    nfile=temp
                self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[sample_n](new
                sample_n(nsam,rng.thisptr.get(),nfile,sfile,removeFixed,recordSamples,recordDetails,locus_boundaries,append)))
    def __iter__(self):
        for i in range(self.vec.size()):
            yield (<sample_n*>self.vec[i].get()).final()
    def __next__(self):
        return next(self)
    def __getitem__(self, int i):
        if i>= self.vec.size():
            raise IndexError("index out of range")
        return (<sample_n*>self.vec[i].get()).final()
    def __len__(self):
        return self.vec.size()

cdef class VASampler(TemporalSampler):
    """
    A :class:`fwdpy.fwdpy.TemporalSampler` that estimates the relationship between mutation frequency and total additive
    genetic variance.

    This type is a model of an iterable container.  Return values may be either yielded
    or accessed via [i].

    .. note:: This is not useful for the standard fwdpy population.  It only actually records anything meaningful in the qtrait and qtrait_mloc modules.  This will change in a future release.
    """
    def __cinit__(self,unsigned n):
        """
        Constructor
        
        :param n: A length.  Must correspond to number of simulations that will be run simultaneously.
        """
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[additive_variance](new additive_variance()))
    def __iter__(self):
        for i in range(self.vec.size()):
            yield (<additive_variance*>(self.vec[i].get())).final()
    def __next__(self):
        return next(self)
    def __len__(self):
        return self.vec.size()
    def __getitem__(self,i):
        if i>=self.vec.size():
            raise IndexError("index out of range")
        return (<additive_variance*>(self.vec[i].get())).final()
    def get(self):
        """
        Retrieve the data from the sampler.

        .. note:: This returns all data as a list.  It is more RAM-friendly to iterate over the object.
        """
        cdef vector[vector[VAcum]] rv
        cdef size_t i=0
        for i in range(self.vec.size()):
            rv.push_back((<additive_variance*>(self.vec[i].get())).final())
        return rv

cdef class FreqSampler(TemporalSampler):
    """
    A :class:`fwdpy.fwdpy.TemporalSampler` to track the frequencies of selected mutations over time.

    This type is a model of an iterable container.  Return values may be either yielded
    or accessed via [i].
    """
    def __cinit__(self,unsigned n):
        """
        Constructor
        
        :param n: A length.  Must correspond to number of simulations that will be run simultaneously.
        """
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[selected_mut_tracker](new selected_mut_tracker()))
    def __iter__(self):
        for i in range(self.vec.size()):
            yield (<selected_mut_tracker*>self.vec[i].get()).final()
    def __next__(self):
        return next(self)
    def __getitem__(self,i):
        if i>=self.vec.size():
            raise IndexError("index out of range")
        return (<selected_mut_tracker*>self.vec[i].get()).final()
    def __len__(self):
        return self.vec.size()

def apply_sampler(PopVec pops,TemporalSampler sampler):
    """
    Apply a temporal sampler to a container of populations.

    :param pops: A :class:`fwdpy.fwdpy.PopVec`
    :param sampler: A :class:`fwdpy.fwdpy.TemporalSampler`

    :return: Nothing
    """

    if not isinstance(pops,PopVec):
        raise TypeError("Expecting PopVec.")

    if isinstance(pops,SpopVec):
        apply_sampler_cpp[singlepop_t]((<SpopVec>pops).pops,sampler.vec)
    elif isinstance(pops,MetaPopVec):
        apply_sampler_cpp[metapop_t]((<MetaPopVec>pops).mpops,sampler.vec)
    elif isinstance(pops,MlocusPopVec):
        apply_sampler_cpp[multilocus_t]((<MlocusPopVec>pops).pops,sampler.vec)
    else:
        raise RuntimeError("PopVec/PopType type not supported")

def apply_sampler_single(PopType pop,TemporalSampler sampler):
    """
    Apply a temporal sampler to an indivudal :class:`fwdpy.fwdpy.PopType`

    :param pop: A :class:`fwdpy.fwdpy.PopType`
    :param sampler: A :class:`fwdpy.fwdpy.TemporalSampler`

    The use case for this function is applying very expensive temporal samplers
    at the end of a simulation.  It is assumed that len(sampler)==1.
    """
    if not isinstance(pop,PopType):
        raise TypeError("Expecting PopType.")
    if isinstance(pop,Spop):
        apply_sampler_single_cpp[singlepop_t]((<Spop>pop).pop.get(),sampler.vec)
    elif isinstance(pop,MlocusPop):
        apply_sampler_single_cpp[multilocus_t]((<MlocusPop>pop).pop.get(),sampler.vec)
    elif isinstance(pop,MetaPop):
        apply_sampler_single_cpp[metapop_t]((<MetaPop>pop).mpop.get(),sampler.vec)
    else:
        raise NotImplementedError("Not implemented for this type")
