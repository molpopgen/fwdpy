from libcpp.string cimport string as cppstring
# distutils: language = c++
cdef class TemporalSampler:
    cpdef size_t size(self):
        return self.vec.size()


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
    def get(self):
        """
        Retrieve the data from the sampler.
        """
        cdef vector[vector[qtrait_stats_cython]] rv
        cdef size_t i=0
        for i in range(self.vec.size()):
            rv.push_back((<pop_properties*>(self.vec[i].get())).final())
        return rv

cdef class PopSampler(TemporalSampler):
    """
    A :class:`fwdpy.fwdpy.TemporalSampler` that takes a sample of size :math:`n \leq N` from the population.
    """
    def __cinit__(self, unsigned n, unsigned nsam,GSLrng
            rng,removeFixed=True,neutral_file=None,selected_file=None,boundaries=None,append=False):
        """
        Constructor
        
        :param n: A length.  Must correspond to number of simulations that will be run simultaneously.
        :param nsam: The sample size to take
        :param rng: A :class:`fwdpy.fwdpy.GSLrng`
        :param removeFixed: (False) Whether or not to include fixations in the output.
        :param neutral_file: (None) File name prefix where neutral data will be written in "ms" format.
        :param selected_file: (None) File name prefix where selected data will be written in "ms" format.
        :param boundaries: (None) For a multi-locus simulation, this must be a list of tuples specifying the positional boundaries of each locus
        :param append: (False) Whether or not to append to output files, or over-write them.

        ..note:: For each of the i threads, the ouput file names will be selected_file.i.gz, etc.
        """
        cdef cppstring sfile,nfile
        cdef vector[pair[double,double]] locus_boundaries
        if boundaries is not None:
            locus_boundaries=boundaries
        for i in range(n):
            sfile.clear()
            nfile.clear()
            if selected_file is not None:
                temp=selected_file+'.'+str(i)+'.gz'
                sfile=temp
            if neutral_file is not None:
                temp=neutral_file+'.'+str(i)+'.gz'
                nfile=temp
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[sample_n](new
                sample_n(nsam,rng.thisptr.get(),nfile,sfile,removeFixed,locus_boundaries,append)))
    def get(self):
        rv=[]
        cdef size_t i=0
        for i in range(self.vec.size()):
            rv.append((<sample_n*>(self.vec[i].get())).final())
        return rv

cdef class VASampler(TemporalSampler):
    """
    A :class:`fwdpy.fwdpy.TemporalSampler` that estimates the relationship between mutation frequency and total additive
    genetic variance.

    .. note:: This is not useful for the standard fwdpy population.  It only actually records anything meaningful in the qtrait and qtrait_mloc modules.  This will change in a future release.
    """
    def __cinit__(self,unsigned n):
        """
        Constructor
        
        :param n: A length.  Must correspond to number of simulations that will be run simultaneously.
        """
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[additive_variance](new additive_variance()))
    def get(self):
        """
        Retrieve the data from the sampler.
        """
        cdef vector[vector[VAcum]] rv
        cdef size_t i=0
        for i in range(self.vec.size()):
            rv.push_back((<additive_variance*>(self.vec[i].get())).final())
        return rv

cdef class freqTrajectories:
    def __cinit__(self):
        self.thisptr=freqTraj(NULL)
    cdef assign(self,freqTraj t):
        self.thisptr = t

cdef class FreqSampler(TemporalSampler):
    """
    A :class:`fwdpy.fwdpy.TemporalSampler` to track the frequencies of selected mutations over time.
    """
    def __cinit__(self,unsigned n):
        """
        Constructor
        
        :param n: A length.  Must correspond to number of simulations that will be run simultaneously.
        """
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[selected_mut_tracker](new selected_mut_tracker()))
    def get(self,rep=None):
        """
        Retrieve the data from the sampler.

        :param rep: If None (the default), then data are returned for all replicates.  Otherwise, rep is the index
        of a replicate, and that replicate's data are returned.

        :raises: RuntimeError if rep is out of range.

        ..note:: This sampler can be *very* RAM-intensive.  For big simulations, we recommend using the 'rep' parameter
        instead of returning everything in a big list.
        """
        cdef freqTraj temp
        if rep is not None:
            if int(rep) > self.vec.size() or int(rep)<0:
                raise RuntimeError("index out of range")
            temp=(<selected_mut_tracker*>self.vec[rep].get()).final()
            t = freqTrajectories() 
            t.assign(temp)
            return t
        else:
            i=0
            rv=[]
            for i in range(self.vec.size()):
                temp=(<selected_mut_tracker*>self.vec[i].get()).final()
                t = freqTrajectories()
                t.assign(temp)
                rv.append(t)
            return rv

def apply_sampler(PopVec pops,TemporalSampler sampler):
    """
    Apply a temporal sampler to a container of populations.

    :param pops: A :class:`fwdpy.fwdpy.PopVec`
    :param sampler: A :class:`fwdpy.fwdpy.TemporalSampler`

    :return: Nothing
    """
    if isinstance(pops,SpopVec):
        apply_sampler_cpp[singlepop_t]((<SpopVec>pops).pops,sampler.vec)
    elif isinstance(pops,MetaPopVec):
        apply_sampler_cpp[metapop_t]((<MetaPopVec>pops).mpops,sampler.vec)
    elif isinstance(pops,MlocusPopVec):
        apply_sampler_cpp[multilocus_t]((<MlocusPopVec>pops).pops,sampler.vec)
    else:
        raise RuntimeError("PopVec type not supported")
        
