# distutils: language = c++
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
    def __cinit__(self, unsigned n, unsigned nsam,GSLrng rng):
        """
        Constructor
        
        :param n: A length.  Must correspond to number of simulations that will be run simultaneously.
        :param nsam: The sample size to take
        :param rng: A :class:`fwdpy.fwdpy.GSLrng`
        """
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[sample_n](new sample_n(nsam,rng.thisptr.get())))
    def get(self):
        cdef vector[vector[pair[uint,detailed_deme_sample]]] rv
        cdef size_t i=0
        for i in range(self.vec.size()):
            rv.push_back((<sample_n*>(self.vec[i].get())).final())
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
    def get(self,unsigned minsojourn = 0, double minfreq = 0.0):
        """
        Retrieve the data from the sampler.
        """
        cdef vector[vector[pair[selected_mut_data, vector[pair[uint,double]]]]] rv
        cdef size_t i=0
        
        for i in range(self.vec.size()):
            rv.push_back((<selected_mut_tracker*>self.vec[i].get()).final())

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
        
