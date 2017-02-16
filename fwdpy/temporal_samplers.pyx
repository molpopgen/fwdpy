from libcpp.string cimport string as cppstring
from cython.operator cimport dereference as deref
import pandas

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

cdef class TrajFilter:
    """
    Base class for filtering trajectories.

    An instance of this object allows all trajectories to be written to file.

    .. note:: See :py:meth:`~fwdpy.fwdpy.FreqSampler.to_sql`
    """
    def __cinit__(self):
        self.tf.reset(new trajFilter())

cdef bool traj_existed_past(const vector[pair[uint,double]] & t,const unsigned & g) nogil:
    if t.empty():
        return False
    if t.back().first >= g:
        return True
    return False

cdef class TrajExistedPast(TrajFilter):
    """
    A type of :class:`fwdpy.fwdpy.TrajFilter`.  Tests that a mutation existed
    in the simulation past a certain generation.
    """
    def __cinit__(self,unsigned g):
        """
        :param g: Generation.  A trajectory is kept if it existed until generation >= g.
        """
        self.tf.reset(new trajFilterData[unsigned](g))
        (<trajFilterData[unsigned]*>self.tf.get()).register_callback(&traj_existed_past)

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
    def __convert_data__(self,dict raw,origin_filter=None,pos_esize_filter=None,freq_filter=None):
        temp=[] #list of dicts with named stuff for pands
        for origin in raw:
            if origin_filter is None or origin_filter(origin) is True:
                for ps in raw[origin]:
                    if pos_esize_filter is None or pos_esize_filter(ps) is True:
                        if freq_filter is None or freq_filter(raw[origin][ps]) is True:
                            temp.extend([{'origin':origin,'pos':ps[0],'esize':ps[1],'generation':i[0],'freq':i[1]} for i in raw[origin][ps]])
        rv=pandas.DataFrame(temp)
        rv.sort_values(by=['origin'])
        rv.drop_duplicates(inplace=True)
        rv.reset_index(inplace=True,drop=True)
        return rv
    def __iter__(self):
        for i in range(self.vec.size()):
            yield self.__convert_data__((<selected_mut_tracker*>self.vec[i].get()).final())
    def __next__(self):
        return next(self)
    def __getitem__(self,i):
        if i>=self.vec.size():
            raise IndexError("index out of range")
        return self.__convert_data__((<selected_mut_tracker*>self.vec[i].get()).final())
    def __len__(self):
        return self.vec.size()
    def fetch(self,i,origin_filter=None,pos_esize_filter=None,freq_filter=None):
        """
        Fetch a filtered data set on allele frequency trajectories.

        :param i: The index of the data set to retrieve.
        :param origin_filter: (None) A callable to filter on the origin time of mutations.
        :param pos_esize_filter: (None) A callable to filter on the position and effect size of a mutation.
        :param freq_filter: (None) A callable to filter on the frequency trajectory itself.

        Each callable must return True or False.  The callable "origin_filter" will receive a single,
        non-negative integer for an argument.  "pos_esize_filter" will recieve a tuple with two elements,
        position and effect size, respectively. Finally, "freq_filter" will recieve a list of tuples.  Each
        tuple will contain (generation, freq) and will be sorted by generation in ascending order.
        """
        if i >= self.vec.size():
            raise IndexError("index out of range")
        raw=(<selected_mut_tracker*>self.vec[i].get()).final()
        return self.__convert_data__(raw,origin_filter,pos_esize_filter,freq_filter)
    def to_sql(self,dbname,TrajFilter traj_filter=None,threshold=1000000,label=0,onedb=False,append=False):
        """
        Write output directly to SQLite database files.  Unlike
        :py:meth:`~fwdpy.fwdpy.FreqSampler.fetch`, this function
        requires an object of type :class:`fwdpy.fwdpy.TrajFilter`
        to filter out unwanted trajectories.

        This function skips the copy of data from C++ to Python, which may make it 
        more efficient than using :py:meth:`~fwdpy.fwdpy.FreqSampler.fetch`.

        :param dbname: Either the name of a database file (when onedb is True), or the prefix for file names (when onedb is False).
        :param traj_filter: (None)  If None, :class:`fwdpy.fwdpy.TrajFilter` is used, which means all trajectories are written to file.  Otherwise, a custom object is used to filter.
        :param threshold: (1,000,000) When onedb is True, this is the number of records to write to the in-memory database before writing to file.
        :param label: (0) The starting value of the replicate id. When onedb is True, data from different replicates will have a "rep" column in the database, with rep goring from threshold to threshold + len(self)-1.
        :param onedb: (False)  If False, each trajectory is written to a separate file.  If true, data are first written to in-memory databases and then flushed to disk at time intervals depending on the value of "threshold".
        :param append: (False) If false, the output file will be deleted if it exsists.  Otherwise, it will be assumed to be a valid SQLite database and appended to.

        .. note:: The schema of the resulting files can be checked with the sqlite3 command-line tool.
        """
        if traj_filter is None:
            traj_filter=TrajFilter()
        cdef shared_ptr[mutex] dblock
        dblock.reset(new mutex())
        traj2sql(self.vec,dblock,
                traj_filter.tf.get(),
                dbname,threshold,label,onedb,append)


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
