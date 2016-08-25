# distutils: language = c++
# distutils: sources = fwdpy/fwdpyio/serialize.cc

from fwdpy.zlib cimport *
from libcpp.string cimport string 
from libc.stdint cimport int64_t
from libcpp.vector cimport vector

ctypedef pair[uint,int64_t] data_t 
ctypedef vector[data_t] gzfinal_t

ctypedef custom_sampler_data[gzfinal_t,string] gzserializer_t

cdef void gzwrite_singlepop(const singlepop_t * pop, const unsigned generation, gzfinal_t & data, string & s) nogil:
    rv=pop.tofile(s.c_str(),True)
    data.push_back(data_t(generation,rv))

cdef void gzwrite_metapop(const metapop_t * pop, const unsigned generation, gzfinal_t & data, string & s) nogil:
    rv=pop.tofile(s.c_str(),True)
    data.push_back(data_t(generation,rv))

cdef void gzwrite_multilocus(const multilocus_t * pop, const unsigned generation, gzfinal_t & data, string & s) nogil:
    rv=pop.tofile(s.c_str(),True)
    data.push_back(data_t(generation,rv))

cdef class gzSerializer(TemporalSampler):
    """
    This class is a :class:`fwdpy.fwdpy.TemporalSampler`, allowing the state of the 
    population to be written to a gzipped file at regular time points.

    ..note:: This is a good way to fill up a hard drive.  Use with caution.
    """
    def __cinit__(self,unsigned n,string basename):
        """
        Constructor.

        :param n: A length.  Must correspond to number of simulations that will be run simultaneously.
        :param basename: A prefix for file names.  For a length n, output file names will be basename.i.gz where 0<=i<n.
        """
        bn=basename
        cdef string temp_string
        cdef gzFile gz
        for i in range(n):
            bni=str(basename)+'.'+str(i)+'.gz'
            temp_string = bni
            gz=gzopen(temp_string.c_str(),"wb")
            gzclose(gz)
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[gzserializer_t](new
                gzserializer_t(&gzwrite_singlepop,temp_string)))
            (<gzserializer_t*>self.vec[i].get()).register_callback(&gzwrite_metapop)
            (<gzserializer_t*>self.vec[i].get()).register_callback(&gzwrite_multilocus)
    def get(self):
        """
        Returns a list for each replicate.  For each replicate, the list consists of the generation and the location of
        that generation in the output file.  The latter value may be used as an offset for later reading the population
        back in.
        """
        rv=[]
        for i in range(self.vec.size()):
            temp=(<gzserializer_t*>self.vec[i].get()).final()
            rvi=[]
            offset=0
            for j in temp:
                rvi.append((j.first,offset))
                offset+=j.second
            rv.append(rvi)
        return rv

##Undocumented fxns are implementation details
def serialize_single(Spop pop):
    return serialize_singlepop(pop.pop.get())

def serialize_meta(MetaPop mpop):
    return serialize_metapop(mpop.mpop.get())

def serialize_mlocus(MlocusPop pop):
    return serialize_multilocus(pop.pop.get())

def serialize(PopType pop):
    """
    Return a binary representation of an evolved population

    :param pop: A list of :class:`fwdpy.fwdpy.PopType`
    
    Example:

    >>> import fwdpy
    >>> import fwdpy.fwdpyio as fpio
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> popsizes=np.tile(popsizes,100)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> strings = [fpio.serialize(i) for i in pops]
    """
    if isinstance(pop,Spop):
        return serialize_single(pop)
    elif isinstance(pop,MetaPop):
        return serialize_meta(pop)
    elif isinstance(pop,MlocusPop):
        return serialize_mlocus(pop)
    else:
        raise RuntimeError("fwdpyio.serialize: unsupported PopType "+str(type(pop)))

def deserialize_singlepops(list strings):
    """
    Convert binary representation back to a :class:`fwdpy.fwdpy.PopVec`

    :param strings: A list of populations in binary format.  This should be the value returned by :func:`fwdpy.fwdpyio.fwdpyio.serialize`

    :returns: :func:`fwdpy.fwdpy.PopVec`

    .. note:: len(strings) determines the length of the return value, and therefore the number of threads to use if the population is evolved further.
        
    Example:
    
    >>> import fwdpy
    >>> import fwdpy.fwdpyio as fpio
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> popsizes=np.tile(popsizes,100)
    >>> pops = fwdpy.evolve_regions(rng,4,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> strings = [fpio.serialize(i) for i in pops]
    >>> len(strings)
    4
    >>> pops2 = fpio.deserialize_singlepops(strings)
    """
    cdef vector[shared_ptr[singlepop_t]] temp = deserialize_singlepop(strings)
    pops=SpopVec(0,0)
    pops.reset(temp)
    return pops

def deserialize_metapops(list strings):
    """
    Convert binary representation of populations back to a :class:`fwdpy.fwdpy.MetaPopVec`

    :param strings: A list of populations in binary format.  This should be the value returned by :func:`fwdpy.fwdpyio.fwdpyio.serialize`

    :returns: :func:`fwdpy.fwdpy.MetaPopVec`

    .. note:: len(strings) determines the length of the return value, and therefore the number of threads to use if the population is evolved further.

    Example:

    TODO
    """
    cdef vector[shared_ptr[metapop_t]] temp = deserialize_metapop(strings)
    mpops = MetaPopVec(0,[0]*1)
    mpops.reset(temp)
    return mpops

def deserialize_mlocus(list strings):
    """
    Convert binary representation of populations back to a :class:`fwdpy.fwdpy.MlocusPopVec`

    :param strings: A list of populations in binary format.  This should be the value returned by :func:`fwdpy.fwdpyio.fwdpyio.serialize`

    :returns: :func:`fwdpy.fwdpy.MlocusPopVec`

    .. note:: len(strings) determines the length of the return value, and therefore the number of threads to use if the population is evolved further.

    Example:

    TODO
    """
    cdef vector[shared_ptr[multilocus_t]] temp = deserialize_multilocus(strings)
    rv = MlocusPopVec(0,0,0)
    rv.reset(temp)
    return rv

