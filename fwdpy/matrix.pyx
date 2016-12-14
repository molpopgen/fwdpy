from fwdpy cimport *
from cython.operator cimport dereference as deref
from fwdpy cimport uint
import numpy
import array

cdef key_pair remove_fixed_keys(key_pair & keys,const uint n) nogil:
    cdef key_pair rv
    for i in keys.first:
        if i.second < n:
            rv.first.push_back(i)
    for i in keys.second:
        if i.second < n:
            rv.second.push_back(i)
    return rv

cdef key_pair apply_min_sample_daf(key_pair & keys,const double n,const double x) nogil:
    cdef key_pair rv
    for i in keys.first:
        if <double>i.second/n >= x:
            rv.first.push_back(i)
    for i in keys.second:
        if <double>i.second/n >= x:
            rv.second.push_back(i)
    return rv

cdef key_pair apply_min_pop_daf(key_pair & keys,const unsigned twoN, const double q,const ucont_t & mcounts) nogil:
    cdef key_pair rv
    for i in keys.first:
        if <double>mcounts[i.first]/<double>twoN >= q:
            rv.first.push_back(i)
    for i in keys.second:
        if <double>mcounts[i.first]/<double>twoN >= q:
            rv.second.push_back(i)
    return rv

cdef sort_keys_by_position(key_pair & keys,const mcont_t & mutations):
    #positions of neutral mutations
    pos = [mutations[i[0]].pos for i in keys.first]         
    keys.first=[i for j,i in sorted(zip(pos,keys.first),key=lambda pair: pair[0])]
    pos = [mutations[i[0]].pos for i in keys.second]         
    keys.second=[i for j,i in sorted(zip(pos,keys.second),key=lambda pair: pair[0])]

cdef class DataMatrix(object):
    """
    Matrix representation of a sample.
    All array types are array.array.
    
    .. versionadded:: 0.0.4
    """
    def __cinit__(self,const data_matrix & d):
        self.neutral=array.array('h',d.neutral)
        self.selected=array.array('h',d.selected)
        self.neutral_positions=array.array('d',d.neutral_positions)
        self.selected_positions=array.array('d',d.selected_positions)
        self.neutral_popfreq=array.array('d',d.neutral_popfreq)
        self.selected_popfreq=array.array('d',d.selected_popfreq)
        self.nrow=d.nrow
        self.nn=len(self.neutral_positions)
        self.ns=len(self.selected_positions)
    def nmatrix(self,dtype = numpy.dtype('i1')):
        """
        Returns matrix of neutral sites

        :param dtype: (numpy.dtype('i1')) Data type for matrix.  Default is 8-bit integer.
        :rtype: numpy.array or None if the matrix is empty
        """
        if self.nn==0:
            return None
        return numpy.reshape(numpy.array(self.neutral,dtype=dtype),[self.nrow,self.nn])
    def smatrix(self,dtype = numpy.dtype('i1')):
        """
        Returns matrix of selected sites

        :param dtype: (numpy.dtype('i1')) Data type for matrix.  Default is 8-bit integer.
        :rtype: numpy.array or None if the matrix is empty
        """
        if self.ns==0:
            return None
        return numpy.reshape(numpy.array(self.neutral,dtype=dtype),[self.nrow,self.nn])
    def matrix(self,dtype = numpy.dtype('i1')):
        """
        Returns matrix of neutral and selected sites.  If either matrix is empty, the other
        is returned.  Otherwise, the neutral and selected matrix are appended (neutral columns
        first, then selected).

        :param dtype: (numpy.dtype('i1')) Data type for matrix.  Default is 8-bit integer.
        :rtype: numpy.array or None if the matrix is empty
        """
        n = self.nmatrix(dtype)
        s = self.smatrix(dtype)
        if n is None:
            if s is None: return None
            else: return s
        if s is None: return n
        return numpy.append(n,s,1)

cdef class GenotypeMatrix(DataMatrix):
    """
    Interface is identical to :class:`fwdpy.fwdpy.DataMatrix`
    """
    pass

cdef class HaplotypeMatrix(DataMatrix):
    """
    An extension of :class:`fwdpy.fwdpy.DataMatrix`
    that allows for conversion to data formats compatible
    wiht pylibseq/libsequence.
    """
    def as_sample(self):
        """
        Return a representation of the data in a form 
        usable in pylibseq.

        :return: A dictionary with elements 'neutral' and 'selected'

        :rtype: dict

        """
        nm = self.nmatrix()
        sm = self.smatrix()
        #These will be lists of lists of characters,
        #to be converted into list of tuple of float,string
        ss=[]

        rv={'neutral':[],'selected':[]}
        if nm is not None:
            ns=[]
            for col in range(nm.shape[1]):
                nsi=numpy.array(['0']*nm.shape[0])
                nsi[numpy.where(nm[:,col]==1)[0]]='1'
                nsi[numpy.where(nm[:,col]==2)[0]]='2'
                ns.append(nsi)
            rv['neutral']=[(i,"".join(j)) for i,j in zip(self.neutral_positions,ns)]
        if sm is not None:
            ss=[]
            for col in range(sm.shape[1]):
                ssi=numpy.array(['0']*sm.shape[0])
                ssi[numpy.where(sm[:,col]==1)[0]]='1'
                ssi[numpy.where(sm[:,col]==2)[0]]='2'
                ss.append(ssi)
            rv['selected']=[(i,"".join(j)) for i,j in zip(self.selected_positions,ss)]
        return rv

def get_mutation_keys(pop,list individuals,include_neutral=True,include_selected=True,remove_fixed=False,sort=True,min_sample_daf=None,min_pop_daf=None,deme=None):
    """
    Generate keys to the set of mutations present in the sample defined by a set of individuals.
    The keys may be used to generate a :class:`fwdpy.matrix.DataMatrix` via a call to 
    :func:`fwdpy.matrix.haplotype_matrix` and/or :func:`fwdpy.matrix.genotype_matrix`.

    :param pop: A :class:`fwdpy.fwdpy.PopType`
    :param individuals: A list of indexes of individuals in the sample
    :param include_neutral: (True) If True, generate keys for neutral mutations
    :param include_selected: (True) If True, generate keys for selected mutations
    :param remove_fixed: (False) If True, prune the return value of keys to mutations that are fixed in the sample
    :param sort: (True) If True, keys are sorted according to mutation positions
    :param min_sample_daf: (None) If not None, then it must be a float (0 <= min_sample_daf < 1) and only mutations with derived allele frequency **in the sample** :math:`q \geq  \mathrm{min_sample_daf}` will be kept.
    :param min_pop_daf: (None) If not None, then it must be a float (0 <= min_pop_daf < 1) and only mutations with derived allele frequency **in the (meta-)population** :math:`q \geq  \mathrm{min_pop_daf}` will be kept.
    :param deme: If pop is a metapopulation, this is the index of the deme.

    :return: Keys to neutral and selected mutations stored in separate lists.
    
    :rtype: tuple of two lists

    .. versionadded:: 0.0.4

    .. note::
        
        For meta-populations, min_pop_daf applies to the entire meta-population, and not the specific deme.
        Other types of key sorting must be done by the user.  This is currently difficult without writing
        a custom Cython plugin based on fwdpy.  This will be improved in the future with a new "views" functionality.
    """
    deme_=deme
    if deme is None:
        deme_=0
    if min_sample_daf is not None:
        if float(min_sample_daf) < 0. or float(min_sample_daf) >= 1.:
            raise RuntimeError("Minimum derived allele frequency in sample must be a float 0.0 <= 1 < 1.0")
    if min_pop_daf is not None:
        if float(min_pop_daf) < 0. or float(min_pop_daf) >= 1.:
            raise RuntimeError("Minimum derived allele frequency in sample must be a float 0.0 <= 1 < 1.0")

    #Get the mutation keys first
    cdef key_pair keys
    if isinstance(pop,Spop):
        keys = mutation_keys[singlepop_t](deref((<Spop>pop).pop.get()),individuals,include_neutral,include_selected,<size_t>deme_)
        if sort is True:
            sort_keys_by_position(keys,(<Spop>pop).pop.get().mutations)
        if min_pop_daf is not None:
            keys = apply_min_pop_daf(keys,2*<Spop>pop.popsize(),min_pop_daf,<Spop>pop.pop.get().mcounts)
    elif isinstance(pop,MlocusPop):
        keys = mutation_keys[multilocus_t](deref((<MlocusPop>pop).pop.get()),individuals,include_neutral,include_selected,<size_t>deme_)
        if sort is True:
            sort_keys_by_position(keys,(<MlocusPop>pop).pop.get().mutations)
        if min_pop_daf is not None:
            keys = apply_min_pop_daf(keys,2*<MlocusPop>pop.popsize(),min_pop_daf,<MlocusPop>pop.pop.get().mcounts)
    elif isinstance(pop,MetaPop):
        if deme is None:
            raise RuntimeError("deme cannot be none")
        keys = mutation_keys[metapop_t](deref((<MetaPop>pop).mpop.get()),individuals,include_neutral,include_selected,<size_t>deme_)
        if sort is True:
            sort_keys_by_position(keys,(<MetaPop>pop).mpop.get().mutations)
        if min_pop_daf is not None:
            keys = apply_min_pop_daf(keys,2*sum(<MetaPop>pop.popsizes()),min_pop_daf,<MetaPop>pop.mpop.get().mcounts)
    #Apply filters to keys
    if remove_fixed is True:
        keys = remove_fixed_keys(keys,2*len(individuals))
    if min_sample_daf is not None:
        keys = apply_min_sample_daf(keys,2*len(individuals),min_sample_daf)
    return keys

def haplotype_matrix(pop,list individuals,include_neutral=True,include_selected=True,remove_fixed=False,sort=True,min_sample_daf=None,min_pop_daf=None,deme=None,keys=None):
    """
    Generate a haplotype matrix stored in a :class:`fwdpy.matrix.HaplotypeMatrix`
    
    All arguments are identical to :func:`fwdpy.matrix.get_mutation_keys` with the exception of:

    :param keys: (None) Pre-calculated mutation keys.  Must be a tuple of lists of positive integers. Exceptions will be raised if keys are invalid or out of range.

    :rtype: :class:`fwdpy.matrix.HaplotypeMatrix`
    """
    deme_=deme
    if deme is None: deme_ = 0 
    if keys is None:
        keys = get_mutation_keys(pop,individuals,include_neutral,include_selected,remove_fixed,sort,min_sample_daf,min_pop_daf,deme)
    if isinstance(pop,Spop):
        return HaplotypeMatrix(fwdpp_haplotype_matrix[singlepop_t](deref((<Spop>pop).pop.get()),individuals,keys[0],keys[1],<size_t>deme_))
    if isinstance(pop,MlocusPop):
        return HaplotypeMatrix(fwdpp_haplotype_matrix[multilocus_t](deref((<MlocusPop>pop).pop.get()),individuals,keys[0],keys[1],<size_t>deme_))
    if isinstance(pop,MetaPop):
        if deme is None:
            raise RuntimeError("deme cannot be None")
        return HaplotypeMatrix(fwdpp_haplotype_matrix[metapop_t](deref((<MetaPop>pop).mpop.get()),individuals,keys[0],keys[1],<size_t>deme_))
    else:
        raise NotImplementedError("population type not supported")

def genotype_matrix(pop,list individuals,include_neutral=True,include_selected=True,remove_fixed=False,sort=True,min_sample_daf=None,min_pop_daf=None,deme=None,keys=None):
    """
    Generate a genotype matrix stored in a :class:`fwdpy.matrix.GenotypeMatrix`
    
    All arguments are identical to :func:`fwdpy.matrix.get_mutation_keys` with the exception of:

    :param keys: (None) Pre-calculated mutation keys.  Must be a tuple of lists of positive integers. Exceptions will be raised if keys are invalid or out of range.

    :rtype: :class:`fwdpy.matrix.GenotypeMatrix`
    """
    deme_=deme
    if deme is None:
        deme_=0
    if keys is None:
        keys = get_mutation_keys(pop,individuals,include_neutral,include_selected,remove_fixed,sort,min_sample_daf,min_pop_daf,deme)
    if isinstance(pop,Spop):
        return GenotypeMatrix(fwdpp_genotype_matrix[singlepop_t](deref((<Spop>pop).pop.get()),individuals,keys[0],keys[1],<size_t>deme_))
    if isinstance(pop,MlocusPop):
        return GenotypeMatrix(fwdpp_genotype_matrix[multilocus_t](deref((<MlocusPop>pop).pop.get()),individuals,keys[0],keys[1],<size_t>deme_))
    if isinstance(pop,MetaPop):
        if deme is None:
            raise RuntimeError("deme cannot be None")
        return GenotypeMatrix(fwdpp_genotype_matrix[metapop_t](deref((<MetaPop>pop).mpop.get()),individuals,keys[0],keys[1],<size_t>deme_))
    else:
        raise NotImplementedError("population type not supported")
