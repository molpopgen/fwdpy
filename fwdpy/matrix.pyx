from fwdpy cimport *
from cython.operator cimport dereference as deref
from fwdpy cimport uint
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

cdef key_pair apply_min_daf(key_pair & keys,const double n,const double x) nogil:
    cdef key_pair rv
    for i in keys.first:
        if <double>i.second/n >= x:
            rv.first.push_back(i)
    for i in keys.second:
        if <double>i.second/n >= x:
            rv.second.push_back(i)
    return rv

cdef class DataMatrix(object):
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

cdef sort_keys_by_position(key_pair & keys,const mcont_t & mutations):
    #positions of neutral mutations
    pos = [mutations[i[0]].pos for i in keys.first]         
    keys.first=[i for j,i in sorted(zip(pos,keys.first),key=lambda pair: pair[0])]
    pos = [mutations[i[0]].pos for i in keys.second]         
    keys.second=[i for j,i in sorted(zip(pos,keys.second),key=lambda pair: pair[0])]

def get_mutation_keys(pop,list individuals,include_neutral=True,include_selected=True,remove_fixed=False,sort=True,min_daf=None,deme=None):
    deme_=deme
    if deme is None:
        deme_=0
    if min_daf is not None:
        if float(min_daf) < 0. or float(min_daf) >= 1.:
            raise RuntimeError("Minimum derived allele frequency must be a float 0.0 <= 1 < 1.0")
    #Get the mutation keys first
    cdef key_pair keys
    if isinstance(pop,Spop):
        keys = mutation_keys[singlepop_t](deref((<Spop>pop).pop.get()),individuals,include_neutral,include_selected,<size_t>deme_)
        if sort is True:
            sort_keys_by_position(keys,(<Spop>pop).pop.get().mutations)
    elif isinstance(pop,MlocusPop):
        keys = mutation_keys[multilocus_t](deref((<MlocusPop>pop).pop.get()),individuals,include_neutral,include_selected,<size_t>deme_)
        if sort is True:
            sort_keys_by_position(keys,(<Spop>pop).pop.get().mutations)
    elif isinstance(pop,MetaPop):
        if deme is None:
            raise RuntimeError("deme cannot be none")
        keys = mutation_keys[metapop_t](deref((<MetaPop>pop).mpop.get()),individuals,include_neutral,include_selected,<size_t>deme_)
        if sort is True:
            sort_keys_by_position(keys,(<Spop>pop).pop.get().mutations)
    #Apply filters to keys
    if remove_fixed is True:
        keys = remove_fixed_keys(keys,2*len(individuals))
    if min_daf is not None:
        keys = apply_min_daf(keys,2*len(individuals),min_daf)
    return keys

def haplotype_matrix(pop,list individuals,include_neutral=True,include_selected=True,remove_fixed=False,sort=True,min_daf=None,deme=None,keys=None):
    deme_=deme
    if keys is None:
        keys = get_mutation_keys(pop,individuals,include_neutral,include_selected,remove_fixed,sort,min_daf,deme)
    if isinstance(pop,Spop):
        return DataMatrix(fwdpp_haplotype_matrix[singlepop_t](deref((<Spop>pop).pop.get()),individuals,keys[0],keys[1],<size_t>deme_))

def genotype_matrix(pop,list individuals,include_neutral=True,include_selected=True,remove_fixed=False,sort=True,min_daf=None,deme=None,keys=None):
    deme_=deme
    if deme is None:
        deme_=0
    if keys is None:
        keys = get_mutation_keys(pop,individuals,include_neutral,include_selected,remove_fixed,sort,min_daf,deme)
    if isinstance(pop,Spop):
        return DataMatrix(fwdpp_genotype_matrix[singlepop_t](deref((<Spop>pop).pop.get()),individuals,keys[0],keys[1],<size_t>deme_))
