from fwdpy cimport *
from cython.operator cimport dereference as deref
from fwdpy cimport uint

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

def get_mutation_keys(pop,list individuals,include_neutral=True,include_selected=True,remove_fixed=False,min_daf=None,deme=None):
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
    elif isinstance(pop,MlocusPop):
        keys = mutation_keys[multilocus_t](deref((<MlocusPop>pop).pop.get()),individuals,include_neutral,include_selected,<size_t>deme_)
    elif isinstance(pop,MetaPop):
        if deme is None:
            raise RuntimeError("deme cannot be none")
        keys = mutation_keys[metapop_t](deref((<MetaPop>pop).mpop.get()),individuals,include_neutral,include_selected,<size_t>deme_)
        
    #Apply filters to keys
    if remove_fixed is True:
        keys = remove_fixed_keys(keys,2*individuals.size())
    if min_daf is not None:
        keys = apply_min_daf(keys,2*individuals.size(),min_daf)

def hap_matrix(pop,list individuals,include_neutral=True,include_selected=True,remove_fixed=False,min_daf=None,deme=None,keys=None):
    deme_=deme
    if deme is None:
        deme_=0
    if keys is None:
        keys_ = get_mutation_keys(pop,individuals,include_neutral,include_selected,remove_fixed,min_daf,deme)

    if isinstance(pop,Spop):
        return haplotype_matrix[singlepop_t](deref((<Spop>pop).pop.get()),individuals,keys.first,keys.second,<size_t>deme_)

def geno_matrix(pop,list individuals,include_neutral=True,include_selected=True,remove_fixed=False,min_daf=None,deme=None,keys=None):
    deme_=deme
    if deme is None:
        deme_=0
    if keys is None:
        keys_ = get_mutation_keys(pop,individuals,include_neutral,include_selected,remove_fixed,min_daf,deme)

    if isinstance(pop,Spop):
        return genotype_matrix[singlepop_t](deref((<Spop>pop).pop.get()),individuals,keys.first,keys.second,<size_t>deme_)
