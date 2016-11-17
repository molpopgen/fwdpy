from fwdpy.fwdpp cimport popgenmut
from fwdpy.fwdpy cimport *

##Define some low-level functions that may be useful for others
cdef struct popgen_mut_data:
    double pos,s,h
    unsigned n,g,ftime
    bint neutral
    uint16_t label

cdef struct gamete_data:
    vector[popgen_mut_data] neutral,selected
    unsigned n

cdef struct diploid_data:
    gamete_data chrom0,chrom1
    double g,e,w,sh0,sh1
    int n0,n1

cdef struct diploid_mloc_data:
    vector[gamete_data] chrom0,chrom1
    double g,e,w
    vector[double] sh0,sh1
    vector[int] n0,n1

cdef class MutationView(object):
    cdef readonly float pos,s,h
    cdef readonly int n,g,label,key
    cdef readonly object ftime
    cdef readonly bint neutral
    cdef object __weakref__

cdef class GameteView(object):
    cdef readonly list neutral,selected
    cdef readonly int n,key 
    cdef object __weakref__

cdef class DiploidView(object):
    cdef readonly GameteView first,second
    cdef readonly float g,e,w
    cdef readonly int key
    cdef object __weakref__

cdef class MultiLocusDiploidView(object):
    cdef readonly list first,second
    cdef readonly float g,e,w
    cdef readonly int key
    cdef object __weakref__

cdef MutationView get_mutation(const popgenmut & m, size_t n,int key) 
cdef GameteView get_gamete(const gamete_t & g, const mcont_t & mutations, const mcounts_cont_t & mcounts,int key) 
cdef DiploidView get_diploid(const diploid_t & dip, const gcont_t & gametes, const mcont_t & mutations, const mcounts_cont_t & mcounts,int key) 
cdef MultiLocusDiploidView get_diploid_mloc(const dipvector_t & dip, const gcont_t & gametes, const mcont_t & mutations, const mcounts_cont_t & mcounts,int key) 
