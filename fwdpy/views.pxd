from fwdpy.fwdpp cimport popgenmut
from fwdpy.fwdpy cimport *

##Define some low-level functions that may be useful for others
#cdef struct popgen_mut_data:
#    double pos,s,h
#    unsigned n,g,ftime
#    bint neutral
#    uint16_t label
#
#cdef struct gamete_data:
#    vector[popgen_mut_data] neutral,selected
#    unsigned n
#
#cdef struct diploid_data:
#    gamete_data chrom0,chrom1
#    double g,e,w,sh0,sh1
#    int n0,n1
#
#cdef struct diploid_mloc_data:
#    vector[gamete_data] chrom0,chrom1
#    double g,e,w
#    vector[double] sh0,sh1
#    vector[int] n0,n1

cdef class MutationView(object):
    cdef readonly object pos
    """Position"""
    cdef readonly object s
    """Effect size/selection coefficient"""
    cdef readonly object h
    """Dominance"""
    cdef readonly object n
    """Number of occurrences in the population"""
    cdef readonly object g
    """Origin time (generation mutation first appeared)"""
    cdef readonly object label
    """Mutation label"""
    cdef readonly object mut_key
    """Index of the mutation in the mutation container.  Value is only meaningful for the generation when view was taken"""
    cdef readonly object ftime
    """Fixation time"""
    cdef readonly object neutral
    """True if mutation neutral, False otherwise"""
    cdef object __weakref__

cdef class GameteView(object):
    cdef readonly list neutral
    """:class:`fwdpy.views.MutationView` for neutral variants"""
    cdef readonly list selected
    """:class:`fwdpy.views.MutationView` for selected variants"""
    cdef readonly int n
    """Number of occurrences of this gamete in the population"""
    cdef readonly int gam_key 
    """Index of gamete in gamete container. Value is only meaningful for the generation when view was taken"""  
    cdef object __weakref__

cdef class DiploidView(object):
    cdef readonly GameteView first
    """A :class:`fwdpy.views.GameteView`"""
    cdef readonly GameteView second
    """A :class:`fwdpy.views.GameteView`"""
    cdef readonly float g
    """Genetic value"""
    cdef readonly float e
    """Random/environmental value"""
    cdef readonly float w
    """Fitness"""
    cdef readonly int dip_key
    """Index of individual in diploid container. Value is only meaningful for the generation when view was taken"""  
    cdef object __weakref__

cdef class MultiLocusDiploidView(object):
    cdef readonly GameteView first
    """A list of :class:`fwdpy.views.GameteView`"""
    cdef readonly GameteView second
    """A list of A :class:`fwdpy.views.GameteView`"""
    cdef readonly float g
    """Genetic value"""
    cdef readonly float e
    """Random/environmental value"""
    cdef readonly float w
    """Fitness"""
    cdef readonly int dip_key
    """Index of individual in diploid container. Value is only meaningful for the generation when view was taken"""  
    cdef object __weakref__

cdef MutationView get_mutation(const popgenmut & m, size_t n,int key) 
cdef MutationView empty_MutationView()
cdef GameteView get_gamete(const gamete_t & g, const mcont_t & mutations, const mcounts_cont_t & mcounts,int key) 
cdef DiploidView get_diploid(const diploid_t & dip, const gcont_t & gametes, const mcont_t & mutations, const mcounts_cont_t & mcounts,int key) 
cdef MultiLocusDiploidView get_diploid_mloc(const dipvector_t & dip, const gcont_t & gametes, const mcont_t & mutations, const mcounts_cont_t & mcounts,int key) 
