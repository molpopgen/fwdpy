from fwdpy.fwdpp cimport popgenmut
from fwdpy.fwdpy cimport *

cdef popgen_mut_data get_mutation( const popgenmut & m, size_t n) nogil
cdef gamete_data get_gamete( const gamete_t & g, const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil
cdef diploid_data get_diploid( const diploid_t & dip, const gcont_t & gametes, const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil
cdef diploid_mloc_data get_diploid_mloc( const dipvector_t & dip, const gcont_t & gametes, const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil

