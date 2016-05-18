from libcpp.vector cimport vector
from fwdpy.fwdpy cimport singlepop_t

cdef extern from "gwas_genotype_matrix.hpp" namespace "fwdpy::gwas" nogil:
    cdef struct genotype_matrix:
        vector[double] G
        vector[double] E
        vector[double] neutral
        vector[double] causative
        size_t N
        size_t n_neutral
        size_t n_causative

    genotype_matrix make_geno_matrix(const singlepop_t * pop,
                                     bint maf) except +
        
        
        
