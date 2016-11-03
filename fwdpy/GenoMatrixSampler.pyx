#This is a custom temporal sampler
#Currently, only singlepop and multilocus pop are supported

from fwdpy.gsl cimport *
from fwdpy.gsl_data_matrix cimport  *
from fwdpy.fwdpy cimport TemporalSampler,sampler_base,custom_sampler_data,singlepop_t,multilocus_t,uint
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from cython_gsl cimport gsl_matrix_alloc,gsl_matrix_set_zero,gsl_matrix_get
import numpy as np

cdef struct geno_matrix:
    vector[double] G,m
    size_t nrow,ncol

ctypedef vector[pair[uint,geno_matrix]] geno_matrix_final_t
ctypedef pair[bint,bint] geno_matrix_data_t
ctypedef custom_sampler_data[geno_matrix_final_t,geno_matrix_data_t] geno_matrix_sampler_t

cdef void  fill_matrix(const gsl_matrix_ptr_t & t, geno_matrix & gm) nogil:
    cdef size_t N=t.get().size1,ncol=t.get().size2
    if t.get().tda == t.get().size2:
        #Then physical allocation size is even
        #across rows
        gm.m.assign(<double*>t.get().data,<double*>(t.get().data+(N*ncol)))
    else:
        for i in range(N):
            for j in range(ncol):
                gm.m.push_back(gsl_matrix_get(t.get(),i,j))
    
cdef void singlepop_geno_matrix(const singlepop_t * pop, const unsigned generation, geno_matrix_final_t & f, geno_matrix_data_t & d) nogil:
    mut_keys = get_mut_keys(pop,d.first,d.second)
    cdef size_t i,j
    cdef size_t N = pop.diploids.size()
    cdef size_t ncol = mut_keys.size()+1
    cdef gsl_matrix_ptr_t t=gsl_matrix_ptr_t(<gsl_matrix*>gsl_matrix_alloc(N,ncol))
    update_matrix_counts(pop,mut_keys,t.get())
    gm.nrow=N
    gm.ncol=ncol
    cdef geno_matrix gm
    #Fill in genetic values
    for dip in range(pop.diploids.size()):
        gm.G.push_back(pop.diploids[dip].g)
    fill_matrix(t,gm)
    f.push_back(pair[uint,geno_matrix](generation,gm))

cdef void multilocus_geno_matrix(const multilocus_t * pop, const unsigned generation, geno_matrix_final_t & f, geno_matrix_data_t & d) nogil:
    mut_keys = get_mut_keys(pop,d.first,d.second)
    cdef size_t i,j
    cdef size_t N = pop.diploids.size()
    cdef size_t ncol = mut_keys.size()+1
    cdef gsl_matrix_ptr_t t=gsl_matrix_ptr_t(<gsl_matrix*>gsl_matrix_alloc(N,ncol))
    update_matrix_counts[multilocus_t](pop,mut_keys,t.get())
    gm.nrow=N
    gm.ncol=ncol
    cdef geno_matrix gm
    #Fill in genetic values
    for dip in range(pop.diploids.size()):
        gm.G.push_back(pop.diploids[dip][0].g)
    fill_matrix(t,gm)
    f.push_back(pair[uint,geno_matrix](generation,gm))

cdef class GenoMatrixSampler(TemporalSampler):
    """
    A :class:`fwdpy.fwdpy.TemporalSampler` to record a 0,1,2-encoded matrix of genotypes at sites 
    affecting fitness/trait values.
    """
    def __cinit__(self,unsigned n,bint sort_freq,bint sort_esize):
        """
        Constructor.

        :param n: A length.  Must correspond to number of simulations that will be run simultaneously.
        :param sort_freq: If True, then columns (mutations) will be sorted in descending order of allele freuqency (left to right)
        :param sort_esize: If True, then columns (mutations) will be sorted in descending order of |esize| (left to right)

        If sort_freq is True and sort_esize is True, then mutations are first sorted by frequency, then sorted by |esize| within each 
        frequency bin.
        """
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[geno_matrix_sampler_t](new geno_matrix_sampler_t(&singlepop_geno_matrix,geno_matrix_data_t(sort_freq,sort_esize))))
            (<geno_matrix_sampler_t*>self.vec[i].get()).register_callback(&multilocus_geno_matrix) 
    def get(self,bint keep_origin = False):
        rv=[]
        for i in range(self.vec.size()):
            temp=(<geno_matrix_sampler_t*>self.vec[i].get()).final()
            n=np.matrix(temp.m,temp.ncol,temp.nrow)
            if keep_origin is False:
                n=np.delete(n,[0],axis=1)
            n=np.insert(n,0,temp.G)
            rv.append(temp)
        return rv
