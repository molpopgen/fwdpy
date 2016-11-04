# distutils: language = c++
#This is a custom temporal sampler
#Currently, only singlepop and multilocus pop are supported

from fwdpy.gsl cimport *
from fwdpy.gsl_data_matrix cimport  *
from fwdpy.fwdpy cimport TemporalSampler,sampler_base,custom_sampler_data,singlepop_t,multilocus_t,uint
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
from libcpp.utility cimport pair
from cython_gsl cimport gsl_matrix_alloc,gsl_matrix_set_zero,gsl_matrix_get
from cython.operator cimport dereference as deref
from cython.parallel import parallel, prange
import numpy as np

ctypedef vector[pair[uint,shared_ptr[geno_matrix]]] geno_matrix_final_t
ctypedef pair[bint,bint] geno_matrix_data_t
ctypedef custom_sampler_data[geno_matrix_final_t,geno_matrix_data_t] geno_matrix_sampler_t

cdef void singlepop_geno_matrix(const singlepop_t * pop, const unsigned generation, geno_matrix_final_t & f, geno_matrix_data_t & d) nogil:
    mut_keys = get_mut_keys(pop,d.first,d.second)
    cdef size_t i,j
    cdef size_t N = pop.diploids.size()
    cdef size_t ncol = mut_keys.size()+1
    cdef shared_ptr[geno_matrix] gm=shared_ptr[geno_matrix](new geno_matrix(N,ncol))
    gsl_matrix_set_zero(gm.get().m.get())
    update_matrix_counts(pop,mut_keys,gm.get().m.get())
    #Fill in genetic values
    for dip in range(pop.diploids.size()):
        gm.get().G.push_back(pop.diploids[dip].g)
    f.push_back(pair[uint,shared_ptr[geno_matrix]](generation,gm))

cdef void multilocus_geno_matrix(const multilocus_t * pop, const unsigned generation, geno_matrix_final_t & f, geno_matrix_data_t & d) nogil:
    mut_keys = get_mut_keys(pop,d.first,d.second)
    cdef size_t i,j
    cdef size_t N = pop.diploids.size()
    cdef size_t ncol = mut_keys.size()+1
    cdef shared_ptr[geno_matrix] gm=shared_ptr[geno_matrix](new geno_matrix(N,ncol))
    gsl_matrix_set_zero(gm.get().m.get())
    update_matrix_counts[multilocus_t](pop,mut_keys,gm.get().m.get())
    #Fill in genetic values
    for dip in range(pop.diploids.size()):
        gm.get().G.push_back(pop.diploids[dip][0].g)
    f.push_back(pair[uint,shared_ptr[geno_matrix]](generation,gm))

cdef class GenoMatrixSampler(TemporalSampler):
    """
    A :class:`fwdpy.fwdpy.TemporalSampler` to record a 0,1,2-encoded matrix of genotypes at sites 
    affecting fitness/trait values.
    """
    def __cinit__(self,unsigned n,bint sort_freq=False,bint sort_esize = False):
        """
        Constructor.

        :param n: A length.  Must correspond to number of simulations that will be run simultaneously.
        :param sort_freq: (False) If True, then columns (mutations) will be sorted in descending order of allele freuqency (left to right)
        :param sort_esize: (False) If True, then columns (mutations) will be sorted in descending order of |esize| (left to right)

        If sort_freq is True and sort_esize is True, then mutations are first sorted by frequency, then sorted by |esize| within each 
        frequency bin.
        """
        for i in range(n):
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[geno_matrix_sampler_t](new geno_matrix_sampler_t(&singlepop_geno_matrix,geno_matrix_data_t(sort_freq,sort_esize))))
            (<geno_matrix_sampler_t*>self.vec[i].get()).register_callback(&multilocus_geno_matrix) 
    cdef make_numpy_matrix(self,vector[pair[uint,shared_ptr[geno_matrix]]].iterator beg,bint keep_origin):
        cdef size_t N = deref(beg).second.get().nrow
        cdef size_t S = deref(beg).second.get().ncol
        n=np.array([0.]*(N*S),dtype=np.float64)
        n=np.reshape(n,[N,S])
        cdef size_t i,j
        for i in range(N):
            for j in range(S):
                n[i,j]=gsl_matrix_get(deref(beg).second.get().m.get(),i,j)
        return n
        #n=np.array(deref(beg).second.get().m)
        #n=np.reshape(n,[deref(beg).second.get().nrow,deref(beg).second.get().ncol])
        #if keep_origin is False:
        #    #Replace column of 1. with
        #    #genetic values
        #    n[:,0]=deref(beg).second.get().G
        #else:
            #Insert genetic values as column 0
        #    n=np.insert(n,0,deref(beg).second.get().G)
        #return n 
    def get(self,bint keep_origin = False):
        rv=[]
        cdef vector[pair[uint,shared_ptr[geno_matrix]]].iterator beg,end
        for i in range(self.vec.size()):
            beg = (<geno_matrix_sampler_t*>self.vec[i].get()).f.begin()
            end = (<geno_matrix_sampler_t*>self.vec[i].get()).f.end()
            while beg<end:
                rv.append((deref(beg).first,self.make_numpy_matrix(beg,keep_origin)))
                beg+=1
        return rv
    cdef void tofile_details_task(self,const geno_matrix_final_t & f,cppstring stub,int repid,bint keep_origin) nogil:
        cdef vector[pair[uint,shared_ptr[geno_matrix]]].const_iterator beg,end
        beg=f.const_begin()
        end=f.const_end()
        while beg<end:
            write_geno_matrix(deref(beg).second.get(),deref(beg).first,stub,repid,keep_origin)
            beg+=1
    cdef void tofile_details(self,cppstring stub,int repstart,bint keep_origin) nogil:
        cdef int task
        for task in prange(self.vec.size(),nogil=True,schedule='static'):
            self.tofile_details_task((<geno_matrix_sampler_t*>self.vec[task].get()).f,stub,repstart+task,keep_origin)
    def tofile(self,stub,repstart=0,keep_origin=False):
        self.tofile_details(stub,repstart,keep_origin)
