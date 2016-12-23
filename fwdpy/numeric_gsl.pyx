from cython_gsl.gsl_blas cimport *
from cython_gsl.gsl_linalg cimport * 
from cython_gsl.gsl_math cimport * 

#These are unique_ptr wrappers around
#the GSL pointer types.
from fwdpy.gsl cimport gsl_vector_ptr_t,gsl_matrix_ptr_t 

cdef pair[double,vector[double]] sum_of_squares(const gsl_vector * v,
                                                gsl_matrix * m) nogil:
    cdef gsl_vector_ptr_t TAU,SUMS
    cdef gsl_matrix_ptr_t Q,R
    TAU.reset(gsl_vector_alloc(min(m.size1,m.size2)))
    SUMS.reset(gsl_vector_alloc(m.size1))
    Q.reset(gsl_matrix_alloc(m.size1,m.size1))
    R.reset(gsl_matrix_alloc(m.size1,m.size2))

    gsl_linalg_QR_decomp(m,TAU.get())
    gsl_linalg_QR_unpack(m,TAU.get(),Q.get(),R.get())
    gsl_blas_dgemv(CblasTrans,1.0,Q.get(),v,0.0,SUMS.get())

    cdef pair[double,vector[double]] rv

    cdef size_t i
    for i in range(0,m.size2): 
        rv.second.push_back(gsl_pow_2(gsl_vector_get(SUMS.get(),i+1)))
   
    cdef size_t DF = m.size2-1
    RSS=0.0
    for i in range(DF+1,SUMS.get().size):
        RSS += gsl_pow_2(gsl_vector_get(SUMS.get(),i))

    cdef double sqi=0.
    rv.first = RSS
    for sqi in rv.second: rv.first += sqi
    return rv

