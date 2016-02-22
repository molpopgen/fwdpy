cdef extern from "gsl/gsl_rng.h" nogil:
    ctypedef struct gsl_rng
    ctypedef struct gsl_ran_discrete_t
    size_t gsl_ran_discrete ( gsl_rng * ,  gsl_ran_discrete_t *)
    gsl_ran_discrete_t * gsl_ran_discrete_preproc(size_t k,const double * p)
