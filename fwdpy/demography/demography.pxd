from fwdpy.fwdpp cimport *
from fwdpy.fwdpp_internal cimport gsl_ran_discrete_t_ptr
from fwdpy.fwdpy cimport metapop_t,singlepop_t
from libcpp.vector cimport vector

cdef extern from "metapop.hpp" namespace "fwdpy" nogil:
    void re_init_mpop( metapop_t * mpop, const singlepop_t * pop)
    void copy_deme( metapop_t * mpop, const size_t i ) except +
    void remove_deme( metapop_t * mpop, const size_t i ) except +
    void merge_demes(metapop_t  * mpop, const size_t i, const size_t j) except +
    void split_deme(const gsl_rng * r, metapop_t * mpop, const size_t i, const unsigned N_new, const bint replacement ) except +
    void admix_demes(const gsl_rng * r, metapop_t * mpop, const size_t i, const size_t j, const double prop_i,const unsigned N_new, const bint replacement) except +
    void swap_demes(metapop_t * mpop, const size_t i, const size_t j) except +

cdef extern from "demography_migrates.hpp" namespace "fwdpy::demography" nogil:
    cdef cppclass migrates:
        vector[gsl_ran_discrete_t_ptr] lookup
        migrates()
        migrates(const vector[vector[double]] & weights)
        migrates(migrates &)

    migrates make_migrates(const vector[vector[double]] & weights)
