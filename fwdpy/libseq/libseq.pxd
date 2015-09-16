from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

cdef extern from "summstats.hpp" namespace "fwdpy::libseq":
    double tajd( const vector[pair[double,string]] & __data )
