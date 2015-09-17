from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.map cimport map

cdef extern from "summstats.hpp" namespace "fwdpy::libseq":
    map[string,double] libseq_basic_stats( const vector[pair[double,string]] & __data )
