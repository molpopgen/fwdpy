from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.map cimport map

cdef extern from "libseq.hpp" namespace "fwdpy::libseq":
    map[string,double] libseq_basic_stats( const vector[pair[double,string]] & __data ) except +
    vector[vector[pair[double,string]]] sliding_windows_cpp( const vector[pair[double,string]] & sample, const double window_size, const double steplen, const double starting_pos, const double ending_pos ) except +
