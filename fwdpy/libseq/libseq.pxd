from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.map cimport map

cdef extern from "libseq/libseq.hpp" namespace "fwdpy::libseq":
    map[string,double] libseq_basic_stats( const vector[pair[double,string]] & __data )
    map[string,double] libseq_extra_ld_stats( const vector[pair[double,string]] & __data, const double minfreq, const double binsize) except +
    map[string,double] libseq_extra_ld_stats( const vector[pair[double,string]] & __data, const double minfreq, const double binsize, const vector[double] & gmap) except +
    vector[vector[pair[double,string]]] sliding_windows_cpp( const vector[pair[double,string]] & sample, const double window_size, const double steplen, const double starting_pos, const double ending_pos ) except +
