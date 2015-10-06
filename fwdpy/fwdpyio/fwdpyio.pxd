from libcpp.string cimport string
from fwdpy.fwdpy cimport *

cdef extern from "serialize.hpp" namespace "fwdpy::serialize":
    string serialize_singlepop(const singlepop_t * pop)
    vector[shared_ptr[singlepop_t]] deserialize_singlepop(const vector[string] & strings)