from fwdpy.structs cimport VAcum
from libcpp.vector cimport vector

cdef extern from "sampler_additive_variance.hpp" namespace "fwdpy" nogil:
    vector[VAcum] additive_variance_wrapper[POPTYPE](const POPTYPE * p)
