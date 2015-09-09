# distutils: language = c++
# distutils: sources = fwdpy/src/sample.cpp fwdpy/src/neutral.cpp fwdpy/src/deps.cc fwdpy/src/callbacks.cc fwdpy/src/evolve_regions.cc
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

## for DataFrame
import pandas

include "classes.pyx"
include "evolve_simple.pyx"
include "sampling.pyx"
include "sregionCallbacks.pyx"
include "evolve_regions.pyx"

def pkg_dependencies():
    """
    Return version numbers of dependencies.

    This function is very handy when reporting bugs!

    :returns: dict
    """
    cdef vector[string] deps = dependencies()
    return ({'fwdpp':deps[0],'libsequence':deps[1],'GSL':deps[2]})

