# distutils: language = c++
# distutils: sources = fwdpy/fwdpy/poptypes.cpp fwdpy/fwdpy/sample.cpp fwdpy/fwdpy/neutral.cpp fwdpy/fwdpy/deps.cc fwdpy/fwdpy/evolve_regions.cc
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

## for DataFrame
import pandas

include "classes.pyx"
include "evolve_simple.pyx"
include "sampling.pyx"
include "evolve_regions.pyx"
include "regions.pyx"

def pkg_dependencies():
    """
    Return version numbers of dependencies.

    This function is very handy when reporting bugs!

    :returns: dict
    """
    cdef vector[string] deps = fwdpy_dependencies()
    return ({'fwdpp':deps[0],'libsequence':deps[1],'GSL':deps[2]})

def pkg_version():
    """
    Return version numbers of this package

    This function is very handy when reporting bugs!

    :returns: dict
    """
    cdef vector[string] v = fwdpy_version()
    return ({'fwdpy':v[0]})
