# distutils: language = c++
# distutils: sources = fwdpy/fwdpy/sample.cc fwdpy/fwdpy/deps.cc fwdpy/fwdpy/evolve_regions.cc 
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

## for DataFrame
import pandas

include "classes.pyx"
include "sampling.pyx"
include "evolve_regions.pyx"
include "regions.pyx"
include "copy.pyx"
include "views.pyx"
include "debug.pyx"
include "ages.pyx"

def pkg_version():
    """
    Return version numbers of this package

    This function is very handy when reporting bugs!

    :returns: dict
    """
    cdef vector[string] v = fwdpy_version()
    return ({'fwdpy':v[0]})

def cite():
    """
    Returns how to cite this package
    """
    fwdpy_citation()
