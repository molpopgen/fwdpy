# distutils: language = c++
# distutils: sources = fwdpy/fwdpy/sample.cc fwdpy/fwdpy/deps.cc 
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

include "classes.pyx"
include "sampling.pyx"
include "evolve_regions.pyx"
include "regions.pyx"
include "copy.pyx"
include "views.pyx"
include "view_fixations.pyx"
include "debug.pyx"
include "temporal_samplers.pyx"
include "add_mutations.pyx"

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

def get_includes():
    """
    Returns absolute path to location of fwdpy headers
    """
    import os,fwdpy
    return os.path.abspath(fwdpy.__file__).split('lib')[0]+'include/fwdpy'

