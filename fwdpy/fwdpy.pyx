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
    return os.path.dirname(fwdpy.__file__)+'/headers'

def get_fwdpp_includes():
    """
    Returns absolute path to location of the fwdpp headers 
    installed along with fwdpy.
    """
    return get_includes()+'/fwdpp'

def make_pyxbld(pkgname):
    """
    When writing Cython extensions to be used as "plugins",
    each .pyx file needs a .pyxbld file.  This function
    auto-generates that file.

    :param pkgname: For plugin.pyx, pkgname=plugin

    .. note:: The output file will be over-written if it exists.
    """
    ofile = pkgname + '.pyxbld'
    text="""
import fwdpy as fp
fwdpy_includes=fp.get_includes()
fwdpp_includes=fp.get_fwdpp_includes()
def make_ext(modname, pyxfilename):
    from distutils.extension import Extension
    return Extension(name=modname,
                     sources=[pyxfilename],
                     language='c++',
		     include_dirs=[fwdpy_includes,fwdpp_includes],
		     extra_compile_args=['-std=c++11'],
		     libraries=['sequence','gsl','gslcblas','pthread'])
"""
    with open(ofile,'w') as output:
        output.write(text.format())
