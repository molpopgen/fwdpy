#!/usr/bin/env python

from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(name="fwdpy",
        version="0.0.1",
        author='Kevin R. Thornton',
        author_email='krthornt@uci.edu',
        license = "GSL",
        url = "http://github.com/molpopgen/fwdpy",
        ext_modules = cythonize(Extension(
            "fwdpy",                                # the extesion name
            sources=["fwdpy.pyx"], # the Cython source and additional C++ source files
            language="c++",                        # generate and compile C++ code
            extra_compile_args=["-std=c++11","-I.","-I.."],
            extra_link_args=["-std=c++11"],
            libraries=["gsl","gslcblas","tcmalloc","sequence"]))
      )
