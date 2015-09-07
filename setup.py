#!/usr/bin/env python

# setup.py.in.distutils
#
# Copyright 2012, 2013 Brandon Invergo <brandon@invergo.net>
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.


from distutils.core import setup, Extension
from Cython.Build import cythonize
import platform


if platform.system() == 'Linux':
    doc_dir = '/usr/local/share/doc/fwdpy'
else:
    try:
        from win32com.shell import shellcon, shell
        homedir = shell.SHGetFolderPath(0, shellcon.CSIDL_APPDATA, 0, 0)
        appdir = 'fwdpy'
        doc_dir = os.path.join(homedir, appdir)
    except:
        pass

long_desc = \
"""
"""

setup(name='fwdpy',
      version='0.0.1',
      author='',
      author_email='krthornt@uci.edu',
      maintainer='',
      maintainer_email='',
      url='',
      description="""""",
      long_description=long_desc,
      download_url='',
      classifiers=[''],
      platforms=[''],
      license='',
      requires=['pandas'],
      provides=['fwdpy'],
      obsoletes=['none'],
      packages=[],
      py_modules=[],
      scripts=[],
      data_files=[(doc_dir, ['COPYING', 'README.md'])],
      package_data={},
      ext_modules=cythonize(Extension("fwdpy.fwdpy",
          sources=["fwdpy/fwdpy.pyx"], # the Cython source and additional C++ source files
          language="c++",                        # generate and compile C++ code
          include_dirs=['.'],
          extra_compile_args=["-std=c++11","-I.","-I..","-Ifwdpy"],
          extra_link_args=["-std=c++11"],
          libraries=["gsl","gslcblas","tcmalloc","sequence"]),
          version='0.0.1'
          )
     )
