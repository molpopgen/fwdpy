#!/usr/bin/env python

from __future__ import print_function
from distutils.core import setup,Extension
from Cython.Build import cythonize
    
import platform, glob, sys, subprocess

if 'install' in sys.argv:
    print("Sorry, you shouldn't be trying to install these!")
    sys.exit(0)
          
pyxnames = glob.glob('*.pyx')

modules=[]
for i in pyxnames:
    i=i.rsplit('.')[0]
    modules.append(i)

extensions=[]
provided=[]

for i in modules:
    extensions.append(Extension(i,
                                sources=[i+".pyx"],
                                language="c++",                  
                                extra_compile_args=["-std=c++11","-fopenmp"],  
                                extra_link_args=["-std=c++11","-fopenmp"],
                                libraries=["sequence"])
                                )
    
extensions=cythonize(extensions)

setup(license='GPL >= 2',
      ext_modules=extensions
)
