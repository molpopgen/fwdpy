#!/usr/bin/env python

from __future__ import print_function
from distutils.core import setup,Extension
from Cython.Build import cythonize
import importlib
import platform, glob, sys, subprocess

if 'install' in sys.argv:
    print("Sorry, you shouldn't be trying to install these!")
    sys.exit(0)

if '--run' in sys.argv:
    run = glob.glob('run_test_extensions/run*.py')
    for i in run:
        i=i.replace('/','.').replace(".py",'')
        print("running test "+i)
        importlib.import_module(i)
    sys.exit(0)

DIR='test_fwdpy_extensions'
print (DIR+'/'+'*.pyx')
pyxnames = glob.glob('test_fwdpy_extensions/*.pyx')
print(pyxnames)
modules=[]
for i in pyxnames:
    print(i)
    i=i.rsplit('.')[0]
    print (i)
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


