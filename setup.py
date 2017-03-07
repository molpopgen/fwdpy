#!/usr/bin/env python

# setup.py.in.distutils
#
# Copyright 2012, 2013 Brandon Invergo <brandon@invergo.net>
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.

from __future__ import print_function
from distutils.core import setup, Extension
#@CYTHONIMPORT@
import platform, glob, sys, subprocess, os


##Check C/C++ dependencies
FWDPP_V=""
GSL_V=""
LIBSEQ_V=""
CUSTOM_DIPLOID_BASE=1
##fwdpp
try:
    proc = subprocess.Popen(['fwdppConfig','--version'],stdout=subprocess.PIPE)
    (out,err) = proc.communicate()
    version = out.decode('utf-8').rstrip()
    print ("fwdpp version",version," found.")
    FWDPP_V=version
    if version < '0.5.4':
        sys.exit("fwdpp >= ,'0.5.4' required, but ",version, "found.")
    if version >= '0.5.5':
        CUSTOM_DIPLOID_BASE=0
except:
    sys.exit("fwdppConfig not found.  Please install fwdpp (http://github.com/molpopgen/fwdpp)")

##libseqeuence
try:
    proc = subprocess.Popen(['libsequenceConfig','--version'],stdout=subprocess.PIPE)
    (out,err) = proc.communicate()
    version = out.decode('utf-8').rstrip()
    print ("libsequence version",version," found.")
    LIBSEQ_V=version
except:
    sys.exit("libsequenceConfig not found.  Please install libsequence (http://github.com/molpopgen/libsequence)")

##GSL
try:
    proc = subprocess.Popen(['gsl-config','--version'],stdout=subprocess.PIPE)
    (out,err) = proc.communicate()
    version = out.decode('utf-8').rstrip()
    GSL_V=version
    print ("GSL version ",version," found.")
except:
    sys.exit("gsl-config not found.  Please install the GNU Scientific Library")

#Are we gonna build using Cython or not?  Default is not to,
#which allows us to ship this in a standard way.
if '--use-cython' in sys.argv:
    USE_CYTHON = True
    sys.argv.remove('--use-cython')
else:
    USE_CYTHON = False

if '--qtrait' in sys.argv:
    QTRAIT=True
    sys.argv.remove('--qtrait')
else:
    QTRAIT=False

##Set up our dependent libraries
GSLLIBS=["gsl","gslcblas"]
#MEMLIBS=None
LIBS=["sequence","gsl","gslcblas","z"]

#USETCMALLOC=False

#if str('@TCMALLOC@') == str('tcmalloc'):
#    USETCMALLOC=True
#    MEMLIBS=['@TCMALLOC@']

#if str('@TBBPRX@') == str('tbbmalloc_proxy'):
#    if USETCMALLOC == True:
#        raise RuntimeError("cannot use both tcmalloc and tbbmalloc_proxy")
#    MEMLIBS=['@TBBPRX@']

#if MEMLIBS is not None:
#    for i in MEMLIBS:
#        LIBS.append(i)

##Can we compile program based on dependencies?
print("Attempting to compile and link a test program using fwdpy's dependencies")
try:
    proc = subprocess.check_output(['make','-f','check_deps/Makefile','clean'])
#    if MEMLIBS is not None:
#        os.environ['TESTMEMLIBS'] = '-l'+str(MEMLIBS[0])
    proc = subprocess.check_output(['make','-f','check_deps/Makefile'])
    proc = subprocess.check_output(['make','-f','check_deps/Makefile','clean'])
except subprocess.CalledProcessError as e:
    print (e.returncode)
    sys.exit("Could not compile and link test programs!")

print("Success!")

##These check dependencies:
try:
    import pandas
except ImportError:
    sys.exit("import pandas failed.  Please install pandas")
try:
    import numpy
except ImportError:
    sys.exit("import numpy failed.  Please install numpy")

long_desc = open("README.rst").read()

#EXTENSION="@EXTENSION@"

GLOBAL_COMPILE_ARGS=['-std=c++11','-fopenmp',
                     str('-DPACKAGE_VERSION=')+'"0.0.4-rc3"',
                     '-DHAVE_INLINE'
]
if CUSTOM_DIPLOID_BASE != 0:
    GLOBAL_COMPILE_ARGS.append('-DCUSTOM_DIPLOID_BASE')

LINK_ARGS=["-std=c++11",'-fopenmp']
GLOBAL_INCLUDES=['.','..','fwdpy/headers','fwdpy/headers/fwdpp']

#Cython generates .cpp files from the .pyx files
#In general, we wisht to complies the .cpp files
#that are already here.  But, if we are modifying
#the package, we need to re-generate them, and thus
#re-process the .pyx files.
#This variable handes that choice, based on
#use of --use-cython or not (see above).
if len(glob.glob("fwdpy/*cpp"))==0:
    USE_CYTHON=1
EXTENSION = '.pyx' if USE_CYTHON else '.cpp'

FWDPY_MAIN_LIBS=LIBS
FWDPY_MAIN_LIBS.append('sqlite3')
extensions = [
    Extension("fwdpy.fwdpy",
              sources=["fwdpy/fwdpy"+EXTENSION]+glob.glob("fwdpy/fwdpy/*.cc"), # the Cython source and additional C++ source files
              language="c++",                        # generate and compile C++ code
              include_dirs=GLOBAL_INCLUDES,
              extra_compile_args=GLOBAL_COMPILE_ARGS,
              extra_link_args=LINK_ARGS,
              libraries=LIBS),
    Extension("fwdpy.internal.internal",
              sources=["fwdpy/internal/internal"+EXTENSION]+glob.glob("fwdpy/internal/*.cc"),
              language="c++",
              include_dirs=GLOBAL_INCLUDES,
              extra_compile_args=GLOBAL_COMPILE_ARGS,
              extra_link_args=LINK_ARGS,
              libraries=LIBS),
    Extension("fwdpy.fwdpyio.fwdpyio",
              sources=["fwdpy/fwdpyio/fwdpyio"+EXTENSION]+glob.glob("fwdpy/fwdpyio/*.cc"),
              language="c++",
              include_dirs=GLOBAL_INCLUDES,
              extra_compile_args=GLOBAL_COMPILE_ARGS,
              extra_link_args=LINK_ARGS,
              libraries=LIBS),
    Extension("fwdpy.demography.demography",
              sources=["fwdpy/demography/demography"+EXTENSION]+glob.glob("fwdpy/demography/*.cc"),
              language="c++",
              include_dirs=GLOBAL_INCLUDES,
              extra_compile_args=GLOBAL_COMPILE_ARGS,
              extra_link_args=LINK_ARGS,
              libraries=LIBS),
]

extensions.extend(
    [Extension("fwdpy.fitness",
               sources=["fwdpy/fitness"+EXTENSION],
               language="c++",
               include_dirs=GLOBAL_INCLUDES,
               extra_compile_args=GLOBAL_COMPILE_ARGS,)]
    )

extensions.extend(
    [Extension("fwdpy.matrix",
            sources=["fwdpy/matrix"+EXTENSION],
            language="c++",
            include_dirs=GLOBAL_INCLUDES,
            extra_compile_args=GLOBAL_COMPILE_ARGS,)]
    )

extensions.extend(
    [Extension("fwdpy.numeric_gsl",
            sources=["fwdpy/numeric_gsl"+EXTENSION],
            language="c++",
            include_dirs=GLOBAL_INCLUDES,
            extra_compile_args=GLOBAL_COMPILE_ARGS,
            libraries=['gsl','gslcblas'],)
            ]
    )

##This is the list of extension modules that are 
##in distinct sub-directories with their own __init__.py
PKGS=['fwdpy','fwdpy.internal','fwdpy.fwdpyio','fwdpy.demography']

if QTRAIT is True:
    extensions.extend([
        Extension("fwdpy.qtrait.qtrait",
                  sources=["fwdpy/qtrait/qtrait"+EXTENSION]+glob.glob("fwdpy/qtrait/*.cc"), # the Cython source and additional C++ source files
                  language="c++",                        # generate and compile C++ code
                  include_dirs=GLOBAL_INCLUDES,
                  extra_compile_args=GLOBAL_COMPILE_ARGS,
                  extra_link_args=LINK_ARGS,
                  libraries=LIBS)
    ])
    PKGS.append('fwdpy.qtrait')
    extensions.extend([
        Extension("fwdpy.qtrait_mloc.qtrait_mloc",
                  sources=["fwdpy/qtrait_mloc/qtrait_mloc"+EXTENSION]+glob.glob("fwdpy/qtrait_mloc/*.cc"),
                  language="c++",
                  include_dirs=GLOBAL_INCLUDES,
                  extra_compile_args=GLOBAL_COMPILE_ARGS,
                  extra_link_args=LINK_ARGS,
                  libraries=LIBS)
    ])
    PKGS.append('fwdpy.qtrait_mloc')
    PKGS.append('fwdpy.gwas')
    extensions.extend([
        Extension("fwdpy.gwas.gwas",
                  sources=["fwdpy/gwas/gwas"+EXTENSION]+glob.glob("fwdpy/gwas/*.cc"),
                  language="c++",
                  include_dirs=GLOBAL_INCLUDES,
                  extra_compile_args=GLOBAL_COMPILE_ARGS,
                  extra_link_args=LINK_ARGS,
                  libraries=LIBS)
    ])


#If using Cython, edit extensions here:
if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

generated_package_data={}
for root, dirnames, filenames in os.walk('fwdpy'):
    if 'fwdpp' not in root:
        g=glob.glob(root+'/*.pxd')
        if len(g)>0:
            replace=root.replace('/','.')
            generated_package_data[replace]=['*.pxd']
    g=glob.glob(root+'/*.hpp')
    if len(g)>0:
        replace=root.replace('/','.')
        PKGS.append(replace) #If there's a header file, we add the directory as a package
        generated_package_data[replace]=['*.hpp']
    
#add in the fwdpp headers, excluding directories other
#than the library itself
for root, dirnames, filenames in os.walk('fwdpy/headers/fwdpp/fwdpp'):
    g=glob.glob(root+'/*.hh')
    if len(g)>0:
        replace=root.replace('/','.')
        if replace not in PKGS:
            PKGS.append(replace) #If there's a header file, we add the directory as a package
        generated_package_data[replace]=['*.hh']
    g=glob.glob(root+'/*.hh')
    if len(g)>0:
        replace=root.replace('/','.')
        if replace not in PKGS:
            PKGS.append(replace) #If there's a header file, we add the directory as a package
        try:
            if '*.hpp' not in generated_package_data[replace]:
                generated_package_data[replace].append('*.hpp')
        except:
            generated_package_data[replace]=['*.hpp']
    g=glob.glob(root+'/*.tcc')
    if len(g)>0:
        replace=root.replace('/','.')
        if replace not in PKGS:
            PKGS.append(replace) #If there's a template implementation file, we add the directory as a package
        try:
            if '*.tcc' not in generated_package_data[replace]:
                generated_package_data[replace].append('*.tcc')
        except:
            generated_package_data[replace]=['*.tcc']

setup(name='fwdpy',
      version='0.0.4-rc3',
      author='Kevin R. Thornton',
      author_email='krthornt@uci.edu',
      maintainer='Kevin R. Thornton',
      maintainer_email='krthornt@uci.edu',
      url='http://www.molpopgen.org',
      description="Forward-time population genetic simulation in Python",
      long_description=long_desc,
      download_url='http://github.com/molpopgen/fwdpy',
      classifiers=['Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'],
      platforms=['Linux','OS X'],
      license='GPL >= 3',
      requires=['pandas','numpy'],
      provides=['fwdpy'],
      obsoletes=['none'],
      packages=PKGS,
      py_modules=[],
      scripts=[],
      data_files=[('fwdpy',['COPYING', 'README.rst'])],
      ##Note: when installing the git repo, headers will be put somewhere like /usr/local/include/pythonVERSION/fwdpy
      #headers=glob.glob("include/*.hpp"),
      package_data=generated_package_data,
      ext_modules=extensions,
      )
