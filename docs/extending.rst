Extending *fwdpy*
=======================================

It is possible to use *fwdpy* as the basis for other Python/Cython packages.  Once installed, you have access to both Cython and C++ declarations to use in your own code.  The Cython definitions are found in the "pxd" files that come with the source code and are installed with the package.  You use these in your own package via Cython's "cimport" mechanism.  See the Cython documentation for details.

The "pxd" files give you access to C++ types (classes, structs, etc.), Cythons "extension types" (very similar to Python's classes, but with the ability to have function implemented in another language), and functions.  With these files, you can write your own Cython code, should you ever find yourself needing direct access to the details of the C++ types.

The C++ headers give you direct access to different types of C++ constructs:

1. The raw types exposed by the "pxd" files.
2. C++ types and functions that cannot be exposed to Python via Cython.  Cython has limited support for "modern" C++11/C++14.
3. Low-level implementation details.

You may use *fwdpy*'s headers directly in your own C++ code in order to write new routines for your own Python extension modules.

Including the C++ headers
-------------------------------------------

In order to compile any extension based on *fwdpy* will require knowledge of where your Python installation installs header files that come with extension modules.  **In general, this will be the trickiest part of making a module based on fwdpy.**  The section is quite detailed, in order to help you out.

In general, there will be some ROOT directory where Python puts things.  What that ROOT is depend on how you installed a package.

For example, when I use python3 to install the package on my OS X machine using the following command:

.. code-block:: bash
		
   $ #We must use an OMP-enabled clang on OS X:
   $ CC=clang-omp CXX=clang-omp++ python3 setup.py install

I see this in the output:

.. code-block:: bash
		
   $ copying fwdpy/allele_ages.hpp -> /usr/local/include/python3.5m/fwdpy

This tells me that I need to add *-I/usr/local/include/python3.5m/fwdpy* when compiling an extension module based on fwdpy.  

Unfortunately, python does not seem to offer any automatic way of finding out where files are installed.  That said, python can give you hints:

.. code-block:: python
		
   >>> import imp
   >>> imp.find_module('fwdpy')

Or:

.. code-block:: python
		
   >>> import fwdpy,os
   >>> os.path.dirname(fwdpy.__file__)

On my OS X machine, the latter gives:

.. code-block:: python
		
   >>> '/usr/local/lib/python3.5/site-packages/fwdpy'

That path refers to where the *module itself* is installed.  The ROOT is /usr/local, you should be able to find the path to the headers with the help of our good friend bash:

.. code-block:: bash
		
   $ find /usr/local -name "*.hpp" | grep fwdpy

   


   
