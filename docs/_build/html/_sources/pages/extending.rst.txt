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

fwdpy installs its own headers and the fwdpp_ headers that it was built with.  The locations of these headers are given by:

.. code-block:: python

    import fwdpy
    fwdpy_headers=fwdpy.get_includes()
    fwdpp_headers=fwdpy.get_fwdpp_includes()

API Documentation
---------------------------------------------------------------------

Cython's .pxd files have no analog to Python's dosctrings.  However, those pxd files do tell you the name of the header where each C++ declaration can be found.  Those headers are (starting to be) documented into a reference manual using doxygen_.  Currently, this is a work in progress.

To build the API documentation:

.. code-block:: bash

   ./configure
   doxygen fwdpy.doxygen

Then, point your browser to html/index.html.

.. _fwdpp: http://molpopgen.github.io/fwdpp
.. _doxygen: http://doxygen.org

   


   
