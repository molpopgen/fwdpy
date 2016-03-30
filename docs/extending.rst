Extending *fwdpy*
=======================================

It is possible to use *fwdpy* as the basis for other Python/Cython packages.  Once installed, you have access to both Cython and C++ declarations to use in your own code.  The Cython definitions are found in the "pxd" files that come with the source code and are installed with the package.  You use these in your own package via Cython's "cimport" mechanism.  See the Cython documentation for details.


