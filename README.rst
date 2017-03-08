fwdpy: forward simulation in Python using fwdpp
*****************************************************

This is the main README for the fwdpy software.

This package is a testing ground for providing access to efficient forward-time population simulation machinery in Python.

This package is implemented in terms of:

1. Cython_, which is a package allowing C++ and Python to work together
2. fwdpp_, which is a C++11 template library for implementing efficient population genetic simulations
3. libsequence_, which is a C++11 library for population genetic calculations.
4. gsl_, which is a C-language library for numerical methods.  This package uses the GSL random number generation system plus several other features.

.. note:: libsequence_ is only used internally for some format conversions and other operations.

Please note that this package is likely to be quite unstable/actively developed.

The package is usable now (in fact, we are currently using it for research), but the API, etc., could change without warning.

Citation
===========

See the project home page for details
(http://molpopgen.github.io/fwdpy).

Features:
===========

So far, there is support for:

* Simulation of a recombining region with arbitrary variation in neutral mutation rate, recombination rate, and distribution of selective effects and their dominance along the region.
* Support for arbitrary changes in population size.
* The ability to sample from simulated populations.
* Calculate some standard summary statistics from samples taken from simulated populations.
* Selfing
* The ability to vary model parameters over time (recombination rates, genetic maps, selfing, selection, etc.)
* Sampling populations at various time points
* Parallel executiom of simulations.  Multiple replicates may be run simultaenously via C++11's threading mechanism.  This is a "sneaky" end-run around Python's Global Interpreter Lock, or GIL.

The following distributions of selection coefficients are supported:

* constant (*i.e.*, *s* takes on a fixed value)
* uniform
* exponential
* gamma
* gaussian

The following distributions of dominance are supported:

* constant (*i.e.*, *h* takes on a fixed value)

Google Group
=================

For help with various issues, check out the `fwdpy Google Group`_

Examples
=============

Several examples exist in the form of output from "iPython"/Jupyter notebooks:

* `Background selection`_
* `Viewing simulated populations`_
* `Sliding windows`_
* `Tracking mutation frequencies`_

Availability
===============

This package is distributed at the following github repository: https://github.com/molpopgen/fwdpy.

Dependencies
===============

This section assumes that all packages are installed in fairly standard locations, such as /usr/local.  See the troubleshooting section for more complex setups.

This package *minimally* depends on:

* GSL_
* libsequence_

.. note:: fwdpp_ is not a build dependency. The library is included with fwdpy and its headers are installed as package data.

.. note:: If installing from GitHub, then you also must have Cython_ >= 0.24.0 and cythonGSL_ installed on your system.

You also need a C++11-compliant compiler.  For linux users, GCC 4.8 or
newer should suffice.  

Notes for OS X users
---------------------------------

OS X users are recommended to use Anaconda_ instead of brew_.  Further, use the Anaconda_ version of GCC instead of the system (Xcode) Clang if you compile fwdpy from source.  This package requires OpenMP, which is not supported via the clang provided with Xcode.

See next section for some details.

Anaconda
------------------------------------

Anaconda_ may be the easiest way to install this package for many users.  The dependencies are available via the Bioconda_ "channel".

Note that using Anaconda_ means over-riding some things that may be provided with your system.  For example, if you install dependencies via Bioconda_ and then wish to install fwdpy from source, you will need GCC from Anaconda_:

.. code-block:: bash

    conda install gcc

The GCC version in Anaconda_ is 4.8.5, which is a bit old but sufficient for the C++11 features needed for all dependencies and for this package.

In order to make sure that the Anaconda_ GCC is used, you will need to make sure that the bin directory of your Anaconda installation is prepended to your users's PATH variable.

If we define CONDAROOT as the location of your Anaconda_ installation, then you should define the following environment variables for your user in the dotfile appropriate for your favorite shell.  For example, for the bash shell:

.. code-block:: bash

    export PATH=$CONDAROOT/bin:$PATH
    export CPPFLAGS="-I$CONDAROOT/include $CPPFLAGS"
    export CFLAGS="-I$CONDAROOT/include $CFLAGS"
    export LDFLAGS=-L$CONDAROOT/lib $LDFLAGS"
    export LD_LIBRARY_PATH="$CONDAROOT/lib:$LD_LIBRARY_PATH"

.. note::

    The above exports *prepend* Anaconda_ paths to existing paths (if they exist).  If you use the system GCC for your own work, then the PATH export may not be something you want set all of the time.

What Python version?
==================================

I'm developing the package using Python 2.7.6 on an Ubuntu machine.  However, I do occasionally run the tests using Python 3, and all appears to work!  Reports of problems using python3 are appreciated!

Installation
==============

The latest release of the package is available via PyPi_, and can be installed with your favorite Python package manager:

.. code-block:: bash

   pip install --upgrade fwdpy

OS X users must first install a compiler that supports the -fopenmp option.  I recommend GCC from Anaconda_ (see above).

Installation from GitHub
----------------------------------------

You may also use pip to install from GitHub.  However, doing so requires that Cython_ be installed.

.. code-block:: bash

   pip install git+git://github.com/molpopgen/fwdpy --install-option="--use-cython"

The above command installs the latest version of the 'master' branch.  Users wanting latest and buggiest may find this useful.  OS X users should follow the instructions for using clang-omp shown above.

Do this at your own risk. While the version number of the master branch may be the same as the version on PyPi_, there may be bugs, API changes, etc.

To install a specific branch:

   pip install git+git://github.com/molpopgen/fwdpy@branchname --install-option="--use-cython"

Installation from source
----------------------------------------

First, install the dependencies (see above).

The best way to install the package is to use 'pip'.  Once you have cloned the source repo and 'cd' into it:

.. code-block:: bash

    git submodule init
    git submodule update
    pip install . --upgrade --intall-option=--use-cython

To build the package in place and run the unit tests:

.. code-block:: bash

    git submodule init
    git submodule update
    #build package locally:
    python setup.py build_ext -i
    #run the unit tests:
    python -m unittest discover fwdpy/tests

Dependencies in non-standard locations
----------------------------------------------------------------------------------------

The instructions above assume that dependencies (fwdpp_ and GSL_) are
found in "standard" locations, which means in /usr/local on a typical
system.

Many users, especially those on clusters, may not have the privileges
needed to install to the standard system locations.  Thus, it may be
necessary to manually tell fwdpy where the dependencies are located.

For example, let us assume that fwdpp_ and GSL_ are installed into
your home folder. On Unix-like systems, $HOME is a variable representing
the location of your home folder.  Thus, the header files for these
libraries will be found in $HOME/include and any run-time libraries
will be found in $HOME/lib.

To tell pip where to find these dependencies, you need to manually set
CPPFLAGS and LDFLAGS:

.. code-block:: bash

   CPPFLAGS="-I$HOME/include" LDFLAGS="-L$HOME/lib" pip install fwdpy

Testing
======================================

Testing occurs via docstring tests and unit tests.  Here is how to test using both methods:

.. code-block:: bash

   #build the package
   python setup.py build_ext -i
   #build the manual--requires Sphinx
   make -f Makefile.sphinx html
   #run the tests
   make -f Makefile.sphinx doctest
   #run the unit tests
   python -m unittest discover fwdpy/tests
   

Note for developers
=================================

Cython is a static compiler.  Code written in Cython is compiled into C or, in the case of this package, C++.  Finally, the system's C/C++ compiler is used to compile the final Python module.

In order to modify the package, you will need Cython installed:

.. code-block:: bash

   pip install Cython

You need Cython >= 0.24.0, so upgrade if you need to:

.. code-block:: bash

   pip install --upgrade Cython

If you wish to modify the package, then you will want setup.py to "re-Cythonize" when you make changes to the package source code.

To do this, use the setup.py script as follows:

.. code-block:: bash

   python setup.py build_ext -i --use-cython

Now, Cython will be a compilation depdendency, and any changes to .pyx/.pyd/.cc files in this package will trigger Cython to regenerate the .cpp files that make up the core of the package.

Compiling in an aggressive debug mode
-----------------------------------------------

To get rid of optimizations, and -DNDEBUG, you need to reset the OPT
flag set by Python's distutils:

.. code-block:: bash

   OPT= python setup.py build_ext -i

Doing this will mean that the fwdpp back-end will *not* be compiled
with -DNDEBUG, which will enable aggressive run-time correctness
testing.  By "aggressive", I mean that an error will trigger a failed
assertion and the Python interpreter will be exited
less-than-gracefully!  Only to this when testing.

It is better to enable some optimizations, though, else things run too
slowly:

.. code-block:: bash

   OPT=-O2 python setup.py build_ext -i
   
Troubleshooting the installation
-----------------------------------------

Dependencies in non-standard locations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Your system's compiler has a default set of paths where it will look for header files, libraries, etc.  Typically, these paths will include /usr and /usr/local.  If you have installed the dependencies somewhere else (your home directory, for example), then the ./configure script may not be able to find them automatically.

**NOTE:** I sometimes get requests for installation help from users who have installed every dependency in a separate folder in their $HOME.  In other words, they have some setup that looks like this:

* $HOME/software/gsl
* $HOME/software/fwdpp

If you insist on doing this, then you are on your own.  You have to manually pass in all of the -I and -L flags to all of these locations.   This setup is problematic because it violates the POSIX Filesystem Hierarchy Standard (http://en.wikipedia.org/wiki/Filesystem_Hierarchy_Standard), and you cannot reasonably expect things to "just work" any more.  It would be best to start over, and simply install all of the dependencies into the following prefix:

.. code-block:: bash

   $HOME/software

Doing so will allow $HOME/software/include, etc., to be populated as they were intended to be.

Better yet, use a system like Anaconda_ (see above).

Documentation
===================

The manual_ is available online in html format at the project web page.  The manual always corresponds to the version of *fwdpy* found on PyPi_.

The API documentation may also be build using doxygen_:

.. code-block:: bash

   ./configure
   doxygen fwdpy.doxygen

Then, load html/index.html in your browser.

.. _fwdpp: http://molpopgen.github.io/fwdpp
.. _Cython: http://www.cython.org/
.. _GSL:  http://gnu.org/software/gsl
.. _brew: http://brew.sh
.. _manual: http://molpopgen.github.io/fwdpy
.. _Background selection: http://molpopgen.github.io/fwdpy/_build/html/examples/BGS.html
.. _Viewing simulated populations: http://molpopgen.github.io/fwdpy/_build/html/examples/views.html
.. _Sliding windows: http://molpopgen.github.io/fwdpy/_build/html/examples/windows.html
.. _Tracking mutation frequencies: http://molpopgen.github.io/fwdpy/_build/html/examples/trajectories.html
.. _PyPi: https://pypi.python.org
.. _fwdpy Google Group: https://groups.google.com/forum/#!forum/fwdpy-users
.. _doxygen: http://doxygen.org
.. _cythonGSL: https://pypi.python.org/pypi/CythonGSL
.. _libsequence: http://molpopgen.github.io/libsequence
.. _Anaconda: https://www.continuum.io/why-anaconda
.. _Bioconda: https://bioconda.github.io
