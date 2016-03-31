fwdpy: forward simulation in Python using fwdpp
*****************************************************

This is the main README for the fwdpy software.

This package is a testing ground for providing access to efficient forward-time population simulation machinery in Python.

This package is implemented in terms of:

1. Cython_, which is a package allowing C++ and Python to work together
2. fwdpp_, which is a C++11 template library for implementing efficient population genetic simulations

Please note that this package is likely to be quite unstable/actively developed.

The package is usable now (in fact, we are currently using it for research), but the API, etc., could change without warning.

Citation
===========

See the project home page for details
(http://molpopgen.github.io/fwdpy).

Changelog (rough)
=====================

0.0.3
-----------------
* Change from std::thread to std::async for concurrency.
* The asynchronous futures allow for the same "evolve" function to be
  used in different contexts.
* The different contexts include calculating things from the
  population every "k" generation or doing nothing.
* These things are implemented as classes with call operators and a
  minimal set of API requirements.
* Fixed a bug in "mutation views"
* Better parameter checking for various "evolve" functions
* Source code re-organized so that all header files are installed

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
* fwdpp_

The configure script will enforce minimum version numbers of these dependencies, if necessary.

**Note:** fwdpy may require the 'dev' branch of fwdpp.  The configure script checks for *both* the correct dependency version number *and* specific header files within each depdency.  If the version number check passes, but a subsequent header check fails, then that is a sign that you need a development version of the relevant dependency.  The reason for this situation is that the development of fwdpy has generated ideas for how to make fwdpp more accessible.  This situation will remain until fwdpy stabilizes.

You also need a C++11-compliant compiler.  For linux users, GCC 4.8 or
newer should suffice.  OS X users must use the clang-omp package from brew_.

You may use one or the other of these libraries, but not both.  See the Performance subsection of the Installation section below for how to use these libraries.

Notes for OS X users
---------------------------------

Apple is making life difficult for OS X users.  The release of El Capitan made installing third-party Unix tools into /usr/local more difficult.  A lot of the instructions below ask you to use brew_ to install depdendencies.  Please make sure that you have a working brew_ setup before trying any of the below.  If your setup is not working, please do research online about fixing it, which is beyond the scope of this document.

OS X users are recommended to use brew_ to install the various dependencies:

.. code-block:: bash

   $ brew install clang-omp
   $ brew install gsl
   $ ##Risky:
   $ brew install fwdpp
   $ brew install google-perftools

**Important**: you need to install clang-omp on OS X!  This package
uses openmp for parallelizing some tasks.  Sadly, OS X's compiler does
not come with openmp support, and so you need a third-party compiler
that does.

For brew users, you may or may not have luck with their version of fwdpp.  That package can change rapidly, and thus the brew version may get out-of-sync with the version required for this package.

The required Python package dependencies are in the requirements.txt file that comes with the source.

Anaconda (and OS X, again...)
------------------------------------

Users have run into issues getting fwdpy working with Anaconda-based Python installations.  In fact, I've been unable to get the package to compile on OS X using Anaconda.  I recommend that OS X users use Python3 installed bia Homebrew in lieu of Anaconda.


What Python version?
==================================

I'm developing the package using Python 2.7.6 on an Ubuntu machine.

Currently, the package is not 100% compatible with Python 3.  The goal is to make it work, though.

Installation
==============

The latest release of the package is available via PyPi_, and can be installed with your favorite Python package manager:

.. code-block:: bash

   $ pip install --upgrade fwdpy

OS X users must first install clang-omp from brew_ and use the
following command:

.. code-block:: bash

   $ CC=clang-omp CXX=clang-omp++ pip install fwdpy

Installation from source
----------------------------------------

This section describes "vanilla" installation using the minimal dependencies.

First, install the dependencies (see above).

**Special instructions for OS X users**

All compiler commands below must be prefixed with:

.. code-block:: bash

   $ CC=clang-omp CXX=clang-omp++

This is currently necessary on OS X in order to use a version of clang that supports OpenMP protocols.

Generic instructions:

To install system-wide:

.. code-block:: bash
		
   $ sudo python setup.py install

To install for your user:

.. code-block:: bash

   $ python setup.py install --prefix=$HOME

To uninstall:

.. code-block:: bash

   $ #use 'sudo' here if it is installed system-wide...
   $ pip uninstall fwdpy

To build the package in place and run the unit tets:

.. code-block:: bash

   $ #build package locally:
   $ python setup.py build_ext -i
   $ #run the unit tests:
   $ python -m unittest discover fwdpy/tests

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

   $ CPPFLAGS="-I$HOME/include" LDFLAGS="-L$HOME/lib" pip install fwdpy


Note for developers
=================================

Cython is a static compiler.  Code written in Cython is compiled into C or, in the case of this package, C++.  Finally, the system's C/C++ compiler is used to compile the final Python module.

In order to modify the package, you will need Cython installed:

.. code-block:: bash

   $ pip install Cython

You need Cython >= 0.22.2, so upgrade if you need to:

.. code-block:: bash

   $ pip install --upgrade Cython


If you wish to modify the package, then you will want setup.py to "re-Cythonize" when you make changes to the package source code.

To do this, use the setup.py script as follows:

.. code-block:: bash

   $ python setup.py build_ext -i --use-cython

Now, Cython will be a compilation depdendency, and any changes to .pyx/.pyd/.cc files in this package will trigger Cython to regenerate the .cpp files that make up the core of the package.

Compiling in an aggressive debug mode
-----------------------------------------------

To get rid of optimizations, and -DNDEBUG, you need to reset the OPT
flag set by Python's distutils:

.. code-block:: bash

   $ OPT= python setup.py build_ext -i

Doing this will mean that the fwdpp back-end will *not* be compiled
with -DNDEBUG, which will enable aggressive run-time correctness
testing.  By "aggressive", I mean that an error will trigger a failed
assertion and the Python interpreter will be exited
less-than-gracefully!  Only to this when testing.

It is better to enable some optimizations, though, else things run too
slowly:

.. code-block:: bash

   $ OPT=-O2 python setup.py build_ext -i
   

Rough guide to installation on UCI HPC
-----------------------------------------

Use the following module:

.. code-block:: bash

   $ module load krthornt/thorntonlab

That command loads the proper dependencies for compiling much of the tools that we use.

**Note**: this module replaces/over-rules some modules already on HPC.  The "thorntonlab" modules are all consistently compiled with a GCC version that we've deemed suitable.

Troubleshooting the installation
-----------------------------------------

Incorrect fwdpp version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This package is compatible with fwdpp >= 0.4.7, which means that you should have a binary installed on your systems called fwdppConfig.  You can check if you have it:

.. code-block:: bash

   $ which fwdppConfig


If the above command returns nothing, then it is very likely that fwdpp is either too old, missing entirely from your system, or it is installed somewhere non-standard.  For example, if you installed fwdpp locally for your user, and did not edit PATH to include ~/bin, then fwdppConfig cannot be called without referring to its complete path.

Dependencies in non-standard locations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Your system's compiler has a default set of paths where it will look for header files, libraries, etc.  Typically, these paths will include /usr and /usr/local.  If you have installed the dependencies somewhere else (your home directory, for example), then the ./configure script may not be able to find them automatically.

**NOTE:** I sometimes get requests for installation help from users who have installed every dependency in a separate folder in their $HOME.  In other words, they have some setup that looks like this:


* $HOME/software/gsl
* $HOME/software/fwdpp


If you insist on doing this, then you are on your own.  You have to manually pass in all of the -I and -L flags to all of these locations.   This setup is problematic because it violates the POSIX Filesystem Hierarchy Standard (http://en.wikipedia.org/wiki/Filesystem_Hierarchy_Standard), and you cannot reasonably expect things to "just work" any more.  It would be best to start over, and simply install all of the dependencies into the following prefix:

.. code-block:: bash

   $ $HOME/software

Doing so will allow $HOME/software/include, etc., to be populated as they were intended to be.

Documentation
===================

The manual_ is available online in html format at the project web page.

The API documentation may also be build using doxygen_:

.. code-block:: bash

   $ ./configure
   $ doxygen fwdpy.doxygen

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
