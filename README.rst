fwdpy: forward simulation in Python using fwdpp
*****************************************************

This is the main README for the fwdpy software.

This package is a testing ground for providing access to efficient forward-time population simulation machinery in Python.

This package is implemented in terms of:

1. Cython_, which is a package allowing C++ and Python to work together
2. fwdpp_, which is a C++11 template library for implementing efficient population genetic simulations
3. libsequence_, which is a C++11 library for population genetic calculations.
4. gsl_, which is a C-language library for numerical methods.  This package uses the GSL random number generation system plus several other features.

..note:: libsequence_ is only used internally for some format conversions and other operations.

Please note that this package is likely to be quite unstable/actively developed.

The package is usable now (in fact, we are currently using it for research), but the API, etc., could change without warning.

Citation
===========

See the project home page for details
(http://molpopgen.github.io/fwdpy).

TODO 
=================

* Document custom temporal samplers.  Target version 0.0.4
* Add safety check to (de)serialization.  Current method trusts the input string completely.  This likely requires a
  saftey-checking template that is specialized for each population type.  Target version 0.0.5.
* Remove fwdpy.fwdpyio module in favor of object-level serialization.  Target version 0.0.5.
* Serialization should be at the level of PopType and PopVec objects.  The latter can be done in parallel.  Target
  version 0.0.5.
* Serialization should support direct to/from file, with gzip as option via zlib. Target version 0.0.5.
* Custom rules classes. This will allow "stateful" fitness models, and many other things.  Target version 0.0.5.

Changelog (rough)
=====================

0.0.4
----------------

Changes to the Python side:

* "Evolve" functions are now much more generic due to fitness function
  objects and temporal sampler objects (see below)
* The application of temporal samplers is now consistent for all types of simulation ("pop-gen", "quant-trait", etc.)
* Added ability to use custom fitness functions!!! On the Python side,
  these work via :class:`fwdpy.fwdpy.SpopFitness` and :class:`fwdpy.fwdpy.MlocusFitness`
* Class names now more "Pythonic".  This will break existing scripts.
* Add fwdpy.demography module.
* Add :class:`fwdpy.fwdpy.MlocusPop`
* Add :class:`fwdpy.fwdpy.MlocusPopVec`
* Add concept of a temporal sampler via
  :class:`fwdpy.fwdpy.TemporalSampler`.
* Add temporal sampler objects :class:`fwdpy.fwdpy.NothingSampler`,
  :class:`fwdpy.fwdpy.QtraitStatsSampler`,
  :class:`fwdpy.fwdpy.PopSampler`,
  :class:`fwdpy.fwdpy.VASampler`,
  :class:`fwdpy.fwdpy.FreqSampler`
* Add function :func:`fwdpy.fwdpy.apply_sampler`
* Add :func:`fwdpy.fwdpy.tidy_trajectories`, which really speeds up
  coercion of mutation frequency trajectories to a pandas DataFrame.
* Add :func:`fwdpy.fwdpy.hapmatrix` and :func:`fwdpy.fwdpy.genomatrix`
* Added views of fixed mutations via :func:`fwdpy.fwdpy.view_fixations`
* Better Python3 compatibility
* Add support to serialize/deserialize :class:`fwdpy.fwdpy.MlocusPop`

Changes to the Cython/C++ back end:

* diploid fitness now defaults to 1 instead of 0
* Bug fixed in get_gamete in views.pyx.  This affected the output of almost all "views" functions except those viewing just mutations.
* cythonGSL_ is now required. We expect to use more GSL in this package, and so it makes sense to not reinvent the wheel.
* Massive reduction in code base
* Update to Cython_ 0.24.0
* Generic temporal samplers and fitness functions are now supported.
* Expose more fwdpp types for multi-locus/region simulations
* Expose fwdpp's fitness function objects site_dependent_fitness,
  additive_diploid, and multiplicative_diploid.  Call operators
  (e.g. operator()) are only exposed for custom diploids.
* More unit tests of sampling and "views"
* Update how samples are taken from populations, reflecting a bug fix
  in fwdpp 0.4.9 that made the Cython wrappers in this package
  incorrect.
* Population objects in types.hpp now have serialization/deserialization functions.
* Single-parameter constructors for population objects in types.hpp are now "explicit".

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

.. note:: If installing from GitHub, then you also must have Cython_ >= 0.24.0 and cythonGSL_ installed on your system.

.. note:: fwdpy may require the 'dev' branch of fwdpp.  The configure script checks for *both* the correct dependency version number *and* specific header files within each depdency.  If the version number check passes, but a subsequent header check fails, then that is a sign that you need a development version of the relevant dependency.  The reason for this situation is that the development of fwdpy has generated ideas for how to make fwdpp more accessible.  This situation will remain until fwdpy stabilizes.

You also need a C++11-compliant compiler.  For linux users, GCC 4.8 or
newer should suffice.  OS X users must use the clang-omp package from brew_, or gcc6 from brew_ if they use Anaconda.

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

**Important**: you need to install clang-omp on OS X!  This package
uses openmp for parallelizing some tasks.  Sadly, OS X's compiler does
not come with openmp support, and so you need a third-party compiler
that does.

For brew users, you may or may not have luck with their version of fwdpp_.  That package can change rapidly, and thus the brew version may get out-of-sync with the version required for this package.

The required Python package dependencies are in the requirements.txt file that comes with the source.

Anaconda (and OS X, again...)
------------------------------------

**As things stabilize, the dependencies, etc., will be installable from Bioconda.**

Anaconda_ allows user to manage their own Python installations for themselves. It is expecially useful on systems where you don't have root access.  Anaconda_ works by installing intself into $HOME/anaconda2 (for Python2) or $HOME/anaconda3 (for Python3).  Further, it works best if all dependencies are also installed in the same location.

An advantage of Anaconda_ is that you can avoid some of the complications involved in managing dependencies. By using Anaconda_'s installation of the GCC compiler, you can can guarantee that dependencies are compatible with one another (*i.e.* no "ABI compatibility" issues).  However, this means you must manage the installation of the dependencies yourself.  Here, I outline how to do, based on what worked for me on my Ubuntu 16.04 system.

.. note:: The "recipes" below were all tested on new user accounts, thus avoiding any complications due to settings in my main accounts.

Ubuntu
=====================

I am assuming that Anaconda_ for Python2 is installed. I'm further assuming that this all works for their Python3 installation, as fwdpy is compatible with both major versions of Python. Finally, I assume that $HOME/anaconda2/bin is prepended to your $PATH, meaning that Anaconda_ binaries are preferred over system binaries.

First, install Anaconda_'s GCC, GSL_, and zlib:

.. code-block:: bash

    $ conda install gcc
    $ conda install gsl
    $ conda install zlib

Now, we wish to install libsequence_, which depends on Intel's TBB library, which we will install as follows:

.. code-block:: bash

    $ conda install -c dlr-sc tbb=4.3.6

Get libsequence_'s dev branch and install:

.. code-block:: bash

    $ git clone http://github.com/molpopgen/libsequence
    $ cd libsequence
    $ git branch dev
    $ ./configure --prefix=$HOME/anaconda2
    $ #change -j to some number of threads
    $ #appropriate for your machine
    $ make -j 40 && make install

Get fwdpp_'s dev branch and install:

.. code-block:: bash

    $ git clone http://github.com/molpopgen/fwdpp
    $ cd fwdpp
    $ git branch dev
    $ ./configure --prefix=$HOME/anaconda2
    $ make && make install

Install cythonGSL_, which is a dependency for fwdpy:

.. code-block:: bash

    $ pip install cythongsl

Get the dev branch of fwdpy and install:

.. code-block:: bash

    $ pip install git+git://github.com/molpopgen/fwdpy@dev --install-option="--use-cython"

The result of all of the above is:

* All dependencies are compiled with the same version of GCC, which is whatever Anaconda_ is currently using (GCC 4.8.5 at the time of this writing).
* All dependencies get installed into $HOME/anaconda2
* fwdpy is installed and linked against the dependencies in $HOME/anaconda2

OS X
=====================

As is too often the case, the situation on OS X is more complex.  If you want to use Anaconda_ on OS X, then the following worked for me.  I did it using a Python3 installation this time, just for fun.

For OS X, we will rely on installing the Anaconda_ GCC.  It is also possible to use GCC6 from brew_, but I will not document that here, and instead focus on the path of least resistance, which is an "all 'conda" approach.

.. note:: Installing Anaconda_ GCC means that compiler will be preferrred over the Xcode installation of clang, which is aliased to GCC on OS X.  Thus, there may be side effects when you play outside the Anaconda world.

All of the steps shown above for Ubuntu work, with the following modifications:

Any steps involving a "./configure" command need to have $HOME/anaconda[2|3]/include added to CPPFLAGS:

.. code-block:: bash

    $ CPPFLAGS="-I$HOME/anaconda3/include" ./configure --prefix=$HOME/anaconda3

(Use anaconda2 instead of anaconda3 as needed.)

The command to install fwdpy from GitHub must be told which compilers to use.  No idea why, but Anaconda_ on OS X really likes to force the use of clang!

.. code-block:: bash

    $ CC=gcc CXX=g++ pip install git+git://github.com/molpopgen/fwdpy@dev --install-option="--use-cython"

Finally, it is wise for OS X users to add the following to their .bash_profiles:

.. code_block:: bash

    $ LD_LIBRARY_PATH=$HOME/anaconda3/lib
    $ export LD_LIBRARY_PATH

Again, substitute anaconda2 as necessary.  For the record, no idea why this is needed on OS X but not Linux...

What Python version?
==================================

I'm developing the package using Python 2.7.6 on an Ubuntu machine.  However, I do occasionally run the tests using Python 3, and all appears to work!  Reports of problems using python3 are appreciated!

Installation
==============

The latest release of the package is available via PyPi_, and can be installed with your favorite Python package manager:

.. code-block:: bash

   $ pip install --upgrade fwdpy

OS X users must first install clang-omp from brew_ and use the
following command:

.. code-block:: bash

   $ CC=clang-omp CXX=clang-omp++ pip install fwdpy

Installation from GitHub
----------------------------------------

You may also use pip to install from GitHub.  However, doing so requires that Cython_ be installed.

.. code-block:: bash

   $ pip install git+git://github.com/molpopgen/fwdpy --install-option="--use-cython"

The above command installs the latest version of the 'master' branch.  Users wanting latest and buggiest may find this useful.  OS X users should follow the instructions for using clang-omp shown above.

Do this at your own risk. While the version number of the master branch may be the same as the version on PyPi_, there may be bugs, API changes, etc.

To install a specific branch:

   $ pip install git+git://github.com/molpopgen/fwdpy@branchname --install-option="--use-cython"

Installation from source
----------------------------------------

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

To build the package in place and run the unit tests:

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

Testing
======================================

Testing occurs via docstring tests and unit tests.  Here is how to test using both methods:

.. code-block:: bash

   $ #build the package
   $ python setup.py build_ext -i
   $ #build the manual--requires Sphinx
   $ make -f Makefile.sphinx html
   $ #run the tests
   $ make -f Makefile.sphinx doctest
   $ #run the unit tests
   # python -m unittest discover fwdpy/tests
   

Note for developers
=================================

Cython is a static compiler.  Code written in Cython is compiled into C or, in the case of this package, C++.  Finally, the system's C/C++ compiler is used to compile the final Python module.

In order to modify the package, you will need Cython installed:

.. code-block:: bash

   $ pip install Cython

You need Cython >= 0.24.0, so upgrade if you need to:

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

This package is compatible with fwdpp >= 0.4.8, which means that you should have a binary installed on your systems called fwdppConfig.  You can check if you have it:

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

The manual_ is available online in html format at the project web page.  The manual always corresponds to the version of *fwdpy* found on PyPi_.

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
.. _cythonGSL: https://pypi.python.org/pypi/CythonGSL
.. _libsequence: http://molpopgen.github.io/libsequence
.. _Anaconda: https://www.continuum.io/why-anaconda
