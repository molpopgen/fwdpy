fwdpy: forward simulation in Python using fwdpp
*****************************************************

This is the main README for the fwdpy software.

This package is a testing ground for providing access to efficient forward-time population simulation machinery in Python.

This package is implemented in terms of:

1. Cython_, which is a package allowing C++ and Python to work together
2. fwdpp_, which is a C++11 template library for implementing efficient population genetic simulations
3. libsequence_, which is a C++11 library for various population-genetic calculations.

Please note that this package is likely to be quite unstable/actively developed.  The first quasi-stable, but usable, release is 0.0.4.

Previous 'releases' on Github are no longer supported and not guaranteed to even compile.

Citation
===========

The digital object identifier (DOI) for this software is available for each release.  You can find it from the [repo page](http://github.com/molpopgen/fwdpy) for this project.

Please cite it if you use this package for your research.

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
* Parallel executiom of simulations.  Multiple replicates may be run simultaenously via C++11's threading mechanism.  This is a "sneaky" around Pythons' Global Interpreter Lock, or GIL.

The following distributions of selection coefficients are supported:

* constant (*i.e.*, *s* takes on a fixed value)
* beta
* exponential
* gamma
* gaussian

The following distributions of dominance are supported:

* constant (*i.e.*, *h* takes on a fixed value)
* beta
* uniform

The distributions on *s* and *h* can be mixed and matched as the user sees fit, and these parameters can vary along a simulated region.

Example
=============
Let's simulate the following model:

* The population will consist of N=1000 diploids
* The total mutation rate to variants not affecting fitness will be theta = 4Nu = 100 for the region.  Neutral mutations occur uniformly on the interval [0,1)
* The recombination rate will be rho = 4Nr = 1000 for the region
* The *total* mutation rate (per gamete, per generation) to mutations affecting fitness will be 0.01.

Given that a mutation affecting fitness mutation arises:

* There is a 99 percent chance that it is deleterious with sh = -0.1.  These mutations occur uniformly on the interval of [-1,0) or [1,2) with equal probability.
* There is a 1 percent chance that it is beneficial with s=0.001 and h = 1, and occur uniformly along the interval [-1,2).

This model of selected mutations is basically a standard model of background selection with a low rate of weakly-beneficial mutations thrown in.

Recombination is uniform along the interval [-2,2)

The demographic model is the bottleneck for European *Drosophila* inferred by Thornton, K., & Andolfatto, P. (2006). Approximate Bayesian inference reveals evidence for a recent, severe bottleneck in a Netherlands population of Drosophila melanogaster. Genetics, 172(3), 1607-1619. http://doi.org/10.1534/genetics.105.048223

Let's set up the model and run it:

.. code-block:: python

   #not done yet

They use GNU parallel to run jobs in parallel and impose a wall clock limit.

Availability
===============

This package is distributed at the following github repository: https://github.com/molpopgen/fwdpy.

Dependencies
===============

This section assumes that all packages are installed in fairly standard locations, such as /usr/local.  See the troubleshooting section for more complex setups.

A lot of them:

* Cython_, preferably something like 0.22 or more recent.  (I'm not really sure what the minimum version requirement here is.)
* GSL_
* fwdpp_ 
* libsequence_
* tcmalloc_

The configure script will enforce minimum version numbers of these dependencies, if necessary.

You also need a C++11-compliant compiler.  For OS X users, that means Yosemite + current Xcode installation.  For linux users, GCC 4.8 or newer should suffice.

OS X users are recommended to use brew_ to install the various dependencies:

.. code-block:: bash
   
   $ brew install gsl
   $ brew install libsequence
   $ ##Risky:
   $ brew install fwdpp


For brew users, you may or may not have luck with their version of fwdpp.  That package can change rapidly, and thus the brew version may get out-of-sync with the version required for this package.

The required Python package dependencies are in the requirements.txt file that comes with the source.  

Installation
==============

First, install the dependencies (see above).

To install system-wide

.. code-block:: bash

   $ ./configure
   $ sudo python setup.py install

To install for your user:

.. code-block:: bash

   $ ./configure --prefix=$HOME
   $ #yes, the prefix is needed again here...
   $ python setup.py install --prefix=$HOME

To uninstall:

.. code-block:: bash

   $ #use 'sudo' here if it is installed system-wide...
   $ pip uninstall fwdpy


Rough guide to installation on UCI HPC
-----------------------------------------


Troubleshooting the installation
-----------------------------------------

Incorrect fwdpp version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This package is compatible with fwdpp >= 0.3.3, which means that you should have a binary installed on your systems called fwdppConfig.  You can check if you have it:

.. code-block:: bash

   $ which fwdppConfig


If the above command returns nothing, then it is very likely that fwdpp is either too old, missing entirely from your system, or it is installed somewhere non-standard.  For example, if you installed fwdpp locally for your user, and did not edit PATH to include ~/bin, then fwdppConfig cannot be called without referring to its complete path.

Dependencies in non-standard locations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Your system's compiler has a default set of paths where it will look for header files, libraries, etc.  Typically, these paths will include /usr and /usr/local.  If you have installed the dependencies somewhere else (your home directory, for example), then the ./configure script may not be able to find them automatically.

**NOTE:** I sometimes get requests for installation help from users who have installed every dependency in a separate folder in their $HOME.  In other words, they have some setup that looks like this:


* $HOME/software/boost
* $HOME/software/libsequence
* $HOME/software/gsl
* $HOME/software/fwdpp


If you insist on doing this, then you are on your own.  You have to manually pass in all of the -I and -L flags to all of these locations.   This setup is problematic because it violates the POSIX [ilesystem Hierarchy Standard (http://en.wikipedia.org/wiki/Filesystem_Hierarchy_Standard), and you cannot reasonably expect things to "just work" any more.  It would be best to start over, and simply install all of the dependencies into the following prefix:

.. code-block:: bash

   $ $HOME/software

Doing so will allow $HOME/software/include, etc., to be populated as they were intended to be.

Documentation
===================

TBD


.. _fwdpp: http://molpopgen.github.io/fwdpp 
.. _libsequence: http://molpopgen.github.io/libsequence/
.. _Cython: http://www.cython.org/
.. _GSL:  http://gnu.org/software/gsl
.. _tcmalloc: https://code.google.com/p/gperftools/
.. _brew: http://brew.sh
