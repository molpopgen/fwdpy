Changelog (rough)
=====================

0.0.4-rc3
-------------------

* fwdpy_ is now included with fwdpy!  Further, the fwdpp_ headers are installed as package data, meaning that plugins/extensions to fwdpy always have access to the version of fwdpp_ used to build fwdpy.
* :func:`fwdpy.fwdpy.get_includes` and :func:`fwdpy.fwdpy.get_fwdpp_includes` return the locations of the instaled fwdpy and fwdpp_ headers, respectively.
* :func:`fwdpy.fdwpy.tidy_trajectories` removed.
* C++ back-end and Cython_ class definition of :class:`fwdpy.fwdpy.FreqSampler` refactored. New version is much, much faster!
* :class:`fwdpy.fwdpy.FreqSampler` is now able to output directly to SQLite database files.  There is also a new member function called "fetch" that allows filtering of trajectories before returning them as a Pandas DataFrame object.
* fwdpy.numeric_gsl added, providing a Cython_ (nogil) API to some numeric operations implemented in terms of the GSL 
  
0.0.4 (through release candidate 2)
---------------------------------------

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
* Streamline implementation of the various :class:`fwdpy.fwdpy.PopVec` classes.  They no longer contain two containers,
  and they yield :class:`fwdpy.fwdpy.PopType` objects upon iteration.

Changes to the Cython_/C++ back end:

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
  in fwdpp 0.4.9 that made the Cython_ wrappers in this package
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

.. _fwdpp: http://molpopgen.github.io/fwdpp
.. _Cython: http://www.cython.org/
.. _cythonGSL: https://pypi.python.org/pypi/CythonGSL
