Calculating summary statistics
===========================

Users will typically want to calculate summary statistics based on samples drawn from simulated populations.  There are three ways to accomplish that:

* Write your own functions.
* Use the fwdpy.libseq module in this package, which has some minimal support for summary statistics.
* Use the Python package pyseq_, which provides an interface to the libsequence_ C++11 library.

The pyseq_ package provides a complete API for processing samples.  As pyseq_ evolves, the fwdpy.libseq module may disappear from this package without warning.

.. _libsequence: http://molpopgen.github.io/libsequence/
.. _pyseq: http://molpopgen.github.io/pyseq
