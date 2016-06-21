Overview
=========================

There are places where you can change the behavior of *fwdpy*.  Doing this involves writing some Cython_ code yourself.

In order to do this, you'll have to make your own Python package using Cython_.  This is not that difficult, although it may seem daunting at first.  The Cython manual_ is a very thorough reference, but we'll try in these sections to give you most of the tools you need.

Basic protocol
------------------------------

In this section, we'll refer to custom extensions to fwdpy as "widgets".  A custom fitness function would be an example of such a "widget".

The generic recipe is:

* Generate a "pxd" file that *declares* the name(s) of your widget(s).  If you are familiar with C/C++, this is a lot like a header file providing a list of names of types and functions.
* Generate a "pyx" file that *defines* your widgets.  Again C/C++ users can view this as a source file that contains the code that actually makes these types and functions work.
* Generate a script to compile your fancy custom module!

The next sections of the manual cover the following topics:

* :ref:`customFitness`

.. _Cython: http://cython.org
.. _manual: http://docs.cython.org
