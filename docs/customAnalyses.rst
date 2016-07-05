Overview
=========================

There are places where you can change the behavior of *fwdpy*.  Doing this involves writing some Cython_ code yourself.

In order to do this, you'll have to make your own Python package using Cython_.  This is not that difficult, although it may seem daunting at first.  The Cython manual_ is a very thorough reference, but we'll try in these sections to give you most of the tools you need.

Basic protocol
------------------------------

In this section, we'll refer to custom extensions to fwdpy as "widgets".  A custom fitness function would be an example of such a "widget".

The generic recipe is:

* Generate a "pyx" file that *defines* your widgets.  A "pyx" file is a Cython source file--the "py" means Python and the "x" means extension. Again C/C++ users can view this as a source file that contains the code that actually makes these types and functions work.
* Generate a "pyxbld" file that will help to compile your custom module.

The next sections of the manual cover the following topics:

* :ref:`customFitness`

  Compiling a custom module
--------------------------------

Simply put, there's a fast/easy way, and a more involved way.  This document covers the former.  The latter involves writing a "setup.py" script for your custom module.  Documenting that is beyond the scope of this document, and the easy way will suffice in most cases.  The easy way involved using the pyximport module that comes with Cython_.

The generic recipe is:

* Have fwdpy installed!
* Have Cython_ installed!
* Your module code is in the file mymod.pyx. (This is just an example, of course--you can use any valid file name with the .pyx suffix.)
* You create a "build" file called mymod.pyx.  You do not have to write this!  Simply copy one that comes with the fwdpy source code (*e.g.* extension_tests/test_custom_fitness/test_custom_fitness.pyxbld).
* Write a Python file called "compile.py" with the following code in it:

.. code-block:: python

   import pyximport
   pyximport.install
   import mymod
   
* Figure out where the fwdpy headers are.  This is often the limiting step.  Here's a trick:

.. code-block:: bash
   #This will print the location of where the module is installed
   python -c "import fwdpy; print fwdpy.__path__"

The result on my system is:

.. code-block:: bash
		
   /home/kevin/.local/lib/python2.7/site-packages/fwdpy/__init__.pyc

Replace lib with include, delete site-packaged, and get rid of /__init__.pyc to get:

.. code-block:: bash
		
   /home/kevin/.local/include/python2.7/fwdpy

That is where the fwdpy headers are.

With that out of our way, this will compile our custom module:

.. code-block:: bash

   CPPFLAGS=-I/home/kevin/.local/include/python2.7/fwdpy python compile.py

To see a fully-worked out example, see extension_tests/run_custom_fitness.py in the fwdpy source repository.






.. _Cython: http://cython.org
.. _manual: http://docs.cython.org

