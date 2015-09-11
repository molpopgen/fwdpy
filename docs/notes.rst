Notes on the documentation
================================================

This package is split up into several modules.

The "base" functionality is in the module fwdpy.fwdpy.  This module is automatically loaded when fwdpy is imported into the name space 'fwdpy'.

In other words, code like this will allow fwdpy.fwdpy.foo to be referred to as just fwdpy.foo:

.. code-block:: python

   import fwdpy

It is important to keep this in mind when reading the reference manual.  Many objects will be documented as fwdpy.fwdpy.something but example code will refer to fwdpy.something, which is correct.

"Internal" modules
----------------------------

Some of the modules have "internal" as part of their name.  These modules represent the inner workings of the C++ layer.  As such, *they are not intended to be called directly by people using this package*.  However, they are included in the documentation in order to provide docstring tests for the developer and to document their existence for other developers interested in extending the package.
