Notes on the documentation
================================================

This package is split up into several modules.

The "base" functionality is in the module fwdpy.fwdpy.  This module is automatically loaded when fwdpy is imported into the name space 'fwdpy'.

In other words, code like this will allow fwdpy.fwdpy.foo to be referred to as just fwdpy.foo:

.. code-block:: python

   import fwdpy

It is important to keep this in mind when reading the reference manual.  Many objects will be documented as fwdpy.fwdpy.something but example code will refer to fwdpy.something, which is correct.
