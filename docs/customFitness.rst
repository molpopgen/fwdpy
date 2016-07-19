.. _customFitness:

Implementing custom fitness functions
==============================================

This section covers how to write your own "fitness" functions.  First, we cover important concepts.  Then, we describe the set of functions that *fwdpy* provides to help implementing custom fitness models.  Then, we work through some examples.

Please note that I'm using the term "fitness" here somewhat loosely, as sometimes such a function may be used to calculate a genetic value for a diploid, as in the case of simulations of quantitative traits.

Concepts
--------------------------

Types of fitness models
''''''''''''''''''''''''''''

Fundamentally, there are two ways that these types of functions work:

1. Fitness depends on whether or not a mutation is found in heterozyzous (Aa) or homozygous (aa) state in a diploid.  We'll call this a "site-based" fitness model.
2. Fitness depends on properties of separate haplotypes, which are then combined into a single value for fitness.  We'll call these "haplotype-based" fitness models.

Types of fitness model *objects*
''''''''''''''''''''''''''''''''''''''''''''''

Further, if all we need to know about a diploid in order to calculate its fitness is what gametes (and therefore what mutations) a diploid contains, then our fitness model requires *no extra data*.  Such a situation can be represented via a "stateless" object, one requiring no extra information about the population or its history.

If however, fitness depends on comparing a diploid to all other diploids, then our code for calculating fitness must keep track information about the entire population.  Such situations require "stateful" fitness objects.

*fwdpy* supports both stateful and stateless fitness objects.  Further, stateful objects can be implemented in terms of any valid C++ type.

The relevant Cython extension types
'''''''''''''''''''''''''''''''''''''''

These are defined in the source file fwdpy/fitness.pxd.  The types are:

* :class:`fwdpy.fitness.SpopFitness`
* :class:`fwdpy.fitness.MlocusFitness`

Custom types will be extensions of these base types. For example:

.. code-block:: cython

   cdef class MyFitness(SpopFitness):
       pass

Helper functions
'''''''''''''''''''''''''''''''''''''''

*fwdpy* provides several "helper" functions to simplify writing custom fitness functions.  These are defined in fwdpy/fitness.pxd, and that file is heavily commented.

In general, you'll probably just want to import everything from the fitness module:

.. code-block:: cython

   from fwdpy.fitness cimport *

To see how these functions are used, take a look at fwdpy/fitness.pyx, which is where the built-in fitness models are implemented.  The built-in models are implemented using the machinery for implementing custom fitness models.
		 
Examples
'''''''''''''''''''''''''

More examples of both stateless and stateful custom fitness functions are found in the repo `fwdpy_plugins <http://github.com/molpopgen/fwdpy_plugins>`_.  The examples are kept in a different source repo for a few reasons.  First, it is easiest to have a repo that we can test instead of static documentation that may drift from reality.  Second, compiling plugins required that *fwdpy* be installed and that the plugin code is not in the *fwdpy* repo.

.. _fwdpp: http://molpopgen.github.io/fwdpp/
