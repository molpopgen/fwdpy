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

1. As a function of the genotypes of each variable position in the diploid.  In other words, fitness depends on whether or not a diploid is heterozyguous (Aa) or homozygous (aa) for a specific mutation
2. As a function of the two haplotypes in each diploid.

The standard additive and multiplicative models with dominance are an example of the first type.   Let's call these models "site-based" fitness models.  We'll call the latter "haplotype-based" fitness models.

To implement a site-based model, we need three different functions:

1. A function to adjust fitness/trait values based on Aa genotypes.
2. A function to adjust fitness/trait values based on aa genotypes.
3. A function to return the final fitness/trait value.  Below, we'll see why this is necessary.

For haplotype-based models, we need the following:

1. A function to process a haplotype.  This function must return a single number.
2. A function to combined the values for each haplotype into a single, final value.

The relevant C++ types
'''''''''''''''''''''''''''''''''''''

The following C++ types are relevant for implementing custom fitness models:

1. popgenmut is the name of the type of mutation.  It has a position (pos), selection coefficient/effect size (s), and a dominance term (h).  All of these are floating point types (double in C/C++).
2. mcont_t is the type representing a container of popgenmuts.
3. gamete is a gamete.  A gamete contains two containers called mutations and smutations.  These are containers of unsigned integers (*e.g.* strictly non-negative integers).  The integers represent the locations in an mcont_t where mutations are found.  These integers are sorted according to their mutation position (in ascending order).  The mutations container is for "neutral" mutations, and smutations is for "selected" mutations.
4. gcont_t is a container of gametes.

Under the hood, these types are all aliases ("typedefs") to objects that fwdpp_ knows how to work with.  From a Cython/Python, popgenmut and gamete are classes, and mcont_t and gcont_t play the role of lists.

The relevant Cython extension types
'''''''''''''''''''''''''''''''''''''''

These are defined in the source file fwdpy/fitness.pxd.  The types are:

* :class:`fwdpy.fitness.SpopFitness`
* :class:`fwdpy.fitness.MlocusFitness`

Custom types will be extensions of these base types.

Helper functions
'''''''''''''''''''''''''''''''''''''''

*fwdpy* provides several "helper" functions to simplify writing custom fitness functions.  These are defined in fwdpy/fitness.pxd, and that file is heavily commented.

In general, you'll probably just want to import everything from the fitness module:

.. code-block:: cython

   from fwdpy.fitness cimport *

To see how these functions are used, take a look at fwdpy/fitness.pyx, which is where the built-in fitness models are implemented.  The built-in models are implemented using the machinery for implementing custom fitness models.
		 

.. _fwdpp: http://molpopgen.github.io/fwdpp/
