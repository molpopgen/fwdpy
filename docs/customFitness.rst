.. _customFitness:

Implementing custom fitness functions
==============================================

This section covers how to write your own "fitness" functions.  (I'm using the term "fitness" here somewhat loosely, as sometimes such a function may be used to calculate a genetic value for a diploid, as in the case of simulations of quantitative traits.)

Fundamentally, there are two ways that these types of functions work:

1. As a function of the genotypes of each variable position in the diploid.  In other words, fitness depends on whether or not each site is heterozyguous (Aa) or homozygous (aa).
2. As a function of the two haplotypes in each diploid.

The standard additive and multiplicative models with dominance are an example of the first type.   Let's call these models "site-based" fitness models.

The standard additive model with dominance
--------------------------------------------------

Here, we re-implement the additve model as if it were a custom model.  Working through this standard model will help us understand how to implement custom "site-based" fitness models in general.

This is a very common model using in population genetics.  The fitnesses of the three diploids are 1, 1+sh, and 1+2s for genotypes AA, Aa, and aa respectively.  Across variable sites, the final fitness of a diploid is :math:`w = max(0,1+\sum_i w_i)`, where :math:`w_i` is 0, sh, or 2s depending on the genotype at the :math:`i^{th}` site.

The formula for additive fitness tells us that we have to do the following:

1. Initialize fitness, w, to 0.
2. Decide how the Aa and aa genotypes will be used to update w
3. Return max(0,1+w).  It is important that we ensure that the minimum possible fitness is zero, otherwise fwdpp_ will not be happy and our simulations will exit in a **very** ungraceful fashion at run time.

Thus, for "site-based" models, we need:

1. A starting value for w
2. Functions for updating w based on Aa and aa genotypes at each site
3. A "finalizer" function that returns the final fitness.

In general, these functions will be very short and simple.

First, let us define our custom data types in a file called "myCustomFitness.pxd":

.. code-block:: cython

   from fwdpy.fwdpp cimport popgenmut
   from fwdpy.fitness cimport singlepopFitness

   #This is the name of our Cython extension type.
   #It is an object type inheriting from singlepopFitness.
   #We don't need to define any functions here, so we "pass"
   #on its implementation
   cdef myAdditiveFitness(singlepopFitness):
       pass

   #This is the function that updates w when a genotype
   #is heterozygous. Note the "signature" of this function:
   #1. It returns nothing, which is what the void means
   #2. It takes a reference to current fitness (w) as its first argument
   #3. It takes a constant reference to a mutation (m) as its second argument
   #The function add the product s*h to the current value of w.
   cdef inline void wAa( double & w, const popgenmut & m):
       #This manipulation of w should look odd, and it is.
       #Cython has an issue in not allowing non-const references to
       #be updated, which is incorrect for C++.  This syntax is
       #a workaround.  Just memorize it...
       (&w)[0] += m.s*m.h

   #This function updates w for a homozygous genotye,
   #for which case w += 2s
   cdef inline void waa( double & w, const popgenmut & m):
       (&w)[0] += 2.0*m.s

Now, we have to make our type "myFitness" actually work in Python, which means defining its
behavior in "myCustomFitness.pyx":

.. code-block:: cython
		
   from fwdpy.fitness cimport genotype_fitness_updater,fitness_function_finalizer,make_custom_fitness,return_w_plus1

   #This defines the behavior of our type.
   #The __cinit__ function is a "constructor"
   #which creates a fitness function that can be
   #passed to the C++ code behind fwdpy. Where you see
   #the <x>y syntax, that is a signal to "cast"
   #your functions to the correct type of "function pointer"
   #(a C/C++ thing) that can be used internally.
   #The return_w_plus1 is built in to fwdpy.fitness,
   #and will return max(0,w).
   #Finally, the 0.0 is the starting value of w.
   cdef class myAdditiveFitness(singlepopFitness):
		def __cinit__(self):
		self.wfxn = make_custom_fitness(<genotype_fitness_updater>wAa
		<genotype_fitness_updater>waa,
		<fitness_function_finalizer>return_w_plus1,
                0.0)


Let's review.  At the beginning, we realized we needed to do the following:

1. Initialize fitness, w, to 0.
2. Decide how the Aa and aa genotypes will be used to update w
3. Return max(0,1+w).  It is important that we ensure that the minimum possible fitness is zero, otherwise fwdpp_ will not be happy and our simulations will exit in a **very** ungraceful fashion at run time.

Here's how these things were accomplised:

1. Pass 0.0 as the final argument to make_custom_fitness
2. Define the wAa and waa functions
3. Use of the Cython function return_w_plus1, which is part of fwdpy.fitness

The call to make_custom_fitness passes our custom functions along, and returns an object representing a call to stuff in fwdpp_ that will apply the wAa and waa functions to the appropriate sites in a diploid.  You don't need to understand any of the nasty C++ of what this "wfxn" thing really is--that's all taken care of.
   
.. _fwdpp: http://molpopgen.github.io/fwdpp/
