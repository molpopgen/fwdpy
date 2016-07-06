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

Your custom types will by Cython extension types (think classes) that extend these base types.

The standard additive model with dominance
--------------------------------------------------

.. note:: This model is already implemented in *fwdpy* via :class:`fwdpy.fitness.singlepopAdditive`.  This section is to illustrate what it would take to implement it as an extension to *fwdpy*.

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

First, let us define our custom data types in a file called "myCustomFitness.pyx":

.. code-block:: cython

   from fwdpy.fwdpp cimport popgenmut
   from fwdpy.fitness cimport singlepopFitness

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
   cdef inline void waa(double & w, const popgenmut & m):
       (&w)[0] += 2.0*m.s

Now, we have to make our type "myFitness" actually work in Python, which means defining its
behavior in "myCustomFitness.pyx":

.. code-block:: cython
		
   from fwdpy.fitness cimport genotype_fitness_updater,fitness_function_finalizer,make_custom_fitness,return_w_plus1

   #This is our type that we will use in Python.
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

Without all the comments, the custom fitness setup is quite short:

.. code-block:: cython

   from fwdpy.fwdpp cimport popgenmut
   from fwdpy.fitness cimport singlepopFitness
   from fwdpy.fitness cimport genotype_fitness_updater,fitness_function_finalizer,make_custom_fitness,return_w_plus1
   
   cdef inline void wAa( double & w, const popgenmut & m):
       (&w)[0] += m.s*m.h

   cdef inline void waa(double & w, const popgenmut & m):
       (&w)[0] += 2.0*m.s

   cdef class myAdditiveFitness(singlepopFitness):
		def __cinit__(self):
		self.wfxn = make_custom_fitness(<genotype_fitness_updater>wAa
		<genotype_fitness_updater>waa,
		<fitness_function_finalizer>return_w_plus1,
                0.0)
   

The standard multiplicative model with dominance
--------------------------------------------------

.. note:: This model is already implemented in *fwdpy* via :class:`fwdpy.fitness.singlepopMulti`. This section is to illustrate what it would take to implement it as an extension to *fwdpy*.
	  
The multiplicative model would also be easy to implement: :math:`w = \prod_i (1+w_i)`.  The only difference is that we have a starting fitness of 1 instead of 0, and we can return w instead of 1+w at the end. The full setup would look like this:

.. code-block:: cython

   from fwdpy.fwdpp cimport popgenmut
   from fwdpy.fitness cimport singlepopFitness
   from fwdpy.fitness cimport genotype_fitness_updater,fitness_function_finalizer,make_custom_fitness,return_w_plus1
   
   cdef inline void wAa( double & w, const popgenmut & m):
	#Note difference from additive model:
       (&w)[0] *= (1.0+m.s*m.h)

   cdef inline void waa(double & w, const popgenmut & m):
   	#Note difference from additive model:
       (&w)[0] *= (1.0+2.0*m.s)

   cdef class myMultiplicativeFitness(singlepopFitness):
		def __cinit__(self):
		self.wfxn = make_custom_fitness(<genotype_fitness_updater>wAa
		<genotype_fitness_updater>waa,
		#returns max(0,w) instead of max(0,1+w)
		<fitness_function_finalizer>return_w,
		#starting value of 1.0 instead of 0.0
                1.0)

A custom haplotype-based model
------------------------------------------------------

.. note:: This model is already implemented in *fwdpy* via :class:`fwdpy.fitness.singlepopGBR`. This section is to illustrate what it would take to implement it as an extension to *fwdpy*.


Let's implement the recessive haplotype model from Thornton *et al.* (2013) PLoS Genetics.  For this model, the effect size of a haplotype is additive over all mutations.  The final fitness value is the geometric mean of the two haplotype effect sizes, giving :math:`w = \sqrt{e1*e2}`.

Haplotype-based models differ from site-based models:

1. You need a function to operate on a *gamete* rather than on a *mutation*
2. You need a function to combine the values calculated from each gamete.

Here is the model:

.. code-block:: cython

   from fwdpy.fwdpp cimport popgenmut
   from fwdpy.fwdpy cimport gamete_t,mcont_t
   from fwdpy.fitness cimport singlepopFitness
   from fwdpy.fitness cimport make_custom_haplotype_fitness,haplotype_fitness_fxn,haplotype_fitness_fxn_finalizer
   #This is important--we need to import something from the C++ library:
   from libcpp.vector cimport vector

   #This sums effect sizes on a haplotype
   cdef inline double addEsizes(const gamete_t & g, const mcont_t & m):
       cdef size_t i=0,n=g.smutations.size()
       cdef double sum = 0.0
       while i<n:
          sum+=m[g.smutations[i]].s
          i+=1
       return sum

   cdef inline double geomean(double e1, double e2):
       return sqrt(e1*e2)

   cdef class GBRfitness(singlepopFitness):
       def __cinit__(self):
          self.wfxn = make_custom_haplotype_fitness(<haplotype_fitness_fxn>addEsizes,
          <haplotype_fitness_fxn_finalizer>geomean)


Let's work through that "addEsizes" function in more detail.  Let's start with its "function signature" (or prototype):

.. code-block:: cython

   double addEsizes(const gamete_t & g, const mcont_t & m)
   

This function returns a double representing the effect of the gamete (g) on fitness.  Inside of fwdpp_, the C++ library that *fwdpy* uses, a gamete contains keys to the mutations it contains.  Thus, we need the container of mutations in the population (m) in order to do anything useful.  The "const" bit means that the function may not modify the data containg in g nor that in m.

In fwdpp_, a gamete contains is "neutral" and "selected" mutations in separate containers, which speeds up fitness calculations. Thus, a fitness function will iterate over the latter container, whose name is *smutations*.

Now we can look through the function body:

.. code-block:: cython

   #i is a "dummy" counter, n is now many selected mutations g contains
   cdef size_t i=0,n=g.smutations.size()
   #initialize sum of effect sizes to 0
   cdef double sum = 0.0
   #loop over each mutation in g.smutations
   while i<n:
       #increment the sum of effect sizes
       sum+=m[g.smutations[i]].s
       #update our counter.  Important, else we get an infinite loop!
       i+=1
   #return sum over fitness effects
   return sum
		 

.. _fwdpp: http://molpopgen.github.io/fwdpp/
