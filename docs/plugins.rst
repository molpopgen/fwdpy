Overview
---------------

.. note:: This section assumes that you are familiar with running simulations using *fwdpy*.  If not, this section will make no sense at all.  This section is intended for those who want to build their own tools on top of *fwdpy*.  This document should also be read top to bottom, as later material builds on earlier stuff.

One of the more important features of *fwdpy* is that you may extend it by writing your own plugins.  By "plugin" I mean more than just some custom Python code to deal with output from the simulations.  Rather, a plugin means some combination of Cython_ or new C++11 code that adds a feature that you need for your research.

Most plugins will only need to be written using Cython_.  I will take care to mention when a plugin would require new C++11 code.  The goal is to make the latter case a very rare occurrence.  As a programming language, Cython_ may be viewed as "Python with types".  In other words, consider the following Python funtion to add two numbers:

.. code-block:: python

   def add_em_up(x,y):
      return float(x)+float(y)

Using Cython_, such a function could take one of a few forms:

.. code-block:: cython

   #As a Python function, but with the argument types specified:
   def add_em_up(double x,double y):
      return x+y

   #As a function that is only callable as a C/C++ function
   #elsewhere within a Cython module:
   cdef double add_em_up_C_CPP(double x,double y):
      return x+y

   #As a function that can be called both from C/C++ and from Python:
   cpdef double add_em_up_C_CPP_PY(double y,double y):
      return x+y

A full overview of coding in Cython_ is beyond the scope of this manual--please see their documentation for more details.  A quick overview is found `here <http://docs.cython.org/src/reference/language_basics.html>`_.

For most purposes, a plugin will consist of:

* A "cdef" function that works on some of the underlying C++ types
* A "def" function that will take objects from *fwdpy* and apply the "cdef" function to those objects.

The former is the "back end" of your plugin, and the latter is the Python interface that makes it useful.

Glossary
;;;;;;;;;;;;;;;;;;;;;;;;;;;

In one sense Cython_ is a "dialect" of Python, in that it augments Python with C and C++ data types. You may need to familiarize yourself with the basics:

* C_ data types defined
* A "return value" is the name of a type that is returned by a function.  void = noting.  Otherwise, something is returned.
* "const" is a keyword meaning "cannot be modified", or constant. C++ code attempting to modify a const variable will fail to compile.
* "class" is a complex type.  These are often composed of the fundamental types, such as C_ types.  Alternately, they may be entirely new types defined by a library, such as those descbied below.
* "object" is an instance of type or class. In the following example, vector[int] is a class type, and x is an object whose type is vector[int]:

.. code-block:: cython

   from libcpp.vector cimport vector
   cdef vector[int] x

PS, you just learned how to bring C++'s vector class into scope using Cython_.

File extensions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

The Cython_ code for a plugin will go into a ".pyx" file.  This is known as a "Pyrex_" file type, and Pyrex_ is the precuror to Cython_.

..note:: You will want to make sure that your text editor supports the .pyx and .pxd file extensions.  Cython_ provides a plugin_ for Emacs.  Casual Google searches find several options for the vi/vim family of editors.  I have no idea what, if anything, is avaiable for OS X editors like BBEdit, TextWrangler, or Xcode.  If you use them, you need to sort that out on your own.  I use Emacs and vim because they work the same way on all systems and I can edit files remotely via SSH.

For plugins that require extra C++ code, header files should have the .h or .hpp extension, and source files should have the .cc extension.  You should avoid .cpp for the following reason: Cython_ will process a pyx file into a .cpp file.  This behavior can have side-effects.  IF you have Cython_ code in foo.pyx and some additional C++ code in foo.cpp, the latter file will be over-written when Cython_ compiles foo.pyx into a C++ file (which will be called foo.cpp)!  Use .cc to avoid that.

Compiling a Cython-only plugin 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

This section describes the "easy way" to compile a plugin.  This method only applies to plugins consisting only of Cython_ code (.pyx files).

I will assume that your module will be called "foo".

You will need three files:

* foo.pyx contains your Cython_ code
* foo.pyxbuild contains extra info to tell the system that this plugin the C++11 language
* compile_foo.py is a Python script that will handle the compilation.

The contents of foo.pyx contain whatever code you need to write for your module.

foo.pyxbld contains the following:

.. code-block:: cython

   def make_ext(modname, pyxfilename):
       from distutils.extension import Extension
       return Extension(name=modname,
		sources=[pyxfilename],
                language='c++',
		extra_compile_args=['-std=c++11'])

.. note:: The "pyxbld" file will contain the same code for **all** custom modules that only depend on Cython_ code.  You just need to copy/paste that and rename it to match the prefix of your .pyx files

Finally, compile_foo.py contains *at least* the following:

.. code-block:: python

   import pyximport
   pyximport.install()
   #This import command will process foo.pyx,
   #generate a C++ source file based on it,
   #and compile it!  This only needs to happen once,
   #and recompilation will only happen if you make
   #changes to foo.pyx
   import foo

.. note:: Your Python source file can do more than just compile the module.  It could run simuations and apply your custom plugin code.  Or, you could just have one script that imports a lot of modules to compile them.
   
Finally, you need to figure out where the fwdpy headers are.  This is often the limiting step.  Here's a trick:

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

   CPPFLAGS=-I/home/kevin/.local/include/python2.7/fwdpy python compile_foo.py

To see a fully-worked out example, see extension_tests/run_custom_fitness.py in the fwdpy source repository.

Compiling a plugin that contains extra C++ code
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

**TBD**

The C++ types used in *fwdpy*
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

*fwdpy* is implemented in terms of the following:

1. fwdpp_, which is a C++11 library for forward simulation.  It is the "guts" of fwdpy and does most of the harder work that needs to be fast.
2. GSL_, which is a C library for numeric computation.  *fwdpy* and fwdpp_ use GSL_ for random number generation, fast lookup tables, etc.
3. Cython_ is the glue between fwdpp_, GSL_, and Python.
4. cythonGSL_ exposes the GSL_ to Cython_ so that we don't have to reinvent all of those wrappers.

In order to write a plugin, you need all of the above installed on your system.

Let's go through the fwdpp_ types used in *fwdpy*, and how to manipulate them in Cython_.

I will give the names of header files where these types are defined.  For fwdpp_ types, you may get all the gory details about them from that library's manual_.

The next sections go through the relevant C++ types in a "biological" order: mutation, mutation container, gamete, gamete container, diploid, etc.  Finally, we discuss the population objects that hold all of these together.

Mutation types
''''''''''''''''''''''''''''''

popgenmut
+++++++++++++++++++

This is the C++ name of the type of mutation used in *fwdpy*.  It is a mutation type with the following data members:

1. **pos**: the mutation position.  This is a double-precision floating point number.
2. **s**: the "selection coefficient" or "effect size". This is a double-precision floating point number.
3. **h**: The dominance term. This is a double-precision floating point number.
4. **neutral**: A boolean (C++ type bool) that flags the mutation as "neutral" or "selected" (true and false, respectively).
5. **g**: an unsigned (non-negative) 32-bit integer recording the generation when the mutation first appeared.
6. **label**: this is a 16 bit unsigned integer.  In practice, not much is done with it, but you can use it for adding 16 bits of extra info to a mutation type.  It was given a name in order to make use of wasted storage in the C++ type.

This type is defined in the fwdpp_ header file fwdpp/sugar/popgenmut.hpp.  It is exposed to Cython_ via fwdpy/fwdpp.pxd.

Let's show how to access these data members in Cython_.  First, we will consider the case of simply assigning each data member to another variable.  This is a pointless example, but it serves to illustrate some key concepts:

.. code-block:: cython

   #We must bring popgenmut into scope
   from fwdpy.fwdpp cimport popgenmut

   #The "void" return type mean that the function does not return a value
   cdef void do_something(const popgenmut & m) nogil:
       cdef double s = m.s
       cdef double h = m.h
       cdef double pos = m.pos
       cdef bint neutral = m.neutral

Key points:

* "nogil": this function does not act on any Python objects. As a rule of thumb, declare such functions as nogil so that they may be used in parallel programming. See Cython_ docs for more info.
* "&": this means that our function takes a "reference" to a popgenmut.  Withouth the "&", 'm' would be copied and then passed to do_something.  That copy is unnecessary and expensive, and therefore incorrect!
* "const": our function takes a const reference to a poppgenmut.  The const means that we cannot try to modify any of the data members in m.  Attempting to do so will fail to compile.

We can write non-const functions, too.  But please be aware that this gives you the ability to manipulate the population data directly.  In other words, doing the wrong thing can result in undefined behavior and crashes.

Here is a non-const function to change the selection coefficient:

.. code-block:: cython

   from fwdpy.fwdpp cimport popgenmut

   cdef void change_s(popgenmut & m, connst double news) nogil:
       #We CAN modify m, because it is not const!
       m.s = news
  
And, here is why the "&" matters:

.. code-block:: cython

   from fwdpy.fwdpp cimport popgenmut

   cdef void try_2_change_s(popgenmut m, const double news) nogil:
       #m has been passed in as a COPY, and not as a REFERENCE.
       #Thus, the COPY has its selection coefficent changed,
       #which will not have any effect on the population being
       #simulated
       m.s=news

Mutation containers and mutation counts.
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Defined in the fwdpy header "types.hpp".

A population consists of a C++ vector of mutations.  Functionally, this is very similar to the "list" type in Python.

In *fwdpy*, a vector of mutations goes by the name mcont_t (mutation container type), which refers to a **vector** of **popgenmut** objects.

A population consists of a second container containing the number of times each mutation exists in the population.  This is a C++ vector of unsigned (e.g., non-negative) integers, and *fwdpy* uses the alias ucont_t for this type.

Key points:

* These two vectors are the same length.
* The order of elements in each vector is the same, in the sense that the number of occurrences of the i-th mutation is found in the i-th position of the counts container.
* In these vectors, elements are completely unsorted with respect to age, effect size, position, or anything.

.. note:: A mutation container contains both segregating mutations *and* extinct mutations.  **It is therefore important to skip over extinct mutations when calculating things!!!  You may also want to check for, and skip, fixations**. There are a few reasons for this, efficiency being one of them; fwdpp_ will "recycle" the extinct mutations to create new mutations.  Further, the containers may or may not contain fixations.  Most "standard" population genetic simulations will remove both neutral and selected fixations from these containers.  Simulations of models such as Gaussian stabilizing selection around an optimum trait value (currently) retain fixations.  The reason for this difference is that fixations in the standard model do not contribute to differences in relative fitness.  But, in the Gaussian stabilizing selection models, they still contribute to mean trait value.  Future versions of *fwdpy* may change the behavior of Gaussian selection models, removing fixations and simply keeping track of the sum of fixed effect sizes (at least for the case of additive models).

Fixations are stored in an mcont_t and the corresponding fixation_times are stored in a ucont_t.

We will save examples of processing these objects until the section on dealing with whole-population objects

Gametes and gamete containers
'''''''''''''''''''''''''''''

Defined in the fwdpy header "types.hpp", which refers to the fwdpp_ type defined in fwdpp/forward_types.hpp.

A gamete is a simple object.  It contains the following data members:

* **n** is an unsigned integer representing how many diploids are referring to this exact copy of this gamete.
* **mutations** is a C++ vector of unsigned 64-bit integers.  Each integer is an index referring to a location in the mutation container.  This container is for neutral mutations only.  In other words, the "neutral" value of each mutation must be "true".
* **smutations** is the analog of mutations, but for "selected" mutations (e.g., those affecting fitness/trait values).  The value of each mutation referred to has "neutral" set to "false".

In C/C++, the unsigned 64-bit integer type is size_t.

.. note:: **n** is *not* equivalent to how many times a gamete exists in the population.  fwdpp_ makes no attempt to represent each identical gamete once-and-only-once.

.. note:: The integers in **mutations** and **smutations** are *sorted with respect to mutation position, in ascending order*.  Behind the scenes, fwdpp_ makes sure that this sorting order is maintained.  It allows cool things like log-time lookup of mutations based on position, etc.

*fwdpy* exposes the name gamete_t to refer to this type:

.. code-block:: cython

   from fwdpy.fwdpy cimport gamete_t

Gametes are stored in a C++ vector.  The alias for this type is gcont_t:

.. code-block:: cython

   from fwdpy.fwdpy cimport gcont_t

Again, we will save examples of processing these objects until the section on dealing with whole-population objects.

Diploids
''''''''''''''''''''''''''''''

Defined in the *fdwpy* header "types.hpp".  In fwdpp_ lingo, this is a custom_ diploid.

A diploid is a very simple C++ type with the following data members:

* **first** is a size_t (unsigned 64-bit integer) with is the location in a gamete container of the first gamete
* **second** is a size_t (unsigned 64-bit integer) with is the location in a gamete container of the second gamete
* **g** is a double-precision floating point value representing a "genetic" value
* **e** is a double-precision floating point value representing a "non-genetic" value.  For example, random noise applied to a trait
* **w** is a double-precision floating point value representing fitness.

.. note:: **g**, **e**, and **w** are *not* currently set or used by the following functions: :func:`fwdpy.fwdpy.evolve_regions`, :func:`fwdpy.fwdpy.evolve_regions_more`, :func:`fwdpy.fwdpy.evolve_regions_sampler`, and :func:`fwdpy.fwdpy.evolve_regions_sampler_fitness`.  Currently, they are used by simulations of quantitative traits.  This behavior will change in future releases, as it'll obviously be handy to have this info!

We have the following types:

.. code-block:: cython

   #This is a diploid
   from fwdpy.fwdpy cimport diplod_t
   #This is a C++ vector of diploids
   from fwdpy.fwdpy cimport dipvector_t
	 
Population types
'''''''''''''''''''''''''''''''''''''''

This is where the action is.  A population is a C++ object containing the above data types.

singlepop_t
++++++++++++++++++++++

Defined in *fwdpy* header "types.hpp".  This class inherits from the fwdpp_ tempate type singlepop (fwdpp/sugar/singlepop.hpp).

This type is used to model the following situation:

* A single deme
* A contiguous genomic region. Mutation rates, recombination rates, etc., may vary along this region via the use of :class:`fwdpy.fwdpy.Region` objects.

.. code-block:: cython

   from fwdpy.fwdpy cimport singlepop_t

It has the following data members:

* **generation**, an unsigned 32-bit integer representing the current generation. 0 is the starting value.
* **N**, an unsigned 32-bit integer representing current population size
* **mutations**, an mcont_t containing the mutations
* **mcounts**, a ucont_t containg the number of occurrences of each mutation
* **fixations**, an mcont_t containing fixations
* **fixation_times**, a cont_t containing the fixation times.
* **gametes**, a gcont_t containing the gametes
* **diploids**, a dipvector_t containing the diploids.

singlepop_t and Python
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

"Under the hood", a :class:`fwdpy.fwdpy.Spop` is a wrapper around a singlepop_t.  This type is a "Cython extension type", and is a fundamental type in *fwdpy*.  One uses containers of these types in the form of :class:`fwdpy.fwdpy.SpopVec`.

We have to get a gory detail out of the way.  A :class:`fwdpy.fwdpy.Spop` contains a C++11 "shared pointer" to a singlepop_t.  We'll see the implications of this in the recipes below.

Recipes
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

..note:: These recipes will start simply and quickly get advanced.  They'll show you what stuff from fwdpp_ is exposed to Cython and you'll also see some that *fwdpy* defines that may be of use for your own plugins.

First things first: how to go from a :class:`fwdpy.fwdpy.Spop` to a singlepopt_t in a plugin:

.. code-block:: cython

   from fwdpy.fwdpy cimport singlepop_t

   #A very boring plugin indeed!
   cdef void my_plugin_function(const singlepop_t * pop) nogil:
       pass

   #This will be the function that your plugin exposes
   #to Python:
   def foo(Spop p):
      #p is your Spop, pop is the shared pointer,
      #and pop.get() returns the raw pointer
      #to the singlepop_t
      my_plugin_function(p.pop.get())

Count the number of segregating mutations in the entire population:

.. code-block:: cython

   cdef unsigned count_muts(const singlepop_t * pop) nogil:
       cdef size_t i=0
       cdef size_t n = pop.mcounts.size()
       cdef unsigned twoN = 2*pop.popsize() #This is a member function that returns pop.N
       cdef unsigned extant=0
       for i in range(n):
	   #Check that mutation is not extinct and not fixed	
           if pop.mcounts[i] > 0 and pop.mcounts[i] < twoN:
		extant+=1
       #return our count
       return extant

Count the number of segregating *neutral* mutations in the entire population:

.. code-block:: cython

   cdef unsigned count_neutral_muts(const singlepop_t * pop) nogil:
       cdef size_t i=0
       cdef size_t n = pop.mcounts.size()
       cdef unsigned twoN = 2*pop.popsize() #This is a member function that returns pop.N
       cdef unsigned extant=0
       for i in range(n):
	   #Check that mutation is not extinct and not fixed and is neutral	
           if pop.mcounts[i] > 0 and pop.mcounts[i] < twoN and pop.mutations[i].neutral is True:
		extant+=1
       #return our count
       return extant

.. note:: Counting the number of *selected* mutations would be the same, but checking for "neutral is False".

Count the number of neutral and selected mutations per gamete, return a list of tuples to Python with that info.

.. code-block:: cython
   
   from fwdpy.fdwpy cimport singlepop_t
   #The next 2 cimports are from Cython's wrappers for the C++ standard library.
   from libcpp.vector cimport vector
   from libccp.utility cimport pair

   #KEY: a C++ pair auto-converts to a Python tuple.  A C++ vector auto converts to a list.
   #So guess what a vector of pairs converts to?

   #(A list of tuples)

   #This is a helper function.  It will count the number of segregating mutations
   #in each gamete.
   cdef int count_mutations(const vector[size_t] & keys,const ucont_t & mcounts,const unsigned twoN) nogil:
       cdef size_t i=0
       cdef size_t n=keys.size()
       cdef int rv = 0
       for i in range(n):
           #Note this next line: the i-th element in keys is an index
	   #corresponding to a location in mcounts.
           if mcounts[keys[i]] < twoN:
               rv+=1
       return rv
		
   cdef vector[pair[int,int]] mutations_per_gamete(const singlepop_t * pop) nogil:
       cdef vector[pair[int,int]] rv
       cdef size_t i = 0
       cdef size_t n = pop.gametes.size()
       cdef unsigned twoN = 2*pop.popsize()
       cdef int neutral,selected
       #Now, we go through every gamete and:
       #1. Check that it is not extinct
       #2. Go over every mutation in each gamete and make sure that it is not fixed.
       #   We do not need to check that each mutation in each gamete has a nonzero count.
       #   fwdpp ensures that an extant gamete contains extant mutations.
       for i in range(n):
           if pop.gametes[i].n > 0: #gamete is not extinct
	       #"mutations" = container of indexes to neutral mutations
               neutral = count_mutations(pop.gametes[i].mutations,pop.mcounts,twoN)
	       #"smutations" = container of indexes to selected mutations
               selected = count_mutations(pop.gametes[i].smutations,pop.mcounts,twoN)
               rv.push_back(pair[int,int](neutral,selected))
       return rv

Count the number of neutral and deleterious mutations per diploid, and return a list of tuples:

.. code-block:: cython

   cdef pair[int,int] count_mutations_diploid(const diploid_t & dip, const gcont_t & gametes, const ucont_t & mcounts, const unsigned twoN) nogil:
       #Neat: we can re-use the function defined above:
       cdef int neutral = count_mutations(gametes[dip.first].mutations,mcounts,twoN)
       cdef int selected = count_mutations(gametes[dip.first].smutations,mcounts,twoN)
       return pair[int,int](neutral,selected)

   cdef vector[pair[int,int]] mutations_per_diploid(const singlepop_t * pop) nogil:
       cdef vector[pair[int,int]] rv
       cdef size_t i=0
       cdef size_t n=pop.diploids.size()
       cdef unsigned twoN = 2*n
       cdef pair[int,int] temp
       #Now, go through every diploid:
       for i in range(n):
           temp = count_mutations_diploid(pop.diploids[i],pop.gametes,pop.mcounts,twoN)
           rv.push_back(temp)
       return rv

Time to up the complexity level with the next examples.

Population mean fitness under a multiplicative model.  We will calculate the mean fitness of the population by *explicitly calculating the fitness of each diploid*.  We will make this calculation under a multiplicative model, :math:`w = \prod_i(1+I(i))`, where :math:`I(i)` is :math:`sh` or :math:`scaling*s` for hetero- and homo- zygous mutation positions, respectively.

Some comments:

1. We will use fwdpp's multiplicative_diploid class to do this calculation.
2. We will use a numpy_ array to store the fitness of every diploid and retuirn the mean of the array as the calculation of mean fitness.

Thus, this example shows us how to:

1. Use more fwdpp
2. Integrate numpy_ with Cython_ code via "typed array views"

.. code-block:: cython

   import numpy as np;
   from cython.view cimport array as cvarray
   from fwdpy.fwdpp cimport multiplicative_diploid

   cdef void wbar_multiplicative_details(const singlepop_t * pop, double[:] w, const double scaling) nogil:
       cdef multiplicative_diploid wfxn
       cdef size_t i=0, n=pop.diploids.size()
       for i in range(n):
           #Here is the trick.  wfxn is a C++ class, but it is also a function!
           #Further, it is a template function.  Cython is not willing to just let
           #the C++ compiler figure out the types here, so we have to explicitly use typecasts,
           #which is what the <foo>bar is: type cast a bar to a foo.  This has NO RUNTIME PENALTY!!!
           #Yes, we also have to cast the scaling parameter, even though it is not a template parameter.
           w[i] = wfxn(<diploid_t>pop.diploids[i],<gcont_t>pop.gametes,<mcont_t>pop.mutations,<double>scaling)

   def wbar_mutiplicative(Spop p, const double scaling):
       """
       This is our Python function.
       """
       #Create the numpy array
       w=np.array(p.popsize(),dtype=np.float64)
       #Call our Cython function:
       wbar_multiplicative_details(p.pop.get(),w[:],scaling)
       #return mean fitness:
       return w.mean()

.. note:: The above function is only useful if you run it on a population using the same "scaling" that you used to simulate!!!

       
.. _Cython: http://www.cython.org
.. _fwdpp: http://molpopgen.github.io/fwdpp
.. _GSL:  http://gnu.org/software/gsl
.. _cythonGSL: https://pypi.python.org/pypi/CythonGSL
.. _manual: http://molpopgen.github.io/fwdpp/doc/html/index.html
.. _custom: http://molpopgen.github.io/fwdpp/doc/html/d2/dcd/md_md_customdip.html
.. _C: https://en.wikipedia.org/wiki/C_data_types
.. _numpy: http://www.numpy.org
.. _Pyrex: https://www.cosc.canterbury.ac.nz/greg.ewing/python/Pyrex/
.. _plugin: https://github.com/cython/cython/blob/master/Tools/cython-mode.el
