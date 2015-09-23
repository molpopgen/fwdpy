Modeling changes in population size
============================================

Many simulation functions in this package accept an array of integers as an argument.  See, for example, the 'nlist' argument of :func:`fwdpy.fwdpy.evolve_regions`.  This array represents the population size over time, and the length of the array is the number of generations to simulate.  The array is a NumPy_ array of 32-bit unsigned integers, and a "view" of the array is passed to the simulation function.

Simple example
----------------

Let's look at an example:

.. code-block:: python

		import numpy as np
		#population size
		N=1000
		#nlist corresponds to a constant population size for 10N generations
		#note the "dtype" argument.  Without it, we'd be defaulting to int64,
		#which is a 64-bit signed integer.
		nlist=np.array([N]*(10*N),dtype=np.uint32)
		#This is a 'view' of the array starting from the beginning:
		nlist[0:]

Simple bottleneck
--------------------------------

In order to change population size, one simply has to change the values in the "nlist".   For example, here is a population bottleneck:

.. code-block:: python

		import numpy as np
		N=1000
		#Evolve for 10N generations,
		#bottleneck to 0.25N for 100 generations,
		#recover to N for 50 generations
		nlist = np.concatenate(([N]*(10*N),[int(0.25*N)]*100,[N]*50)).astype(np.int32)

Please note the last command, which changes the concatenated array from an array of 64 bit signed integers to 32 bit unsigned integers.

Exponential growth
------------------------------------

Now, let's do population growth, where we evolve for 10N generations, and then grow the population five fold in the next 500 generations.

.. code-block:: python

		import numpy as np
		import math

		N=1000
		N2=5*N
		tgrowth=500
		
		#G is the growth rate
		G = math.exp( (math.log(N2)-math.log(N))/float(tgrowth) )
		
		nlist = np.array([N]*(10*N+tgrowth),dtype=np.uint32)
		
		#Now, modify the list according to expoential growth rate
		for i in range(tgrowth):
		     nlist[10*N+i] = round( N*math.pow(G,i+1) )

Potential caveats:
---------------------

* Forgetting the 'dtype'.  You will get a run-time error if you don't get the integer type right.  Python will raise a ValueError exception about a buffer type mismatch.

Rationale
---------------------

Why do things with NumPy_ arrays?  Lots of reasons:

1. They are fast
2. The uint32 is the same type used in fwdpp_
3. The 32 bit integer takes half the memory as the default 64 bit intger type of a Python list.
4. Cython lets us directly pass the underlying data to C++, eliminating the need for a copy when going from Python to C++.

.. _NumPy: http://www.numpy.org
.. _fwdpp: http://molpopgen.github.io/fwdpp
