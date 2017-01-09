
Example: Sliding windows
========================

There are two basic ways of getting sliding windows from simulated data:

1. Manually
2. Using `pylibseq <https://github.com/molpopgen/pylibseq>`__

Both work, and both are pretty easy.

.. code:: python

    #import our modules
    from __future__ import print_function
    import fwdpy as fp
    import numpy as np
    import datetime
    import time

.. code:: python

    #set up our sim
    rng = fp.GSLrng(101)
    nregions = [fp.Region(0,1,1),fp.Region(2,3,1)]
    sregions = [fp.ExpS(1,2,1,-0.1),fp.ExpS(1,2,0.1,0.001)]
    rregions = [fp.Region(0,3,1)]
    popsizes = np.array([1000]*10000,dtype=np.uint32)

.. code:: python

    #Run the sim
    pops = fp.evolve_regions(rng,4,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)

.. code:: python

    #Take samples from the simulation
    samples = [fp.get_samples(rng,i,20) for i in pops]

Calculating sliding windows
---------------------------

We are going to want non-overlapping widwos of size 0.1.

One thing to keep track of is the total size of our region, which is the
half-open interval :math:`[0,3)`

Manual method
~~~~~~~~~~~~~

Let's just do it using pure Python:

.. code:: python

    for i in samples:
        windows = []
        start = 0
        while start < 3:
            ##We will only look at neutral mutations, which are element 0 of each sampl
            window = [j[0] for j in i[0] if (j[0] >=start and j[0] < start+0.1)]
            windows.append(window)
            start += 0.1
        ##We now have a full set of windows that we can do something with
        print (len(windows))  ##There should be 30, and many will be empy


.. parsed-literal::

    30
    30
    30
    30


Using `pylibseq <https://github.com/molpopgen/pylibseq>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from libsequence.windows import Windows
    from libsequence.polytable import SimData
    for i in samples:
        ##We need to convert our list of tuples
        ##into types that pylibseq/libsequence understand:
        windows = Windows(SimData(i[0]),0.1,0.1,0,3)
        ##Now, you can analyze the windows, etc.
        print(len(windows))


.. parsed-literal::

    30
    30
    30
    30


Well, the pylibseq version is clearly more compact. Of course, you
can/should abstract the pure Python version into a standalone function.

Why would you ever use the manual version? It can save you memory. The
pylibseq version constructs an iterable list of windows, meaning that
there is an object allocated for each window. For the manual version
above, we grew a list of objects, but we could just have easily
processed them and let them go out of scope.
