
Copying populations in memory
=============================

It is possible to make in-memory copies of populations, and then evolve
the two populations independently. Doing so makes use of the
serialization features in the *fwdpy.fwdpyio* module. You do not need to
import that module. The routines used below do so transparently wihtout
polluting the namespace for your session.

.. code:: python

    import fwdpy as fp
    import numpy as np
    #Let's set up a simulation and evolve some populations...
    nregions = [fp.Region(0,1,1),fp.Region(2,3,1)]
    sregions = [fp.ExpS(1,2,1,-0.1),fp.ExpS(1,2,0.01,0.001)]
    rregions = [fp.Region(0,3,1)]
    rng = fp.GSLrng(201)
    popsizes = np.array([1000],dtype=np.uint32)
    popsizes=np.tile(popsizes,10000)
    pops = fp.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.001,0.001,nregions,sregions,rregions)

All of the above is pretty standard. Let's make a copy of *pops*:

.. code:: python

    pops2 = fp.copypops(pops)

.. code:: python

    pops




.. parsed-literal::

    <fwdpy.fwdpy.popvec at 0x7faaace9e2b8>



.. code:: python

    pops2




.. parsed-literal::

    <fwdpy.fwdpy.popvec at 0x7faaace9e310>



We see above that the two pops are occuppying different spaces in
memory.

Let's now evolve them further.

.. code:: python

    #Evolve "pops"
    fp.evolve_regions_more(rng,pops,popsizes[0:],0.001,0.001,0.001,nregions,sregions,rregions)

.. code:: python

    #Now, let's evolve "pops2", and end with a short bottleneck (100 gens at N=500)
    popsizes = np.array([1000]*10000 + [500]*100,dtype=np.uint32)
    fp.evolve_regions_more(rng,pops2,popsizes[0:],0.001,0.001,0.001,nregions,sregions,rregions)

.. code:: python

    for i in range(len(pops)):
        print pops[i].popsize(),' ',pops2[i].popsize()
        print pops[i].gen(),' ',pops2[i].gen()


.. parsed-literal::

    1000   500
    20000   20100


It is really that easy. Now, you can evolve a population to equilibrium
and then evolve copies of it within a for loop to get multiple
quasi-independent replicates of bottlenecks, etc.
