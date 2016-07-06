
Temporal sampling
=================

_fwdpy_ allows various things to be recorded over time during a simulation.  A family of objects referred to as "temporal samplers" perform these tasks.  All such objects are derived from the base class :class:`fwdpy.fwdpy.TemporalSampler`.

Sampling nothing
----------------

Doing nothing is useful for evolving a population to equilibrium.  The relevant class is :class:`fwdpy.fwdpy.NothingSampler`.

For convenience, :func:`fwdpy.fwdpy.evolve_regions` and :func:`fwdpy.fwdpy.evolve_regions_more` and :func:`fwdpy.fwdpy.evolve_regions_fitness` all implicitly use :class:`fwdpy.fwdpy.NothingSampler`.

Let's evolve 40 populations to mutation-drift equilibrium:

.. code:: python

    import fwdpy as fp
    import numpy as np
    import pandas as pd
    nregions=[fp.Region(0,1,1)]
    sregions=[fp.GammaS(0,1,0.1,0.1,0.1,1.0),
              fp.GammaS(0,1,0.9,-0.2,9.0,0.0)
             ]
    recregions=nregions
    N=1000
    nlist=np.array([N]*(10*N),dtype=np.uint32)
    mutrate_neutral=50.0/float(4*N)
    recrate=mutrate_neutral
    mutrate_sel=mutrate_neutral*0.2
    rng=fp.GSLrng(101)
    pops=fp.SpopVec(40,1000)
    sampler=fp.NothingSampler(len(pops))
    #This function implicitly uses a "nothing sampler"
    fp.evolve_regions_sampler(rng,pops,sampler,nlist,
                                mutrate_neutral,
                                0.0,   #No selected mutations....
                                  recrate,
                                  nregions,sregions,recregions,
                                   #Only sample every 10N generations,
                                   #which is fine b/c we're not sampling anything
                                  10*N)

Take samples from population
----------------------------

Example using :class:`fwdpy.fwdpy.PopSampler`

.. code:: python

    #Take sample of size n=20
    sampler=fp.PopSampler(len(pops),20,rng)
    fp.evolve_regions_sampler(rng,pops,sampler,
                              nlist[:N], #Evolve for N generations
                                mutrate_neutral,
                                mutrate_sel,   
                                  recrate,
                                  nregions,sregions,recregions,
                                #Sampler every 100 generations
                                  100)

.. code:: python

    #Get the data from the sampler
    x=sampler.get()
    print(len(x))


.. parsed-literal::

    40


.. note:: len(pops) == len(x) !!!

The output from this sampler type is a bit complex. Each element in x is
itself a list:

.. code:: python

    print type(x[0])


.. parsed-literal::

    <type 'list'>


Each element in x[0] is a tuple:

.. code:: python

    print type(x[0][0])


.. parsed-literal::

    <type 'tuple'>


The first element of each tuple is the generation when the sample was
taken:

.. code:: python

    print x[0][0][0]


.. parsed-literal::

    10100


The rest is the sample info, in the same format as the output from :func:`fwdpy.fwdpy.get_samples`:

.. code:: python

    print x[0][0][1]


.. parsed-literal::

    {'sh': [(-0.12004595381424654, 0.0), (-0.170667130152595, 0.0)], 'genotypes': ([(0.11155974119901657, '10000000000000000000'), (0.24918362661264837, '00000000000000100100'), (0.42917620693333447, '00100000000000000000'), (0.519921897444874, '00000000000000010000'), (0.5717735027428716, '01000000000000000000'), (0.6115675517357886, '00000100000000000000'), (0.6820425151381642, '00100000000000000000'), (0.7076556014362723, '00000000000100000000'), (0.7335569099523127, '00100000000000000000'), (0.7683701689820737, '00100000000000000000'), (0.778331029927358, '00000000010000000000'), (0.7873499824199826, '00000000000001000000'), (0.9194794276263565, '00010000000000000000'), (0.9491714215837419, '00000000000000010000')], [(0.8127563807647675, '00000000100000000000'), (0.9910217146389186, '00010000000000000000')])}


.. code:: python

    #Selection coefficients and dominance for each selected mutation:
    print x[0][0][1]['sh']


.. parsed-literal::

    [(-0.12004595381424654, 0.0), (-0.170667130152595, 0.0)]


.. code:: python

    #Neutral mutations:
    print x[0][0][1]['genotypes'][0]


.. parsed-literal::

    [(0.11155974119901657, '10000000000000000000'), (0.24918362661264837, '00000000000000100100'), (0.42917620693333447, '00100000000000000000'), (0.519921897444874, '00000000000000010000'), (0.5717735027428716, '01000000000000000000'), (0.6115675517357886, '00000100000000000000'), (0.6820425151381642, '00100000000000000000'), (0.7076556014362723, '00000000000100000000'), (0.7335569099523127, '00100000000000000000'), (0.7683701689820737, '00100000000000000000'), (0.778331029927358, '00000000010000000000'), (0.7873499824199826, '00000000000001000000'), (0.9194794276263565, '00010000000000000000'), (0.9491714215837419, '00000000000000010000')]


.. code:: python

    #Selected mutations:
    print x[0][0][1]['genotypes'][1]


.. parsed-literal::

    [(0.8127563807647675, '00000000100000000000'), (0.9910217146389186, '00010000000000000000')]


These "genotypes" blocks can be used to caculate summary statistics. See
the example on using `pylibseq <http://molpopgen.github.io/pylibseq/>`__
for that task.

Tracking mutation frequencies
-----------------------------

See the example on tracking mutation frequencies.

The relevant class is :class:`fwdpy.fwdpy.FreqSampler`.
