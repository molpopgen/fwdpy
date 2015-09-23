
Example: background selection
=============================

Setting up the simulation
-------------------------

-  Neutral mutations will occur on the interval :math:`[0,1)`.
-  Strongly-deleterious mutations will occur on the intervals
   :math:`[-1,0)` and :math:`[1,2)`.
-  Recombination will be uniform throughout the region.

.. code:: python

    #The next two lines let plots show up
    #in the notbook:
    %matplotlib inline 
    %pylab inline 
    #Use Pyhon 3's print a a function.
    #This future-proofs the code in the notebook
    from __future__ import print_function
    #Import fwdpy.  Give it a shorter name
    import fwdpy as fp
    #Import the module for summary statistics. Give it a shorter name
    import fwdpy.libseq as libseq
    ##Other libs we need
    import numpy as np
    import pandas
    import math


.. parsed-literal::

    Populating the interactive namespace from numpy and matplotlib


.. code:: python

    #We're going to do some plots at the end.
    import matplotlib
    import matplotlib.pyplot as plt

Establishing 'regions' for mutation and recombination
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    # Where neutral mutations occur:
    nregions = [fp.Region(beg=0,end=1,weight=1)]

.. code:: python

    # Where selected mutations occur:
    sregions = [fp.ConstantS(beg=-1,end=0,weight=1,s=-0.05,h=1),
                fp.ConstantS(beg=1,end=2,weight=1,s=-0.05,h=1)]

.. code:: python

    # Recombination:
    recregions = [fp.Region(beg=-1,end=2,weight=1)]

Population size and simulation length
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    #Population size
    N=1000
    #We'll evolve for 10N generations.
    #nlist is a list of population sizes over time.
    #len(nlist) is the length of the simulation
    #We use numpy arrays for speed and optimised RAM
    #use.  Note the dtype=np.uint32, which means 32-bit
    #unsigned integer. Failure to use this type will
    #cause a run-time error.
    nlist = np.array([N]*10*N,dtype=np.uint32)

.. code:: python

    #Initalize a random number generator with seed value of 101
    rng = fp.GSLrng(101)

.. code:: python

    #Simulate 4 replicate populations.  This uses C++11 threads behind the scenes:
    pops = fp.evolve_regions(rng,       #The random number generator 
                             4,         #The number of pops to simulate = number of threads to use.
                             N,         #Initial population size for each of the 4 demes
                             nlist[0:], #List of population sizes over time.
                             0.005,     #Neutral mutation rate (per gamete, per generation)
                             0.01,      #Deleterious mutation rate (per gamete, per generation)
                             0.005,     #Recombination rate (per diploid, per generation)
                             nregions,  #Defined above
                             sregions,  #Defined above
                             recregions)#Defined above

.. code:: python

    #Now, pops is a Python list with len(pops) = 4
    #Each element's type is fwdpy.singlepop
    print(len(pops))
    for i in range(len(pops)):
        print(type(pops[i]))
                    


.. parsed-literal::

    4
    <type 'fwdpy.fwdpy.singlepop'>
    <type 'fwdpy.fwdpy.singlepop'>
    <type 'fwdpy.fwdpy.singlepop'>
    <type 'fwdpy.fwdpy.singlepop'>


Taking samples from simulated populations
-----------------------------------------

.. code:: python

    #Use a list comprehension to get a random sample of size
    #n = 20 from each replicate
    samples = [fp.get_samples(rng,i,20) for i in pops]
    
    #Samples is now a list of tuples of two lists.
    #Each list contains tuples of mutation positions and genotypes.
    #The first list represents neutral variants.
    #The second list represents variants affecting fitness ('selected' variants)
    #We will manipulate/analyze these genotypes, etc.,
    #in a later example
    for i in samples:
        print ("A sample from a population is a ",type(i))
        
    print(len(samples))


.. parsed-literal::

    A sample from a population is a  <type 'tuple'>
    A sample from a population is a  <type 'tuple'>
    A sample from a population is a  <type 'tuple'>
    A sample from a population is a  <type 'tuple'>
    4


Getting additional information about samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    #Again, use list comprehension to get the 'details' of each sample
    #Given that each object in samples is a tuple, and that the second
    #item in each tuple represents selected mutations, i[1] in the line
    #below means that we are getting the mutation information only for
    #selected variants
    details = [fp.get_sample_details(i[1],j) for i,j in zip(samples,pops)]

.. code:: python

    #details is now a list of pandas DataFrame objects
    #Each DataFrame has the following columns:
    #  a: mutation age (in generations)
    #  h: dominance of the mutation
    #  p: frequency of the mutation in the population
    #  s: selection coefficient of the mutation
    for i in details:
        print(i)


.. parsed-literal::

        a  h       p     s
    0  35  1  0.0075 -0.05
    1   8  1  0.0040 -0.05
    2  18  1  0.0030 -0.05
        a  h       p     s
    0  19  1  0.0100 -0.05
    1  13  1  0.0025 -0.05
        a  h       p     s
    0  10  1  0.0090 -0.05
    1  29  1  0.0060 -0.05
    2  15  1  0.0035 -0.05
        a  h       p     s
    0   1  1  0.0005 -0.05
    1  20  1  0.0065 -0.05
    2  33  1  0.0080 -0.05
    3   1  1  0.0005 -0.05


.. code:: python

    #The order of the rows in each DataFrame is the
    #same as the order as the objects in 'samples':
    for i in range(len(samples)):
        print("Number of sites in samples[",i,"] = ",
              len(samples[i][1]),". Number of rows in DataFrame ",i,
              " = ",len(details[i].index),sep="")


.. parsed-literal::

    Number of sites in samples[0] = 3. Number of rows in DataFrame 0 = 3
    Number of sites in samples[1] = 2. Number of rows in DataFrame 1 = 2
    Number of sites in samples[2] = 3. Number of rows in DataFrame 2 = 3
    Number of sites in samples[3] = 4. Number of rows in DataFrame 3 = 4


.. code:: python

    #Pandas DataFrames are cool.
    #Let's add a column to each DataFrame
    #specifying the mutation position,
    #count of derived state,
    #and a "replicate ID"
    for i in range(len(details)):
        ##samples[i][1] again is the selected mutations in the sample taken
        ##from the i-th replicate
        details[i]['pos']=[x[0] for x in samples[i][1]]               #Mutation position
        details[i]['count']=[ x[1].count('1') for x in samples[i][1]] #No. occurrences of derived state in sample
        details[i]['id']=[i]*len(details[i].index)                    #Replicate id

.. code:: python

    ##Merge into 1 big DataFrame:
    BigTable = pandas.concat(details)
    
    print("This is the merged table:")
    print(BigTable)


.. parsed-literal::

    This is the merged table:
        a  h       p     s       pos  count  id
    0  35  1  0.0075 -0.05  1.191983      1   0
    1   8  1  0.0040 -0.05  1.629000      1   0
    2  18  1  0.0030 -0.05  1.721135      1   0
    0  19  1  0.0100 -0.05  1.232333      1   1
    1  13  1  0.0025 -0.05  1.710726      1   1
    0  10  1  0.0090 -0.05 -0.743740      1   2
    1  29  1  0.0060 -0.05 -0.513952      1   2
    2  15  1  0.0035 -0.05 -0.283127      1   2
    0   1  1  0.0005 -0.05 -0.734266      1   3
    1  20  1  0.0065 -0.05 -0.004909      1   3
    2  33  1  0.0080 -0.05  1.460910      1   3
    3   1  1  0.0005 -0.05  1.698163      1   3


Summary statistics from samples
-------------------------------

The sub-module fwdpy.libseq (which we have imported as 'libseq') has a
function, 'summstats', which calculates many commonly-used summaries of
variation data.

.. code:: python

    ##This is an example of where you can do a lot in a 1-liner.
    ##We use nested list comprehensions to:
    ##  1. Get summary statistics for each element in samples.  We do neutral mutations (element 0)
    ##     and selected mutations (element 1) separately.
    ##  2. Turn each dict from libseq.summstats into a pandas.DataFrame
    ##  3. Combine all those DataFrame objects into one large DataFrame
    NeutralMutStats=pandas.concat([pandas.DataFrame(i.items(),columns=['stat','value']) 
                                   for i in [libseq.summstats(j[0]) for j in samples]])
    SelectedMutStats=pandas.concat([pandas.DataFrame(i.items(),columns=['stat','value'])
                                   for i in [libseq.summstats(j[1]) for j in samples]])
    print(NeutralMutStats)
    print(SelectedMutStats)


.. parsed-literal::

              stat      value
    0      thetapi  20.447368
    1  dsingletons  14.000000
    2            S  65.000000
    3         tajd   0.471357
    4   singletons  14.000000
    5       thetah  11.447368
    6       hprime   0.877725
    7       thetaw  18.321525
    0      thetapi  11.342105
    1  dsingletons  13.000000
    2            S  49.000000
    3         tajd  -0.718929
    4   singletons  13.000000
    5       thetah  12.763158
    6       hprime  -0.181022
    7       thetaw  13.811611
    0      thetapi  15.378947
    1  dsingletons  14.000000
    2            S  59.000000
    3         tajd  -0.304694
    4   singletons  14.000000
    5       thetah  15.778947
    6       hprime  -0.042769
    7       thetaw  16.630307
    0      thetapi  11.868421
    1  dsingletons  15.000000
    2            S  49.000000
    3         tajd  -0.565706
    4   singletons  15.000000
    5       thetah   6.973684
    6       hprime   0.623519
    7       thetaw  13.811611
              stat     value
    0      thetapi  0.300000
    1  dsingletons  3.000000
    2            S  3.000000
    3         tajd -1.723310
    4   singletons  3.000000
    5       thetah  0.015789
    6       hprime  0.347472
    7       thetaw  0.845609
    0      thetapi  0.200000
    1  dsingletons  2.000000
    2            S  2.000000
    3         tajd -1.512836
    4   singletons  2.000000
    5       thetah  0.010526
    6       hprime  0.299199
    7       thetaw  0.563739
    0      thetapi  0.300000
    1  dsingletons  3.000000
    2            S  3.000000
    3         tajd -1.723310
    4   singletons  3.000000
    5       thetah  0.015789
    6       hprime  0.347472
    7       thetaw  0.845609
    0      thetapi  0.400000
    1  dsingletons  4.000000
    2            S  4.000000
    3         tajd -1.867878
    4   singletons  4.000000
    5       thetah  0.021053
    6       hprime  0.382404
    7       thetaw  1.127478


The average :math:`\pi` under the model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Under the BGS model, the expectation of :math:`\pi` is
:math:`E[\pi]=\pi_0e^{-\frac{U}{2sh+r}},` :math:`U` is the mutation rate
to strongly-deleterious variants, :math:`\pi_0` is the value expected in
the absence of BGS (*i.e.* :math:`\pi_0 = \theta = 4N_e\mu`), :math:`s`
and :math:`h` are the selection and dominance coefficients, and
:math:`r` is the recombination rate.

Note that the definition of :math:`U` is *per diploid*, meaning twice
the per gamete rate. (See Hudson and Kaplan (1995) PMC1206891 for
details).

For our parameters, we have
:math:`E[\pi] = 20e^{-\frac{0.02}{0.1+0.005}},` which equals:

.. code:: python

    print(20*math.exp(-0.02/(0.1+0.005)))


.. parsed-literal::

    16.5313087525


Now, let's get the average :math:`\pi` from 500 simulated replicates. We
already have four replicates that we did above, so we'll run another 124
sets of four populations.

We will use standard Python to grow "pn", which is our list of
:math:`\pi` values calculated from neutral mutations from each
replicate.

.. code:: python

    for i in range(0,124,1):
        pops = fp.evolve_regions(rng,  
                             4,        
                             N,        
                             nlist[0:],
                             0.005,    
                             0.01,     
                             0.005,    
                             nregions, 
                             sregions, 
                             recregions)
        ##This is another heavy one-liner.
        ##We're taking samples of n=20 from each pop,
        ##Getting summstats for each neutral block from each sample,
        ##Turning the dict into pandas DataFrame objects,
        ##and returning a big DataFrame for all the data.
        temp = pandas.concat([pandas.DataFrame(i.items(),columns=['stat','value']) 
                              for i in [libseq.summstats(j[0]) for j in [fp.get_samples(rng,k,20) for k in pops]]])
        NeutralMutStats=pandas.concat([NeutralMutStats,temp])

Getting the mean diversity
^^^^^^^^^^^^^^^^^^^^^^^^^^

We've collected everything into a big pandas DataFrame. We can easily
get the mean using the built-in groupby and mean functions.

For users happier in R, you could write this DataFrame to a text file
and process it using R's
`dplyr <http://cran.r-project.org/web/packages/dplyr/index.html>`__
package, which is a really excellent tool for this sort of thing.

.. code:: python

    NeutralMutStats.groupby(['stat']).mean()




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>value</th>
        </tr>
        <tr>
          <th>stat</th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>S</th>
          <td>58.878000</td>
        </tr>
        <tr>
          <th>dsingletons</th>
          <td>17.016000</td>
        </tr>
        <tr>
          <th>hprime</th>
          <td>0.025761</td>
        </tr>
        <tr>
          <th>singletons</th>
          <td>17.902000</td>
        </tr>
        <tr>
          <th>tajd</th>
          <td>-0.099714</td>
        </tr>
        <tr>
          <th>thetah</th>
          <td>16.366926</td>
        </tr>
        <tr>
          <th>thetapi</th>
          <td>16.358968</td>
        </tr>
        <tr>
          <th>thetaw</th>
          <td>16.595919</td>
        </tr>
      </tbody>
    </table>
    </div>



The 'thetapi' record is our mean :math:`\pi` from all of the
simulations, and it is quite close to the theoretical value.

