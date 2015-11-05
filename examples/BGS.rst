
Example: background selection
=============================

Setting up the simulation
-------------------------

-  Neutral mutations will occur on the interval :math:`[0,1)`.
-  Strongly-deleterious mutations will occur on the intervals
   :math:`[-1,0)` and :math:`[1,2)`.
-  Recombination will be uniform throughout the region.

.. code:: python

    #Use Pyhon 3's print a a function.
    #This future-proofs the code in the notebook
    from __future__ import print_function
    #Import fwdpy.  Give it a shorter name
    import fwdpy as fp
    ##Other libs we need
    import numpy as np
    import pandas
    import math

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

        a  h      p     s
    0  67  1  0.006 -0.05
         a  h       p     s
    0   19  1  0.0025 -0.05
    1    6  1  0.0055 -0.05
    2   17  1  0.0120 -0.05
    3  149  1  0.0235 -0.05
    4   36  1  0.0070 -0.05
    5   24  1  0.0080 -0.05
        a  h       p     s
    0   3  1  0.0030 -0.05
    1   1  1  0.0005 -0.05
    2  19  1  0.0060 -0.05
        a  h       p     s
    0  24  1  0.0145 -0.05
    1   5  1  0.0010 -0.05
    2   2  1  0.0020 -0.05
    3  10  1  0.0105 -0.05
    4   6  1  0.0020 -0.05


.. code:: python

    #The order of the rows in each DataFrame is the
    #same as the order as the objects in 'samples':
    for i in range(len(samples)):
        print("Number of sites in samples[",i,"] = ",
              len(samples[i][1]),". Number of rows in DataFrame ",i,
              " = ",len(details[i].index),sep="")


.. parsed-literal::

    Number of sites in samples[0] = 1. Number of rows in DataFrame 0 = 1
    Number of sites in samples[1] = 6. Number of rows in DataFrame 1 = 6
    Number of sites in samples[2] = 3. Number of rows in DataFrame 2 = 3
    Number of sites in samples[3] = 5. Number of rows in DataFrame 3 = 5


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
    0   67  1  0.0060 -0.05  1.877543      1   0
    0   19  1  0.0025 -0.05 -0.928520      1   1
    1    6  1  0.0055 -0.05 -0.861297      1   1
    2   17  1  0.0120 -0.05 -0.843893      1   1
    3  149  1  0.0235 -0.05 -0.472551      1   1
    4   36  1  0.0070 -0.05 -0.426389      1   1
    5   24  1  0.0080 -0.05  1.162503      1   1
    0    3  1  0.0030 -0.05 -0.309919      1   2
    1    1  1  0.0005 -0.05  1.534583      1   2
    2   19  1  0.0060 -0.05  1.575055      1   2
    0   24  1  0.0145 -0.05 -0.996158      1   3
    1    5  1  0.0010 -0.05 -0.317365      1   3
    2    2  1  0.0020 -0.05 -0.124097      1   3
    3   10  1  0.0105 -0.05  1.342523      1   3
    4    6  1  0.0020 -0.05  1.515371      1   3


Summary statistics from samples
-------------------------------

We will use the `pylibseq <http://molpopgen.github.io/pylibseq/>`__
package to calculate summary statistics. pylibseq is a Python wrapper
around `libsequence <http://molpopgen.github.io/libsequence/>`__.

.. code:: python

    import libsequence.polytable as polyt
    import libsequence.summstats as sstats
    
    #Convert neutral mutations into libsequence "SimData" objects, 
    #which are intended to handle binary (0/1) data like
    #what comes out of these simulations
    n = [polyt.simData(i[0]) for i in samples]
    
    #Create "factories" for calculating the summary stats
    an = [sstats.polySIM(i) for i in n]
    
    ##Collect a bunch of summary stats into a pandas.DataFrame:
    NeutralMutStats = pandas.DataFrame([ {'thetapi':i.thetapi(),'npoly':i.numpoly(),'thetaw':i.thetaw()} for i in an ])
    
    NeutralMutStats




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>npoly</th>
          <th>thetapi</th>
          <th>thetaw</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>62</td>
          <td>18.405263</td>
          <td>17.475916</td>
        </tr>
        <tr>
          <th>1</th>
          <td>47</td>
          <td>12.763158</td>
          <td>13.247872</td>
        </tr>
        <tr>
          <th>2</th>
          <td>74</td>
          <td>19.310526</td>
          <td>20.858351</td>
        </tr>
        <tr>
          <th>3</th>
          <td>65</td>
          <td>13.205263</td>
          <td>18.321525</td>
        </tr>
      </tbody>
    </table>
    </div>



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


Now, let's get the average $\pi$ from 500 simulated replicates.  We already have four replicates that we did above, so we'll run another 124 sets of four populations.  

We will use standard Python to grow our collection of summary statistics.

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
        samples = [fp.get_samples(rng,i,20) for i in pops]
        simdatasNeut = [polyt.simData(i[0]) for i in samples]
        polySIMn = [sstats.polySIM(i) for i in simdatasNeut]
        ##Append stats into our growing DataFrame:
        NeutralMutStats=pandas.concat([NeutralMutStats,
                                       pandas.DataFrame([ {'thetapi':i.thetapi(),
                                                           'npoly':i.numpoly(),
                                                           'thetaw':i.thetaw()} for i in polySIMn ])])


Getting the mean diversity
^^^^^^^^^^^^^^^^^^^^^^^^^^

We've collected everything into a big pandas DataFrame. We can easily
get the mean using the built-in groupby and mean functions.

For users happier in R, you could write this DataFrame to a text file
and process it using R's
`dplyr <http://cran.r-project.org/web/packages/dplyr/index.html>`__
package, which is a really excellent tool for this sort of thing.

.. code:: python

    #Get means for each column:
    NeutralMutStats.mean(0)




.. parsed-literal::

    npoly      57.280000
    thetapi    15.743958
    thetaw     16.145491
    dtype: float64



The 'thetapi' record is our mean :math:`\pi` from all of the
simulations, and it is quite close to the theoretical value.
