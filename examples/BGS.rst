
Example: background selection
=============================

Setting up the simulation
-------------------------

-  Neutral mutations will occur on the interval :math:`[0,1)`.
-  Strongly-deleterious mutations will occur on the intervals
   :math:`[-1,0)` and :math:`[1,2)`.
-  Recombination will be uniform throughout the region.

.. code:: python

    #Use Python 3's print a a function.
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

    #Simulate 40 replicate populations.  This uses C++11 threads behind the scenes:
    pops = fp.evolve_regions(rng,       #The random number generator 
                             40,         #The number of pops to simulate = number of threads to use.
                             N,         #Initial population size for each of the 40 demes
                             nlist[0:], #List of population sizes over time.
                             0.005,     #Neutral mutation rate (per gamete, per generation)
                             0.01,      #Deleterious mutation rate (per gamete, per generation)
                             0.005,     #Recombination rate (per diploid, per generation)
                             nregions,  #Defined above
                             sregions,  #Defined above
                             recregions)#Defined above

.. code:: python

    #Now, pops is a Python list with len(pops) = 40
    #Each element's type is fwdpy.singlepop
    print(len(pops))
    print(type(pops[0]))
                    


.. parsed-literal::

    40
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
    for i in samples[:4]:
        print ("A sample from a population is a ",type(i))
        
    print(len(samples))


.. parsed-literal::

    A sample from a population is a  <type 'tuple'>
    A sample from a population is a  <type 'tuple'>
    A sample from a population is a  <type 'tuple'>
    A sample from a population is a  <type 'tuple'>
    40


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
    #  label: A label applied for mutations for each region.  Here, I use 0 for all regions
    for i in details[:4]:
        print(i)


.. parsed-literal::

          a    h  label      p     s
    0  86.0  1.0      0  0.011 -0.05
          a    h  label       p     s
    0   6.0  1.0      0  0.0055 -0.05
    1  13.0  1.0      0  0.0080 -0.05
    2   8.0  1.0      0  0.0085 -0.05
          a    h  label       p     s
    0  39.0  1.0      0  0.0095 -0.05
    1   2.0  1.0      0  0.0020 -0.05
    2  15.0  1.0      0  0.0045 -0.05
    3   2.0  1.0      0  0.0005 -0.05
    4  29.0  1.0      0  0.0025 -0.05
          a    h  label       p     s
    0  24.0  1.0      0  0.0145 -0.05
    1   4.0  1.0      0  0.0030 -0.05


.. code:: python

    #The order of the rows in each DataFrame is the
    #same as the order as the objects in 'samples':
    for i in range(4):
        print("Number of sites in samples[",i,"] = ",
              len(samples[i][1]),". Number of rows in DataFrame ",i,
              " = ",len(details[i].index),sep="")


.. parsed-literal::

    Number of sites in samples[0] = 1. Number of rows in DataFrame 0 = 1
    Number of sites in samples[1] = 3. Number of rows in DataFrame 1 = 3
    Number of sites in samples[2] = 5. Number of rows in DataFrame 2 = 5
    Number of sites in samples[3] = 2. Number of rows in DataFrame 3 = 2


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
           a    h  label       p     s       pos  count  id
    0   86.0  1.0      0  0.0110 -0.05  1.984707      2   0
    0    6.0  1.0      0  0.0055 -0.05 -0.861297      1   1
    1   13.0  1.0      0  0.0080 -0.05 -0.201939      1   1
    2    8.0  1.0      0  0.0085 -0.05 -0.045890      1   1
    0   39.0  1.0      0  0.0095 -0.05 -0.092556      1   2
    1    2.0  1.0      0  0.0020 -0.05 -0.036264      1   2
    2   15.0  1.0      0  0.0045 -0.05  1.281776      1   2
    3    2.0  1.0      0  0.0005 -0.05  1.631626      1   2
    4   29.0  1.0      0  0.0025 -0.05  1.798943      1   2
    0   24.0  1.0      0  0.0145 -0.05 -0.996158      1   3
    1    4.0  1.0      0  0.0030 -0.05 -0.035800      1   3
    0    3.0  1.0      0  0.0015 -0.05 -0.344501      1   4
    1   66.0  1.0      0  0.0045 -0.05  1.178150      1   4
    2   20.0  1.0      0  0.0040 -0.05  1.894631      1   4
    0   10.0  1.0      0  0.0010 -0.05 -0.235474      1   5
    0   19.0  1.0      0  0.0010 -0.05  1.997798      1   6
    0   35.0  1.0      0  0.0060 -0.05 -0.248412      1   7
    0   50.0  1.0      0  0.0130 -0.05 -0.224335      1   8
    1   15.0  1.0      0  0.0060 -0.05 -0.178064      1   8
    2    6.0  1.0      0  0.0025 -0.05  1.668177      1   8
    3    2.0  1.0      0  0.0005 -0.05  1.718952      1   8
    4   83.0  1.0      0  0.0150 -0.05  1.962663      1   8
    0    8.0  1.0      0  0.0010 -0.05 -0.759213      2   9
    1    2.0  1.0      0  0.0015 -0.05 -0.223529      1   9
    2   14.0  1.0      0  0.0120 -0.05 -0.168920      1   9
    0   13.0  1.0      0  0.0025 -0.05 -0.992009      1  10
    1   17.0  1.0      0  0.0110 -0.05 -0.818799      1  10
    2   63.0  1.0      0  0.0025 -0.05 -0.203582      2  10
    3    1.0  1.0      0  0.0005 -0.05  1.267236      1  10
    0    7.0  1.0      0  0.0035 -0.05 -0.078749      1  11
    ..   ...  ...    ...     ...   ...       ...    ...  ..
    1   19.0  1.0      0  0.0095 -0.05 -0.057538      1  31
    2    8.0  1.0      0  0.0020 -0.05 -0.032058      1  31
    3   31.0  1.0      0  0.0030 -0.05  1.175901      1  31
    4    6.0  1.0      0  0.0025 -0.05  1.922342      1  31
    0    5.0  1.0      0  0.0010 -0.05  1.163766      1  32
    0    4.0  1.0      0  0.0010 -0.05 -0.976930      1  33
    1   14.0  1.0      0  0.0035 -0.05 -0.728838      1  33
    2    3.0  1.0      0  0.0010 -0.05 -0.573928      1  33
    3   29.0  1.0      0  0.0085 -0.05  1.261325      1  33
    0   18.0  1.0      0  0.0010 -0.05 -0.970838      1  34
    1   27.0  1.0      0  0.0070 -0.05 -0.504697      1  34
    2    2.0  1.0      0  0.0005 -0.05  1.580322      1  34
    3   35.0  1.0      0  0.0095 -0.05  1.901750      1  34
    0    5.0  1.0      0  0.0005 -0.05  1.290389      1  35
    0   12.0  1.0      0  0.0070 -0.05 -0.539382      1  36
    1    7.0  1.0      0  0.0050 -0.05  1.304475      1  36
    0    6.0  1.0      0  0.0015 -0.05 -0.840592      1  37
    1    4.0  1.0      0  0.0010 -0.05 -0.796569      1  37
    2    3.0  1.0      0  0.0025 -0.05  1.028634      1  37
    3   15.0  1.0      0  0.0030 -0.05  1.743109      1  37
    4   14.0  1.0      0  0.0055 -0.05  1.851309      1  37
    0   33.0  1.0      0  0.0045 -0.05 -0.861966      1  38
    1   33.0  1.0      0  0.0150 -0.05 -0.259942      1  38
    2   70.0  1.0      0  0.0090 -0.05 -0.009285      1  38
    3   16.0  1.0      0  0.0020 -0.05  1.966089      1  38
    0   22.0  1.0      0  0.0075 -0.05 -0.677911      1  39
    1   12.0  1.0      0  0.0120 -0.05 -0.379367      1  39
    2   17.0  1.0      0  0.0030 -0.05 -0.157381      1  39
    3    9.0  1.0      0  0.0070 -0.05  1.878176      1  39
    4   15.0  1.0      0  0.0025 -0.05  1.894812      1  39
    
    [129 rows x 8 columns]


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
          <td>50</td>
          <td>19.905263</td>
          <td>14.093481</td>
        </tr>
        <tr>
          <th>1</th>
          <td>54</td>
          <td>15.068421</td>
          <td>15.220959</td>
        </tr>
        <tr>
          <th>2</th>
          <td>75</td>
          <td>18.121053</td>
          <td>21.140221</td>
        </tr>
        <tr>
          <th>3</th>
          <td>72</td>
          <td>13.647368</td>
          <td>20.294612</td>
        </tr>
        <tr>
          <th>4</th>
          <td>71</td>
          <td>24.510526</td>
          <td>20.012742</td>
        </tr>
        <tr>
          <th>5</th>
          <td>41</td>
          <td>8.468421</td>
          <td>11.556654</td>
        </tr>
        <tr>
          <th>6</th>
          <td>78</td>
          <td>18.531579</td>
          <td>21.985830</td>
        </tr>
        <tr>
          <th>7</th>
          <td>17</td>
          <td>3.615789</td>
          <td>4.791783</td>
        </tr>
        <tr>
          <th>8</th>
          <td>87</td>
          <td>27.826316</td>
          <td>24.522656</td>
        </tr>
        <tr>
          <th>9</th>
          <td>59</td>
          <td>14.294737</td>
          <td>16.630307</td>
        </tr>
        <tr>
          <th>10</th>
          <td>65</td>
          <td>20.021053</td>
          <td>18.321525</td>
        </tr>
        <tr>
          <th>11</th>
          <td>63</td>
          <td>14.852632</td>
          <td>17.757786</td>
        </tr>
        <tr>
          <th>12</th>
          <td>78</td>
          <td>23.647368</td>
          <td>21.985830</td>
        </tr>
        <tr>
          <th>13</th>
          <td>64</td>
          <td>21.505263</td>
          <td>18.039655</td>
        </tr>
        <tr>
          <th>14</th>
          <td>68</td>
          <td>19.010526</td>
          <td>19.167134</td>
        </tr>
        <tr>
          <th>15</th>
          <td>52</td>
          <td>16.421053</td>
          <td>14.657220</td>
        </tr>
        <tr>
          <th>16</th>
          <td>67</td>
          <td>24.373684</td>
          <td>18.885264</td>
        </tr>
        <tr>
          <th>17</th>
          <td>54</td>
          <td>14.710526</td>
          <td>15.220959</td>
        </tr>
        <tr>
          <th>18</th>
          <td>62</td>
          <td>20.068421</td>
          <td>17.475916</td>
        </tr>
        <tr>
          <th>19</th>
          <td>55</td>
          <td>17.378947</td>
          <td>15.502829</td>
        </tr>
        <tr>
          <th>20</th>
          <td>53</td>
          <td>13.626316</td>
          <td>14.939089</td>
        </tr>
        <tr>
          <th>21</th>
          <td>49</td>
          <td>12.531579</td>
          <td>13.811611</td>
        </tr>
        <tr>
          <th>22</th>
          <td>55</td>
          <td>17.684211</td>
          <td>15.502829</td>
        </tr>
        <tr>
          <th>23</th>
          <td>49</td>
          <td>13.463158</td>
          <td>13.811611</td>
        </tr>
        <tr>
          <th>24</th>
          <td>44</td>
          <td>10.231579</td>
          <td>12.402263</td>
        </tr>
        <tr>
          <th>25</th>
          <td>80</td>
          <td>27.647368</td>
          <td>22.549569</td>
        </tr>
        <tr>
          <th>26</th>
          <td>61</td>
          <td>16.689474</td>
          <td>17.194046</td>
        </tr>
        <tr>
          <th>27</th>
          <td>81</td>
          <td>21.073684</td>
          <td>22.831439</td>
        </tr>
        <tr>
          <th>28</th>
          <td>69</td>
          <td>17.831579</td>
          <td>19.449003</td>
        </tr>
        <tr>
          <th>29</th>
          <td>54</td>
          <td>21.294737</td>
          <td>15.220959</td>
        </tr>
        <tr>
          <th>30</th>
          <td>38</td>
          <td>13.436842</td>
          <td>10.711045</td>
        </tr>
        <tr>
          <th>31</th>
          <td>48</td>
          <td>14.484211</td>
          <td>13.529741</td>
        </tr>
        <tr>
          <th>32</th>
          <td>54</td>
          <td>15.905263</td>
          <td>15.220959</td>
        </tr>
        <tr>
          <th>33</th>
          <td>73</td>
          <td>24.021053</td>
          <td>20.576482</td>
        </tr>
        <tr>
          <th>34</th>
          <td>67</td>
          <td>14.821053</td>
          <td>18.885264</td>
        </tr>
        <tr>
          <th>35</th>
          <td>45</td>
          <td>11.926316</td>
          <td>12.684133</td>
        </tr>
        <tr>
          <th>36</th>
          <td>66</td>
          <td>16.847368</td>
          <td>18.603394</td>
        </tr>
        <tr>
          <th>37</th>
          <td>72</td>
          <td>26.294737</td>
          <td>20.294612</td>
        </tr>
        <tr>
          <th>38</th>
          <td>68</td>
          <td>22.505263</td>
          <td>19.167134</td>
        </tr>
        <tr>
          <th>39</th>
          <td>68</td>
          <td>23.815789</td>
          <td>19.167134</td>
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


Now, let's get the average $\pi$ from 1000 simulated replicates.  We already have 40 replicates that we did above, so we'll run another 24 sets of four populations.  

We will use standard Python to grow our collection of summary statistics.

.. code:: python

    for i in range(0,24,1):
        pops = fp.evolve_regions(rng,  
                             40,        
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

    npoly      58.430000
    thetapi    16.405058
    thetaw     16.469641
    dtype: float64



The 'thetapi' record is our mean :math:`\pi` from all of the
simulations, and it is quite close to the theoretical value.
