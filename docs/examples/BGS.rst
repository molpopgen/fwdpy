
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
    import pandas as pd
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
    <type 'fwdpy.fwdpy.Spop'>


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
    details = [pd.DataFrame(fp.get_sample_details(i[1],j)) for i,j in zip(samples,pops)]

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

       dcount       ftime  generation    h  label  locus  origin      p     s
    0       1  4294967295       10000  1.0      0      0    9974  0.004 -0.05
       dcount       ftime  generation    h  label  locus  origin       p     s
    0       1  4294967295       10000  1.0      0      0    9976  0.0090 -0.05
    1       1  4294967295       10000  1.0      0      0    9939  0.0090 -0.05
    2       1  4294967295       10000  1.0      0      0    9957  0.0030 -0.05
    3       1  4294967295       10000  1.0      0      0    9991  0.0005 -0.05
    4       1  4294967295       10000  1.0      0      0    9946  0.0040 -0.05
    5       1  4294967295       10000  1.0      0      0    9996  0.0035 -0.05
       dcount       ftime  generation    h  label  locus  origin       p     s
    0       1  4294967295       10000  1.0      0      0    9984  0.0035 -0.05
    1       1  4294967295       10000  1.0      0      0    9975  0.0015 -0.05
    2       1  4294967295       10000  1.0      0      0    9983  0.0020 -0.05
    3       1  4294967295       10000  1.0      0      0    9998  0.0010 -0.05
    4       1  4294967295       10000  1.0      0      0    9940  0.0045 -0.05
       dcount       ftime  generation    h  label  locus  origin       p     s
    0       1  4294967295       10000  1.0      0      0    9990  0.0045 -0.05
    1       2  4294967295       10000  1.0      0      0    9953  0.0300 -0.05
    2       1  4294967295       10000  1.0      0      0    9968  0.0035 -0.05
    3       1  4294967295       10000  1.0      0      0    9987  0.0005 -0.05


.. code:: python

    #The order of the rows in each DataFrame is the
    #same as the order as the objects in 'samples':
    for i in range(4):
        print("Number of sites in samples[",i,"] = ",
              len(samples[i][1]),". Number of rows in DataFrame ",i,
              " = ",len(details[i].index),sep="")


.. parsed-literal::

    Number of sites in samples[0] = 1. Number of rows in DataFrame 0 = 1
    Number of sites in samples[1] = 6. Number of rows in DataFrame 1 = 6
    Number of sites in samples[2] = 5. Number of rows in DataFrame 2 = 5
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
    BigTable = pd.concat(details)
    
    print("This is the merged table:")
    print(BigTable)


.. parsed-literal::

    This is the merged table:
        dcount       ftime  generation    h  label  locus  origin       p     s  \
    0        1  4294967295       10000  1.0      0      0    9974  0.0040 -0.05   
    0        1  4294967295       10000  1.0      0      0    9976  0.0090 -0.05   
    1        1  4294967295       10000  1.0      0      0    9939  0.0090 -0.05   
    2        1  4294967295       10000  1.0      0      0    9957  0.0030 -0.05   
    3        1  4294967295       10000  1.0      0      0    9991  0.0005 -0.05   
    4        1  4294967295       10000  1.0      0      0    9946  0.0040 -0.05   
    5        1  4294967295       10000  1.0      0      0    9996  0.0035 -0.05   
    0        1  4294967295       10000  1.0      0      0    9984  0.0035 -0.05   
    1        1  4294967295       10000  1.0      0      0    9975  0.0015 -0.05   
    2        1  4294967295       10000  1.0      0      0    9983  0.0020 -0.05   
    3        1  4294967295       10000  1.0      0      0    9998  0.0010 -0.05   
    4        1  4294967295       10000  1.0      0      0    9940  0.0045 -0.05   
    0        1  4294967295       10000  1.0      0      0    9990  0.0045 -0.05   
    1        2  4294967295       10000  1.0      0      0    9953  0.0300 -0.05   
    2        1  4294967295       10000  1.0      0      0    9968  0.0035 -0.05   
    3        1  4294967295       10000  1.0      0      0    9987  0.0005 -0.05   
    0        1  4294967295       10000  1.0      0      0    9992  0.0010 -0.05   
    1        1  4294967295       10000  1.0      0      0    9908  0.0125 -0.05   
    2        1  4294967295       10000  1.0      0      0    9997  0.0010 -0.05   
    3        1  4294967295       10000  1.0      0      0    9999  0.0005 -0.05   
    4        1  4294967295       10000  1.0      0      0    9981  0.0065 -0.05   
    5        1  4294967295       10000  1.0      0      0    9997  0.0010 -0.05   
    0        1  4294967295       10000  1.0      0      0    9987  0.0025 -0.05   
    1        1  4294967295       10000  1.0      0      0    9864  0.0130 -0.05   
    2        1  4294967295       10000  1.0      0      0    9938  0.0035 -0.05   
    3        1  4294967295       10000  1.0      0      0    9990  0.0015 -0.05   
    4        1  4294967295       10000  1.0      0      0    9977  0.0075 -0.05   
    0        1  4294967295       10000  1.0      0      0    9983  0.0025 -0.05   
    1        1  4294967295       10000  1.0      0      0    9981  0.0120 -0.05   
    2        1  4294967295       10000  1.0      0      0    9993  0.0050 -0.05   
    ..     ...         ...         ...  ...    ...    ...     ...     ...   ...   
    0        1  4294967295       10000  1.0      0      0    9959  0.0060 -0.05   
    1        1  4294967295       10000  1.0      0      0    9980  0.0080 -0.05   
    2        1  4294967295       10000  1.0      0      0    9975  0.0040 -0.05   
    3        1  4294967295       10000  1.0      0      0    9976  0.0095 -0.05   
    4        1  4294967295       10000  1.0      0      0    9990  0.0005 -0.05   
    0        1  4294967295       10000  1.0      0      0    9998  0.0010 -0.05   
    1        1  4294967295       10000  1.0      0      0    9993  0.0030 -0.05   
    2        1  4294967295       10000  1.0      0      0    9958  0.0075 -0.05   
    0        1  4294967295       10000  1.0      0      0    9978  0.0035 -0.05   
    1        1  4294967295       10000  1.0      0      0    9959  0.0125 -0.05   
    2        1  4294967295       10000  1.0      0      0    9984  0.0065 -0.05   
    3        1  4294967295       10000  1.0      0      0    9985  0.0020 -0.05   
    4        1  4294967295       10000  1.0      0      0    9982  0.0065 -0.05   
    0        1  4294967295       10000  1.0      0      0    9988  0.0020 -0.05   
    1        1  4294967295       10000  1.0      0      0    9992  0.0005 -0.05   
    2        1  4294967295       10000  1.0      0      0    9988  0.0070 -0.05   
    3        1  4294967295       10000  1.0      0      0    9982  0.0020 -0.05   
    4        1  4294967295       10000  1.0      0      0    9975  0.0020 -0.05   
    5        1  4294967295       10000  1.0      0      0    9976  0.0075 -0.05   
    0        1  4294967295       10000  1.0      0      0    9991  0.0110 -0.05   
    1        1  4294967295       10000  1.0      0      0    9985  0.0040 -0.05   
    0        1  4294967295       10000  1.0      0      0    9996  0.0035 -0.05   
    1        1  4294967295       10000  1.0      0      0    9992  0.0060 -0.05   
    2        1  4294967295       10000  1.0      0      0    9998  0.0010 -0.05   
    3        1  4294967295       10000  1.0      0      0    9982  0.0165 -0.05   
    0        1  4294967295       10000  1.0      0      0    9978  0.0030 -0.05   
    1        1  4294967295       10000  1.0      0      0    9981  0.0035 -0.05   
    2        1  4294967295       10000  1.0      0      0    9993  0.0020 -0.05   
    3        2  4294967295       10000  1.0      0      0    9986  0.0095 -0.05   
    4        1  4294967295       10000  1.0      0      0    9998  0.0015 -0.05   
    
             pos  count  id  
    0   1.283749      1   0  
    0  -0.320125      1   1  
    1  -0.119514      1   1  
    2  -0.102664      1   1  
    3   1.042025      1   1  
    4   1.443235      1   1  
    5   1.804796      1   1  
    0  -0.632548      1   2  
    1  -0.460367      1   2  
    2  -0.119099      1   2  
    3  -0.055915      1   2  
    4   1.587909      1   2  
    0  -0.480663      1   3  
    1  -0.354161      2   3  
    2   1.351797      1   3  
    3   1.381058      1   3  
    0  -0.897565      1   4  
    1  -0.751551      1   4  
    2  -0.123021      1   4  
    3  -0.088541      1   4  
    4   1.242531      1   4  
    5   1.813220      1   4  
    0  -0.438263      1   5  
    1  -0.243305      1   5  
    2  -0.196773      1   5  
    3   1.483931      1   5  
    4   1.723836      1   5  
    0  -0.818220      1   6  
    1  -0.323753      1   6  
    2  -0.047837      1   6  
    ..       ...    ...  ..  
    0  -0.923395      1  33  
    1  -0.880657      1  33  
    2  -0.501588      1  33  
    3   1.434086      1  33  
    4   1.897750      1  33  
    0  -0.521129      1  34  
    1  -0.083328      1  34  
    2   1.777795      1  34  
    0  -0.703952      1  35  
    1  -0.385990      1  35  
    2   1.091824      1  35  
    3   1.257222      1  35  
    4   1.772582      1  35  
    0  -0.595274      1  36  
    1  -0.499114      1  36  
    2  -0.416441      1  36  
    3  -0.320007      1  36  
    4   1.178538      1  36  
    5   1.530433      1  36  
    0  -0.941568      1  37  
    1  -0.498442      1  37  
    0  -0.221834      1  38  
    1   1.537788      1  38  
    2   1.738061      1  38  
    3   1.994364      1  38  
    0  -0.983060      1  39  
    1  -0.475971      1  39  
    2  -0.105493      1  39  
    3   1.197116      2  39  
    4   1.257867      1  39  
    
    [143 rows x 12 columns]


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
    n = [polyt.SimData(i[0]) for i in samples]
    
    #Create "factories" for calculating the summary stats
    an = [sstats.PolySIM(i) for i in n]
    
    ##Collect a bunch of summary stats into a pandas.DataFrame:
    NeutralMutStats = pd.DataFrame([ {'thetapi':i.thetapi(),'npoly':i.numpoly(),'thetaw':i.thetaw()} for i in an ])
    
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
          <td>48</td>
          <td>12.636842</td>
          <td>13.529741</td>
        </tr>
        <tr>
          <th>1</th>
          <td>49</td>
          <td>12.563158</td>
          <td>13.811611</td>
        </tr>
        <tr>
          <th>2</th>
          <td>97</td>
          <td>35.968421</td>
          <td>27.341352</td>
        </tr>
        <tr>
          <th>3</th>
          <td>50</td>
          <td>12.921053</td>
          <td>14.093481</td>
        </tr>
        <tr>
          <th>4</th>
          <td>68</td>
          <td>17.015789</td>
          <td>19.167134</td>
        </tr>
        <tr>
          <th>5</th>
          <td>32</td>
          <td>8.584211</td>
          <td>9.019828</td>
        </tr>
        <tr>
          <th>6</th>
          <td>65</td>
          <td>18.905263</td>
          <td>18.321525</td>
        </tr>
        <tr>
          <th>7</th>
          <td>63</td>
          <td>20.036842</td>
          <td>17.757786</td>
        </tr>
        <tr>
          <th>8</th>
          <td>90</td>
          <td>28.968421</td>
          <td>25.368265</td>
        </tr>
        <tr>
          <th>9</th>
          <td>45</td>
          <td>12.452632</td>
          <td>12.684133</td>
        </tr>
        <tr>
          <th>10</th>
          <td>47</td>
          <td>12.242105</td>
          <td>13.247872</td>
        </tr>
        <tr>
          <th>11</th>
          <td>71</td>
          <td>20.752632</td>
          <td>20.012742</td>
        </tr>
        <tr>
          <th>12</th>
          <td>69</td>
          <td>19.900000</td>
          <td>19.449003</td>
        </tr>
        <tr>
          <th>13</th>
          <td>59</td>
          <td>10.900000</td>
          <td>16.630307</td>
        </tr>
        <tr>
          <th>14</th>
          <td>52</td>
          <td>14.494737</td>
          <td>14.657220</td>
        </tr>
        <tr>
          <th>15</th>
          <td>51</td>
          <td>16.563158</td>
          <td>14.375350</td>
        </tr>
        <tr>
          <th>16</th>
          <td>58</td>
          <td>14.500000</td>
          <td>16.348437</td>
        </tr>
        <tr>
          <th>17</th>
          <td>83</td>
          <td>29.478947</td>
          <td>23.395178</td>
        </tr>
        <tr>
          <th>18</th>
          <td>70</td>
          <td>25.689474</td>
          <td>19.730873</td>
        </tr>
        <tr>
          <th>19</th>
          <td>49</td>
          <td>15.473684</td>
          <td>13.811611</td>
        </tr>
        <tr>
          <th>20</th>
          <td>76</td>
          <td>23.847368</td>
          <td>21.422090</td>
        </tr>
        <tr>
          <th>21</th>
          <td>53</td>
          <td>7.784211</td>
          <td>14.939089</td>
        </tr>
        <tr>
          <th>22</th>
          <td>46</td>
          <td>14.315789</td>
          <td>12.966002</td>
        </tr>
        <tr>
          <th>23</th>
          <td>83</td>
          <td>27.100000</td>
          <td>23.395178</td>
        </tr>
        <tr>
          <th>24</th>
          <td>59</td>
          <td>17.157895</td>
          <td>16.630307</td>
        </tr>
        <tr>
          <th>25</th>
          <td>51</td>
          <td>14.778947</td>
          <td>14.375350</td>
        </tr>
        <tr>
          <th>26</th>
          <td>85</td>
          <td>20.984211</td>
          <td>23.958917</td>
        </tr>
        <tr>
          <th>27</th>
          <td>48</td>
          <td>12.915789</td>
          <td>13.529741</td>
        </tr>
        <tr>
          <th>28</th>
          <td>43</td>
          <td>9.763158</td>
          <td>12.120393</td>
        </tr>
        <tr>
          <th>29</th>
          <td>50</td>
          <td>12.063158</td>
          <td>14.093481</td>
        </tr>
        <tr>
          <th>30</th>
          <td>53</td>
          <td>8.421053</td>
          <td>14.939089</td>
        </tr>
        <tr>
          <th>31</th>
          <td>86</td>
          <td>26.073684</td>
          <td>24.240787</td>
        </tr>
        <tr>
          <th>32</th>
          <td>38</td>
          <td>10.889474</td>
          <td>10.711045</td>
        </tr>
        <tr>
          <th>33</th>
          <td>53</td>
          <td>13.236842</td>
          <td>14.939089</td>
        </tr>
        <tr>
          <th>34</th>
          <td>57</td>
          <td>9.842105</td>
          <td>16.066568</td>
        </tr>
        <tr>
          <th>35</th>
          <td>51</td>
          <td>11.831579</td>
          <td>14.375350</td>
        </tr>
        <tr>
          <th>36</th>
          <td>91</td>
          <td>29.989474</td>
          <td>25.650135</td>
        </tr>
        <tr>
          <th>37</th>
          <td>68</td>
          <td>18.752632</td>
          <td>19.167134</td>
        </tr>
        <tr>
          <th>38</th>
          <td>79</td>
          <td>23.105263</td>
          <td>22.267699</td>
        </tr>
        <tr>
          <th>39</th>
          <td>92</td>
          <td>27.031579</td>
          <td>25.932004</td>
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


Now, let's get the average :math:`\pi` from 1000 simulated replicates.
We already have 40 replicates that we did above, so we'll run another 24
sets of four populations.

We will use standard Python to grow our collection of summary
statistics.

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
        simdatasNeut = [polyt.SimData(i[0]) for i in samples]
        polySIMn = [sstats.PolySIM(i) for i in simdatasNeut]
        ##Append stats into our growing DataFrame:
        NeutralMutStats=pd.concat([NeutralMutStats,
                                       pd.DataFrame([ {'thetapi':i.thetapi(),
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

    npoly      58.372000
    thetapi    16.251858
    thetaw     16.453293
    dtype: float64



The 'thetapi' record is our mean :math:`\pi` from all of the
simulations, and it is quite close to the theoretical value.
