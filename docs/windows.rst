
Example: Sliding windows
========================

This is an example of running a simulation and getting a set of sliding
windows from the output

.. code:: python

    #import our modules
    from __future__ import print_function
    import fwdpy as fp
    import fwdpy.libseq as lseq
    import pandas
    import numpy as np
    import datetime
    import time

.. code:: python

    ##Info
    dt=datetime.datetime.now()
    print("This example was processed using ",fp.pkg_version(), "on",dt.month,"/",dt.day,"/",dt.year)
    print("The dependency versions are",fp.pkg_dependencies())



.. parsed-literal::

    This example was processed using  {'fwdpy': '0.0.1'} on 9 / 23 / 2015
    The dependency versions are {'libsequence': '1.8.7', 'GSL': '1.16', 'fwdpp': '0.3.8'}


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

.. code:: python

    #For each of the neutral mutations in each sample, we will split
    #the samples up into non-overlapping windows of size 0.1
    windows = [lseq.windows(i[0],0.1,0.1,0.,3) for i in samples]

Summary stats from each window
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    #For each window in each sample, get the basic summary statistics
    stats = [[lseq.summstats(i) for i in j] for j in windows]

Printing these outputs will be messy as the output is a bunch of dict
objects. Let's merge all the output into a giant pandas.DataFrame for
easier handling.

.. code:: python

    allstats=pandas.DataFrame()
    starts = np.arange(0.,3.,0.1)
    stops = starts + 0.1
    
    for i in range(len(stats)):
        temp = pandas.DataFrame.from_dict(stats[i])
        temp['replicate']=[i]*len(temp.index)
        temp['starts']=starts
        temp['stops']=stops
        allstats=pandas.concat([allstats,temp])
    
    #Now, that's cleaner!
    print (allstats.head())
    print (allstats.tail())


.. parsed-literal::

       S  dsingletons    hprime  singletons      tajd    thetah   thetapi  \
    0  0            0       NaN           0       NaN  0.000000  0.000000   
    1  0            0       NaN           0       NaN  0.000000  0.000000   
    2  0            0       NaN           0       NaN  0.000000  0.000000   
    3  1            0 -1.871112           0  0.722614  1.184211  0.394737   
    4  0            0       NaN           0       NaN  0.000000  0.000000   
    
        thetaw  replicate  starts  stops  
    0  0.00000          0     0.0    0.1  
    1  0.00000          0     0.1    0.2  
    2  0.00000          0     0.2    0.3  
    3  0.28187          0     0.3    0.4  
    4  0.00000          0     0.4    0.5  
        S  dsingletons    hprime  singletons      tajd    thetah   thetapi  \
    25  0            0       NaN           0       NaN  0.000000  0.000000   
    26  2            1  0.498665           1  0.063253  0.263158  0.578947   
    27  0            0       NaN           0       NaN  0.000000  0.000000   
    28  1            0 -1.397097           0  1.025883  1.031579  0.442105   
    29  1            1  0.224533           1 -1.164391  0.005263  0.100000   
    
          thetaw  replicate  starts  stops  
    25  0.000000          3     2.5    2.6  
    26  0.563739          3     2.6    2.7  
    27  0.000000          3     2.7    2.8  
    28  0.281870          3     2.8    2.9  
    29  0.281870          3     2.9    3.0  

