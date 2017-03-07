
.. _BGS_multiprocessing:

Advanced example: background selection (revisited)
==================================================

This is the same background selection simulation as in the previous
example, but with the following change to the implementation details:

-  We change the nature of the parallelism. The previous example uses
   fwdpy to run 40 simulations at a time, process them, and then repeat
   the process 25 times, doing all of the analysis in-memory. Here, we
   use the multiprocessing module to spawn 40 separate Python processes.
   Each process runs 25 simulations and records the summary statistics.
   At the end of the 25 replicates, the data are written to an SQLite
   database and get the mean values via an SQL query, which is
   out-of-memory.

The purpose of this example is to show that there are multiple ways to
do things in terms of how to use parallel processing to perform
simulations. Further, the technique of writing results to an SQLite
database is very powerful as it allows many analyses ("aggregations") to
be done without loading all of your simulation results into RAM.

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
    import os
    import sqlite3
    import multiprocessing as mp
    import libsequence.polytable as polyt
    import libsequence.summstats as sstats

Define the function that we will run in separate Python processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The details of setting up the simulation are identical to the prevous
BGS example.

.. code:: python

    def simulate_async(args):
        """
        This function will be run in a separate process
        using the multiprocessing module.  Its argument 
        list is a tuple.
    
        """
        #Assign names to the tuple values
        seed,dbname,tablename = args
        
        # Where neutral mutations occur:
        nregions = [fp.Region(beg=0,end=1,weight=1)]
    
        # Where selected mutations occur:
        sregions = [fp.ConstantS(beg=-1,end=0,weight=1,s=-0.05,h=1),
                    fp.ConstantS(beg=1,end=2,weight=1,s=-0.05,h=1)]
    
        # Recombination:
        recregions = [fp.Region(beg=-1,end=2,weight=1)]
    
        #Population size
        N=1000
        #We'll evolve for 10N generations.
        #nlist is a list of population sizes over time.
        #len(nlist) is the length of the simulation
        #We use numpy arrays for speed and optimised RAM
        #use.  Note the dtype=np.uint32, which means 32-bit
        #unsigned integer. Failure to use this type will
        #cause a run-time error.
        nlist = np.array([N]*(10*N),dtype=np.uint32)
    
        #Initalize a random number generator with seed value of 101
        rng = fp.GSLrng(seed)
    
        summstats=[]
        for replicate in range(0,25,1):
            pops = fp.evolve_regions(rng,  
                                 1,       #Simulate only 1 population at a time     
                                 N,        
                                 nlist[0:],
                                 0.005,    
                                 0.01,     
                                 0.005,    
                                 nregions, 
                                 sregions, 
                                 recregions)
            sample = fp.get_samples(rng,pops[0],20)
            simdatasNeut = polyt.SimData(sample[0])
            polySIMn = sstats.PolySIM(simdatasNeut)
            ##Append stats into our growing DataFrame:
            summstats.append({'thetapi':polySIMn.thetapi(),'npoly':polySIMn.numpoly(),'thetaw':polySIMn.thetaw()})
        DF=pd.DataFrame(summstats)
        con = sqlite3.connect(dbname)
        DF.to_sql(tablename,con,if_exists='append',index=False)
        con.close()


Run the simulations
===================

The following block of code sets up a thread pool to run the above
function using 40 separate processes.

.. code:: python

    if os.path.isfile('BGSmp.db'):
        os.remove('BGSmp.db')
    np.random.seed(101)
    args=[(seed,'BGSmp.db','stats') for seed in np.random.randint(0,42000000,40)]
    #P a thread pool using the number of processors on your machine
    #If you have < 40 cores, it'll spawn new processes as old ones finish.
    #for i in args: simulate_async(i)
    P=mp.Pool() 
    P.imap_unordered(simulate_async,args)
    P.close()
    P.join()

Getting the mean diversity
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    #open database connection:
    c=sqlite3.connect('BGSmp.db')
    #Get means for each column:
    x=pd.read_sql_query('select avg(npoly),avg(thetapi),avg(thetaw) from stats',c)
    c.close()
    os.remove('BGSmp.db')
    print(x)


.. parsed-literal::

       avg(npoly)  avg(thetapi)  avg(thetaw)
    0      57.635     16.033353    16.245555


The 'thetapi' record is our mean :math:`\pi` from all of the
simulations, and it is quite close to the theoretical value.
