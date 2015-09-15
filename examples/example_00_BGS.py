##Example of background selection (BGS)

#Neutral mutations will occur on the interval [0,1)
#Strongly-deleterious mutations will occur on the intervals [-1,0) and [1,2).
#Recombination will be uniform throughout the region.

from __future__ import print_function
import fwdpy as fp
import numpy as np

# Where neutral mutations occur:
nregions = [fp.Region(beg=0,end=1,weight=1)]
# Where selected mutations occur:
sregions = [fp.ConstantS(beg=-1,end=0,weight=1,s=-0.05,h=1),
            fp.ConstantS(beg=1,end=2,weight=1,s=-0.05,h=1)]
# Recombination:
recregions = [fp.Region(beg=-1,end=2,weight=1)]

#Population size
N=1000
#We'll evolve for 10N generations:
nlist = np.array([N]*10*N,dtype=np.uint32)

#Initalize a random number generator with seed value of 101
rng = fp.GSLrng(101)

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

#Now, pops is a Python list with len(pops) = 4
#Each element's type is fwdpy.singlepop
#for i in range(len(pops)):
#    print(type(pops[i]))
