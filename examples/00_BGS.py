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

#Now, pops is a Python list with len(pops) = 4
#Each element's type is fwdpy.singlepop
#for i in range(len(pops)):
#    print(type(pops[i]))
