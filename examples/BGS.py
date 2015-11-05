
# coding: utf-8

# # Example: background selection
# 

# ## Setting up the simulation
# * Neutral mutations will occur on the interval $[0,1)$.
# * Strongly-deleterious mutations will occur on the intervals $[-1,0)$ and $[1,2)$.
# * Recombination will be uniform throughout the region.

# In[20]:

#Use Pyhon 3's print a a function.
#This future-proofs the code in the notebook
from __future__ import print_function
#Import fwdpy.  Give it a shorter name
import fwdpy as fp
##Other libs we need
import numpy as np
import pandas
import math


# ### Establishing 'regions' for mutation and recombination
# 

# In[3]:

# Where neutral mutations occur:
nregions = [fp.Region(beg=0,end=1,weight=1)]


# In[4]:

# Where selected mutations occur:
sregions = [fp.ConstantS(beg=-1,end=0,weight=1,s=-0.05,h=1),
            fp.ConstantS(beg=1,end=2,weight=1,s=-0.05,h=1)]


# In[5]:

# Recombination:
recregions = [fp.Region(beg=-1,end=2,weight=1)]


# ### Population size and simulation length

# In[6]:

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


# In[7]:

#Initalize a random number generator with seed value of 101
rng = fp.GSLrng(101)


# In[8]:

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


# In[9]:

#Now, pops is a Python list with len(pops) = 4
#Each element's type is fwdpy.singlepop
print(len(pops))
for i in range(len(pops)):
    print(type(pops[i]))
                


# ## Taking samples from simulated populations

# In[10]:

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


# ### Getting additional information about samples

# In[11]:

#Again, use list comprehension to get the 'details' of each sample
#Given that each object in samples is a tuple, and that the second
#item in each tuple represents selected mutations, i[1] in the line
#below means that we are getting the mutation information only for
#selected variants
details = [fp.get_sample_details(i[1],j) for i,j in zip(samples,pops)]


# In[12]:

#details is now a list of pandas DataFrame objects
#Each DataFrame has the following columns:
#  a: mutation age (in generations)
#  h: dominance of the mutation
#  p: frequency of the mutation in the population
#  s: selection coefficient of the mutation
for i in details:
    print(i)


# In[13]:

#The order of the rows in each DataFrame is the
#same as the order as the objects in 'samples':
for i in range(len(samples)):
    print("Number of sites in samples[",i,"] = ",
          len(samples[i][1]),". Number of rows in DataFrame ",i,
          " = ",len(details[i].index),sep="")


# In[14]:

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


# In[15]:

##Merge into 1 big DataFrame:
BigTable = pandas.concat(details)

print("This is the merged table:")
print(BigTable)


# ## Summary statistics from samples
# 
# We will use the [pylibseq](http://molpopgen.github.io/pylibseq/) package to calculate summary statistics.  pylibseq is a Python wrapper around [libsequence](http://molpopgen.github.io/libsequence/).

# In[16]:

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


# ### The average $\pi$ under the model
# 
# Under the BGS model, the expectation of $\pi$ is $E[\pi]=\pi_0e^{-\frac{U}{2sh+r}},$ $U$ is the mutation rate to strongly-deleterious variants, $\pi_0$ is the value expected in the absence of BGS (_i.e._ $\pi_0 = \theta = 4N_e\mu$), $s$ and $h$ are the selection and dominance coefficients, and $r$ is the recombination rate.
# 
# Note that the definition of $U$ is _per diploid_, meaning twice the per gamete rate. (See Hudson and Kaplan (1995) PMC1206891 for details).
# 
# For our parameters, we have $E[\pi] = 20e^{-\frac{0.02}{0.1+0.005}},$ which equals:

# In[17]:

print(20*math.exp(-0.02/(0.1+0.005)))

Now, let's get the average $\pi$ from 500 simulated replicates.  We already have four replicates that we did above, so we'll run another 124 sets of four populations.  

We will use standard Python to grow our collection of summary statistics.
# In[18]:

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


# #### Getting the mean diversity
# We've collected everything into a big pandas DataFrame.  We can easily get the mean using the built-in groupby and mean functions.  
# 
# For users happier in R, you could write this DataFrame to a text file and process it using R's [dplyr](http://cran.r-project.org/web/packages/dplyr/index.html) package, which is a really excellent tool for this sort of thing.

# In[19]:

#Get means for each column:
NeutralMutStats.mean(0)


# The 'thetapi' record is our mean $\pi$ from all of the simulations, and it is quite close to the theoretical value. 
