
# coding: utf-8

# # Sliding windows
# 
# This is an example of running a simulation and getting a set of sliding windows from the output

# In[1]:

#import our modules
from __future__ import print_function
import fwdpy as fp
import fwdpy.libseq as lseq
import pandas
import numpy as np


# In[2]:

#set up our sim
rng = fp.GSLrng(101)
nregions = [fp.Region(0,1,1),fp.Region(2,3,1)]
sregions = [fp.ExpS(1,2,1,-0.1),fp.ExpS(1,2,0.1,0.001)]
rregions = [fp.Region(0,3,1)]
popsizes = np.array([1000]*10000,dtype=np.uint32)


# In[3]:

#Run the sim
pops = fp.evolve_regions(rng,4,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)


# In[4]:

#Take samples from the simulation
samples = [fp.get_samples(rng,i,20) for i in pops]


# ## Calculating sliding windows

# In[5]:

#For each of the neutral mutations in each sample, we will split
#the samples up into non-overlapping windows of size 0.1
windows = [lseq.windows(i[0],0.1,0.1,0.,3) for i in samples]


# ### Summary stats from each window

# In[6]:

#For each window in each sample, get the basic summary statistics
stats = [[lseq.summstats(i) for i in j] for j in windows]


# Printing these outputs will be messy as the output is a bunch of dict objects.  Let's merge all the output into a giant pandas.DataFrame for easier handling.

# In[7]:

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

