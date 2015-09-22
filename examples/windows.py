
# coding: utf-8

# # Sliding windows
# 
# This is an example of running a simulation and getting a set of sliding windows from the output

# In[1]:

#import our modules
import fwdpy as fp
import fwdpy.libseq as lseq
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
windows = [fp.windows(i[0],0.1,0.1,0.) for i in samples]


# ### Summary stats from each window

# In[6]:

#For each window in each sample, get the basic summary statistics
stats = [[lseq.summstats(i) for i in j] for j in windows]


# Printing these outputs will be messy.  Let's just look at the first one. We'll help ourselves by printing out the window boundaries.
# 
# Note the limitation in the output:
# * Our simulation positions (defined above) are the half-open interval $[0,3)$.
# * In the output below, the last window is $[2.8,2.9)$.
# * The final empty window, $[2.9,3.0)$, is missing due to a limitation in libsequence that I'll need to fix soon.

# In[7]:


j=0
for i in np.arange(0,3,0.1):
    if j < len(stats[0]):
        print i," to ",(i+0.1),": ",stats[0][j]
    j=j+1

