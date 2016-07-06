
# coding: utf-8

# # Example of taking 'views' from simulated populations

# In[1]:

from __future__ import print_function
import fwdpy as fp
import pandas as pd
from background_selection_setup import *


# Get the mutations that are segregating in each population:

# In[2]:

mutations = [fp.view_mutations(i) for i in pops]


# Look at the raw data in the first element of each list:

# In[3]:

for i in mutations:
    print(i[0])


# Let's make that nicer, and convert each list of dictionaries to a Pandas DataFrame object:

# In[4]:

mutations2 = [pd.DataFrame(i) for i in mutations]


# In[5]:

for i in mutations2:
    print(i.head())


# The columns are:
# 
# * g = the generation when the mutation first arose
# * h = the dominance
# * n = the number of copies of the mutation in the population.  You can use this to get its frequency.
# * neutral = a boolean
# * pos = the position of the mutation
# * s = selection coefficient/effect size
# * label = The label assigned to a mutation.  These labels can be associated with Regions and Sregions.  Here, 1 is a mutation from the neutral region, 2 a selected mutation from the 'left' region and 3 a selected mutation from the 'right' regin.
# 
# We can do all the usual subsetting, etc., using regular pandas tricks.  For example, let's get the neutral mutations for each population:

# In[6]:

nmuts = [i[i.neutral == True] for i in mutations2]
for i in nmuts:
    print(i.head())


# We can also take views of gametes:

# In[7]:

gametes = [fp.view_gametes(i) for i in pops]


# The format is really ugly. v Each gamete is a dict with two elements:
# 
# * 'neutral' is a list of mutations _not_ affecting fitness.  The format is the same as for the mutation views above.
# * 'selected' is a list of mutations that _do_ affect fitness. The format is the same as for the mutation views above.

# In[8]:

for i in gametes:
    print(i[0])


# OK, let's clean that up.  We'll focus on the selected mutations for each individual, and turn everything into a pd.DataFrame.
# 
# We're only going to do this for the first simulated population.

# In[9]:

smuts = [i['selected'] for i in gametes[0]]


# We now have a list of lists stored in 'smuts'.

# In[10]:

smutsdf = pd.DataFrame()
ind=0
##Add the non-empty individuals to the df
for i in smuts:
    if len(i)>0:
        smutsdf = pd.concat([smutsdf,pd.DataFrame(i,index=[ind]*len(i))])
    ind += 1


# In[11]:

smutsdf.head()


# That's much better.  We can use the index to figure out which individual has which mutations, and their effect sizes, etc.
# 
# Finally, we can also take views of diploids.  Let's get the first two diploids in each population:

# In[12]:

dips = [fp.view_diploids(i,[0,1]) for i in pops]


# Again, the format here is ugly.  Each diploid view is a dictionary:

# In[13]:

for key in dips[0][0]:
    print(key)


# The values are:
# 
# * chrom0, chrom1 are gamete views, just like what we dealt with above
# * g = genetic component of phenotype
# * e = random component of phenotype
# * w = fitness
# * n0 and n1 are the number of selected variants on chrom0 and chrom1, respectively.
# * sh0 and sh1 are the sum of $s \times h$ for all selected mutations on chrom0 and chrom1, respectively
# 
# Please note that g, e, and w, may or may not be set by a particular simulation.  Their use is optional.
