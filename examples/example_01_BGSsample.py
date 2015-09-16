##Sample from the populations simulated in 00_BGS.py

##Run the other script
from __future__ import print_function
from example_00_BGS import *


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
