##Get the 'details' for a sample

#Import results from previous example
from __future__ import print_function
from example_01_BGSsample import *
import pandas

#Again, use list comprehension to get the 'details' of each sample
#Given that each object in samples is a tuple, and that the second
#item in each tuple represents selected mutations, i[1] in the line
#below means that we are getting the mutation information only for
#selected variants
details = [fp.get_sample_details(i[1],j) for i,j in zip(samples,pops)]

#details is now a list of pandas DataFrame objects
#Each DataFrame has the following columns:
#  a: mutation age (in generations)
#  h: dominance of the mutation
#  p: frequency of the mutation in the population
#  s: selection coefficient of the mutation
for i in details:
    print(i)

#The order of the rows in each DataFrame is the
#same as the order as the objects in 'samples':
for i in range(len(samples)):
    print("Number of sites in samples[",i,"] = ",len(samples[i][1]),". Number of rows in DataFrame ",i," = ",len(details[i].index),sep="")

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

##Merge into 1 big DataFrame:
BigTable = pandas.concat(details)

print("This is the merged table:")
print(BigTable)

#Pandas has some "dplyr-like" functionality.
#Users may prefer to do downstream analysis in R using dplyr, etc.
#If so, they can write the DataFrames to files that can be read
#into R as data.frames:
#BigTable.to_csv("BGSTable.csv",sep="\t",index=False)

#But, we'll do some simple stuff here:
print("The mean frequency of a deleterious allele in the population is: ",BigTable['p'].mean(),'.',sep="")
print("The mean age of a deleterious allele in the population is:",BigTable['a'].mean(),"generations.")

