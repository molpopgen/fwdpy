import fwdpy
import pandas

rng = fwdpy.GSLrng(100)
pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*int(1e4),50,50)
s = [fwdpy.get_samples(rng,i,[10,]) for i in pop]

###fxn to ask if a site is a singleton
def isSingleton( site ):
    ones=site[1].count('1')
    if ones == 1:
        return True
    return False

##fxn to count derived singletonss in a sample
def countDerived( sample ):
    nsing=0
    for i in range(len(sample)):
        if isSingleton(sample[i]):
            nsing += 1
    return nsing

#list comprehension for automatic vectorizing
D = [fwdpy.TajimasD(si[0]) for si in s]
print "Tajima's D per sample =", D

#number of seg sites per sample
segsites = [len(si[0]) for si in s]
print "S per sample =",segsites

#number of singletons per sample
nsing = [countDerived(si[0]) for si in s]
print "No. singletons per sample =",nsing

## remove the non-singleton sites from each sample
sAllSing = [filter( lambda x: isSingleton(x) == True, j[0] ) for j in s]
print "No. singletons per sample =",[len(i) for i in sAllSing]

##Remove the singletons
sNoSing = [filter( lambda x: isSingleton(x) == False, j[0] ) for j in s]
print "No. non-singletons per sample =",[len(i) for i in sNoSing]

##Get xtra info for each mutation in the sample
sh = [fwdpy.get_sample_details(i[0],j) for i,j in zip(s,pop)]

##Add a column to each DataFrame specifying the mutation position, count of derived state, and a "replicate ID"
for i in range(len(sh)):
    sh[i]['pos']=[x[0] for x in s[i][0]]
    sh[i]['freq']=[ x[1].count('1') for x in s[i][0]]
    sh[i]['id']=[i]*len(sh[i].index)

##Write all DataFrames to a file
pandas.concat(sh).to_csv("test.csv",sep="\t",index=False)

##Now, evolve them some more and end with a bottleneck + recent, partial recovery
fwdpy.evolve_pops_more_t(rng,pop,[1000]*int(1e4) + [500]*100 + [750]*10,50,50)

##Check that all is cool with the data structures...
for i in range(len(pop)):
    print pop[i].gen()," ",pop[i].popsize()," ",pop[i].sane()
