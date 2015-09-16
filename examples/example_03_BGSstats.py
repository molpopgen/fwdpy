#Summary statistics from population samples

from __future__ import print_function
from example_02_BGSdetails import *

#Currently, there is limited support for summary statistic calculation in fwdpy.
#KRT is still deciding whether or not to put it here or to make a lightweight
#Cython interface to the relevant part of libsequence.

#So far, what we can do is calculate Tajima's D:

DvalsNeutral = [fp.TajimasD(i[0]) for i in samples]
DvalsSelected = [fp.TajimasD(i[1]) for i in samples]

print("The distribution of D for neutral variants:",DvalsNeutral)
print("The distribution of D for selected variants:",DvalsSelected)

#We don't need to wait for KRT to help us, though.
#Summary stats are easy!

###fxn to ask if a site is a derived singleton
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

##fxn to calculate "pi" at a site
def pisite( site ):
    ones=site[1].count('1')
    nsam=len(site[1])
    p=float(ones)/float(nsam)
    q=float(ones-1)/float(nsam-1)
    return 2.*p*q

##Calculate pi for whole sample
def pisample(sample):
    pivals = [pisite(i) for i in sample]
    return sum(pivals)

##Apply our new functions:

##No. singletons at neutral sites
nsing = [countDerived(i[0]) for i in samples]
#No. singletons at selected sites
ssing = [countDerived(i[1]) for i in samples]
##Pi at neutral sites
pn = [pisample(i[0]) for i in samples]
##Pi at selected sites
ps = [pisample(i[1]) for i in samples]

print("Dist. of number of derived neutral singletons is:",nsing)
print("Dist. of number of derived selected singletons is:",ssing)
print("Dist. of pi at neutral mutations is:",pn)
print("Dist. of pi at selected mutations is:",ps)
