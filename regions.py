import fwdpy
import numpy as np

nregions = [fwdpy.Region(0,1,1),
            fwdpy.Region(2,3,1)]
sregions = [fwdpy.ExpS(1,2,1,-0.1),
            fwdpy.ExpS(1,2,0.01,0.001)]
rregions = [fwdpy.Region(0,3,1)]

nr = fwdpy.process_regions(nregions)
sr = fwdpy.process_regions(sregions)
rr = fwdpy.process_regions(rregions)

#print type(sr)
rng = fwdpy.GSLrng(100)

#popsizes are NumPy arrays of 32bit unsigned integers
popsizes = np.array([1000],dtype=np.uint32)
popsizes = np.tile(popsizes,10*1000)

pops = fwdpy.evolve_regions(rng,3,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
s = [fwdpy.ms_sample(rng,i,10) for i in pops]

for i in range(len(s)):
    print s[i],"\n\n"

print [i.gen() for i in pops]

fwdpy.evolve_regions_more(rng,pops,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)

s = [fwdpy.ms_sample(rng,i,10) for i in pops]

for i in range(len(s)):
    print s[i],"\n\n"


print [i.gen() for i in pops]
print [i.sane() for i in pops]
print [i.popsize() for i in pops]


