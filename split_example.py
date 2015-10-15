#The first part is the same as the example for fwdpy.evolve_regions
import fwdpy
import numpy as np
nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
sregions = [fwdpy.ExpS(1,2,1,-0.001,0.0),fwdpy.ExpS(1,2,0.01,0.001)]
rregions = [fwdpy.Region(0,3,1)]
rng = fwdpy.GSLrng(100)
popsizes = np.array([1000],dtype=np.uint32)
# Evolve for 5N generations initially
popsizes=np.tile(popsizes,10000)
pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
#Now, "bud" off a daughter population of same size, and evolve both for another 100 generations
mpops = fwdpy.evolve_regions_split(rng,pops,popsizes[0:100],popsizes[0:100],0.001,0.0001,0.001,nregions,sregions,rregions,[0]*2)
#Sample deme 0
samples = [fwdpy.get_samples(rng,i,500,deme=0) for i in mpops]

print samples
gams = [fwdpy.view_gametes(i,0) for i in pops]

n=0
for i in gams[0]:
    n+=i['n']

print n


