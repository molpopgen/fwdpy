#The next three lines process and compile our custom fitness module:
import pyximport
pyximport.install()
import test_fwdpy_extensions.test_custom_fitness as tfp
#import fwdpy and numpy as usual
import fwdpy as fp
import numpy as np

rng = fp.GSLrng(101)
rngs=fp.GSLrng(202)
p = fp.SpopVec(3,1000)
s = fp.NothingSampler(len(p))

n=np.array([1000]*1000,dtype=np.uint32)
nr=[fp.Region(0,1,1)]
sr=[fp.ExpS(0,1,1,0.1)]

#Now, let's do some evolution with our 'custom' fitness functions:
fitness = tfp.AdditiveFitnessTesting()
fp.evolve_regions_sampler_fitness(rng,p,s,fitness,n,0.001,0.001,0.001,nr,sr,nr,1)

fitness = tfp.AaOnlyTesting()
fp.evolve_regions_sampler_fitness(rng,p,s,fitness,n,0.001,0.001,0.001,nr,sr,nr,1)

fitness = tfp.GBRFitness()
fp.evolve_regions_sampler_fitness(rng,p,s,fitness,n,0.001,0.001,0.001,nr,sr,nr,1)

