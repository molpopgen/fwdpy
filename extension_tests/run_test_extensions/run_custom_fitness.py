import fwdpy as fp
import extension_tests.test_fwdpy_extensions.test_custom_fitness as tfp
import numpy as np

rng = fp.GSLrng(101)
rngs=fp.GSLrng(202)
p = fp.popvec(3,1000)
s = fp.nothing_sampler(len(p))

n=np.array([1000]*10000,dtype=np.uint32)
nr=[fp.Region(0,1,1)]
sr=[fp.ExpS(0,1,1,0.1)]

fitness = tfp.additiveFitnessTesting()
fp.evolve_regions_sampler_fitness(rng,p,s,fitness,n,0.001,0.001,0.001,nr,sr,nr,1)

fitness = tfp.Aa_only_testing()
fp.evolve_regions_sampler_fitness(rng,p,s,fitness,n,0.001,0.001,0.001,nr,sr,nr,1)

