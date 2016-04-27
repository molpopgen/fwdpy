import fwdpy as fp
import fwdpy.qtrait_mloc as qtm
import numpy as np

N=1000
NLOCI=10
NREPS=20
x = fp.popvec_mloc(NREPS,N,NLOCI)
print x[0].popsize()
print len(x)
rnge=fp.GSLrng(100)
rngs=fp.GSLrng(200)
nlist=np.array([N]*(10*N),dtype=np.uint32)
mu_n_region=1e-3
mu_del_ttl=1e-3
print mu_del_ttl/float(NLOCI)
samples = qtm.evolve_qtraits_mloc_sample(rnge,rngs,x,nlist,
                                         [mu_n_region]*NLOCI,
                                         [mu_del_ttl/float(NLOCI)]*NLOCI,
                                         [0.1]*NLOCI,
                                         [mu_n_region]*NLOCI,
                                         [0.5]*(NLOCI-1),#loci unlinked
                                         0,0,0.,1,10,20)
print samples[0]
