import fwdpy as fp
import fwdpy.qtrait_mloc as qtm
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import fwdpy.fwdpyio as fpio
import sys
import pickle

N=1000
NLOCI=10
NREPS=40
rnge=fp.GSLrng(100)
rngs=fp.GSLrng(200)
nlist=np.array([N]*(1000),dtype=np.uint32)
theta_neutral_per_locus=0.0
rho_per_locus=100.0
little_r_per_locus=rho_per_locus/(4.0*float(N))
mu_n_region=theta_neutral_per_locus/(4.0*float(N))
mu_del_ttl=1e-3

f=qtm.MlocusAdditiveTrait()
pfile = open("pickle_out.p","wb")
for i in range(2):
    x = fp.MlocusPopVec(NREPS,N,NLOCI)
    sampler=fp.NothingSampler(len(x))

    qtm.evolve_qtraits_mloc_sample_fitness(rnge,x,sampler,f,
            nlist,[mu_n_region]*NLOCI,
            [mu_del_ttl/float(NLOCI)]*NLOCI,
            [fp.GaussianS(0,1,1,0.1)]*NLOCI,
            [little_r_per_locus]*NLOCI,
             [0.5]*(NLOCI-1),#loci unlinked
             1,optimum=0.5)
    print(fp.view_mutations(x[len(x)-1]))
    d=[fpio.serialize(i) for i in x]
    pickle.dump(d,pfile)
pfile.close()
pfile = open("pickle_out.p","rb")
for i in range(2):
    s=pickle.load(pfile)
    x=fpio.deserialize_mlocus(s)
    print (fp.view_mutations(x[len(x)-1]))
#for i,h in zip(s,details):
#    print(i)
#    print(h)
#
