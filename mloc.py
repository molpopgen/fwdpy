import fwdpy as fp
import fwdpy.qtrait_mloc as qtm
import fwdpy.demography as fpd
import numpy as np
import pandas as pd
import fwdpy.fwdpyio as fpio
import sys
import gzip
import pickle

def get_nlist():
    """
    Generates a numpy array of the canges in N over time

    There are 5 epochs, with t=0 being the present.
    
    E1: Ne= 7,310  from t=start(8N?) to t = - 5920 generation (Ancestral sizes until 5920 generations ago)
    
    E2: Ne =14,474 from t = -5920 to t = -2040 (Ancient growth at 5920 g ago)
    
    E3: Ne =1,861 from t= -2040 to t= -920 (OOA, first bottle neck 2040 g ago)
    
    E4: Ne = 1,032 to Ne = 9,300 during t = -920 to t = -205 ( second bottle neck and onset of 715 g of exponential growth at rate 0.31% per gen )  
    
    E5: Ne = 9,300 to Ne = 512,000 during t = -205 to t = -0 ( 205 g of exponential growth at rate 1.95% per gen )  
    """
    n=[7310]*(10*7310) #E1: evolve ancestral size to mutation/selection/drift equilibrium
    n.extend([14474]*(5920-2040)) #E2
    n.extend([1861]*(2040-920)) #E3
    n.extend(fpd.exponential_size_change(1032,9300,920-205)) #E4
    n.extend(fpd.exponential_size_change(9300,51200,205)) #E5
    return n

NLOCI=10
NREPS=10
rnge=fp.GSLrng(100)
nlist=np.array(get_nlist(),dtype=np.uint32)
theta_neutral=100.0
rho_per_locus=100.0
f=qtm.MlocusAdditiveTrait()
delmuts=[1e-3]*NLOCI
mmodels=[fp.GaussianS(0,1,1,0.25)]*len(delmuts)
for i in mmodels:
    print i
nmuts=[theta_neutral/float(4*nlist[0])]*NLOCI
recrates=[theta_neutral/float(4*nlist[0])]*NLOCI

pfile = gzip.open("pickle_out.p","wb")

for i in range(1):
    x = fp.MlocusPopVec(NREPS,nlist[0],NLOCI)
    sampler=fp.NothingSampler(len(x))

    qtm.evolve_qtraits_mloc_sample_fitness(rnge,x,sampler,f,
            nlist[:5000],
            nmuts,
            delmuts,
            mmodels,
            recrates,
            [0.0]*(NLOCI-1),sample=0,
            VS=2,optimum=0.)
    v = fp.view_mutations(x)
    for i in v:
        print i
    for i in x:
        d=fpio.serialize(i)
        pickle.dump(d,pfile)
pfile.close()
