import fwdpy as fp
import fwdpy.qtrait_mloc as qtm
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

N=1000
NLOCI=10
NREPS=40
x = fp.popvec_mloc(NREPS,N,NLOCI)
print x[0].popsize()
print len(x)
rnge=fp.GSLrng(100)
rngs=fp.GSLrng(200)
nlist=np.array([N]*(10*N),dtype=np.uint32)
theta_neutral_per_locus=0.0
rho_per_locus=100.0
little_r_per_locus=rho_per_locus/(4.0*float(N))
mu_n_region=theta_neutral_per_locus/(4.0*float(N))
mu_del_ttl=1e-3

stats = qtm.evolve_qtraits_mloc_popstats(rnge,x,nlist,
                                         [mu_n_region]*NLOCI,
                                         [mu_del_ttl/float(NLOCI)]*NLOCI,
                                         [0.1]*NLOCI,
                                         [little_r_per_locus]*NLOCI,
                                         [0.5]*(NLOCI-1),#loci unlinked
                                         1)

df1 = [pd.DataFrame(i) for i in stats]
for i in range(len(df1)):
    df1[i]['rep']=[i]*len(df1[i].index)
    
stats = qtm.evolve_qtraits_mloc_popstats(rnge,x,nlist,
                                         [mu_n_region]*NLOCI,
                                         [mu_del_ttl/float(NLOCI)]*NLOCI,
                                         [0.1]*NLOCI,
                                         [little_r_per_locus]*NLOCI,
                                         [0.5]*(NLOCI-1),#loci unlinked
                                         1,optimum=0.5)


df2 = [pd.DataFrame(i) for i in stats]
for i in range(len(df2)):
    df2[i]['rep']=[i]*len(df2[i].index)

df = pd.concat([pd.concat(df1),pd.concat(df2)])

g = df.groupby(['stat','generation']).mean()
g.reset_index(inplace=True)
fig = plt.figure()
plt.plot(g[g.stat=='tbar'].generation,g[g.stat=='tbar'].value,label='Trait')
plt.plot(g[g.stat=='tbar'].generation,g[g.stat=='wbar'].value,label='Fitness')
plt.plot(g[g.stat=='tbar'].generation,g[g.stat=='varw'].value.multiply(100.0),label='100 x V(Fitness)')
plt.plot(g[g.stat=='tbar'].generation,g[g.stat=='VG'].value.multiply(10.0),label='10 x VG')
plt.xlabel("Time (generations)")
plt.ylabel("Mean value")
plt.legend(loc="upper left")
plt.savefig("mloc_test.png")
