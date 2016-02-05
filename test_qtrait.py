import fwdpy as fp
import fwdpy.qtrait as qt
import numpy as np
import pandas
rng = fp.GSLrng(101)

mu=0.000625
r=0.5
sigE=0.1
hdf = pandas.HDFStore("foo2.h5",'w')
hdf.open()
k=0
for i in range(100):
    nlist = np.array([1000]*9000,dtype=np.uint32)
    #Evolve to equilibrium
    pops = qt.evolve_qtrait(rng,
                            1,
                            1000,
                            nlist[0:],
                            0,
                            mu,
                            r,
                            [],
                            [fp.GaussianS(0,1,1,0.25)],
                            [fp.Region(0,1,1)],
                            sigE,
                            0.,
                            track=0,trackStats=0)
    nlist = np.array([1000]*1000,dtype=np.uint32)
    qt.evolve_qtrait_more(rng,pops,
                          nlist[0:],
                          0,
                          mu,
                          r,
                          [],
                          [fp.GaussianS(0,1,1,0.25)],
                          [fp.Region(0,1,1)],
                          sigE,
                          0.,
                          track=0,trackStats=1)
    #Evolve another N gens after shift optimum to 0.5
    nlist = np.array([1000]*3000,dtype=np.uint32)
    qt.evolve_qtrait_more(rng,pops,
                          nlist[0:],
                          0,
                          mu,
                          r,
                          [],
                          [fp.GaussianS(0,1,1,0.25)],
                          [fp.Region(0,1,1)],
                          sigE,
                          optimum=0.5,
                          track=0,trackStats=1)
    st=[qt.popstats(i) for i in pops]
    for i in st:
        hdf.append('popstats',pandas.DataFrame(i))

hdf.close()
