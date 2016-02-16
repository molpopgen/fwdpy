import fwdpy as fp
import fwdpy.qtrait as qt
import numpy as np
import pandas
rng = fp.GSLrng(101)

mu=0.000625
r=0.5
sigE=0.1
hdf = pandas.HDFStore("qtrait_test_out.h5",'w',complevel=6,complib="zlib")
hdf.open()
k=0
#2,000 replicates in 4-thread chunks
for i in range(500):
    nlist = np.array([1000]*9000,dtype=np.uint32)
    #Evolve to equilibrium
    pops = qt.evolve_qtrait(rng,
                            4,
                            1000,
                            nlist[0:],
                            0,
                            mu,
                            r,
                            [],
                            [fp.GaussianS(0,1,1,0.25)],
                            [fp.Region(0,1,1)],
                            sigE,
                            0.)
    nlist = np.array([1000]*1000,dtype=np.uint32)
    st=qt.evolve_qtrait_popstats(rng,pops,
                          nlist[0:],
                          0,
                          mu,
                          r,
                          [],
                          [fp.GaussianS(0,1,1,0.25)],
                          [fp.Region(0,1,1)],
                          sigE,
                                 trackStats=1,
                          optimum=0.)
    for j in st:
        hdf.append('popstats',pandas.DataFrame(j))
    #Evolve another N gens after shift optimum to 0.5
    nlist = np.array([1000]*3000,dtype=np.uint32)
    st=qt.evolve_qtrait_popstats(rng,pops,
                                 nlist[0:],
                                 0,
                                 mu,
                                 r,
                                 [],
                                 [fp.GaussianS(0,1,1,0.25)],
                                 [fp.Region(0,1,1)],
                                 sigE,
                                 1,
                                optimum=0.5)
    for j in st:
        hdf.append('popstats',pandas.DataFrame(j))


hdf.close()
