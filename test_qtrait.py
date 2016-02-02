import fwdpy as fp
import fwdpy.qtrait as qt
import numpy as np
import pandas
rng = fp.GSLrng(101)


hdf = pandas.HDFStore("foo2.h5",'w')
hdf.open()
k=0
for i in range(250):
    nlist = np.array([1000]*9000,dtype=np.uint32)
    #Evolve to equilibrium
    pops = qt.evolve_qtrait(rng,
                            4,
                            1000,
                            nlist[0:],
                            0,
                            0.01,
                            0.001,
                            [],
                            [fp.GaussianS(0,1,1,0.25)],
                            [fp.Region(0,1,1)],
                            0.1,
                            0.,
                            track=1,trackStats=0)
    nlist = np.array([1000]*1000,dtype=np.uint32)
    qt.evolve_qtrait_more(rng,pops,
                          nlist[0:],
                          0,
                          0.01,
                          0.001,
                          [],
                          [fp.GaussianS(0,1,1,0.25)],
                          [fp.Region(0,1,1)],
                          0.1,
                          0.,
                          track=1,trackStats=1)
    #Evolve another N gens after shift optimum to 0.5
    qt.evolve_qtrait_more(rng,pops,
                            nlist[0:],
                            0,
                            0.01,
                            0.001,
                            [],
                            [fp.GaussianS(0,1,1,0.25)],
                            [fp.Region(0,1,1)],
                            0.1,
                            optimum=0.5,
                            track=1,trackStats=1)
    st=[qt.popstats(i) for i in pops]
    for i in st:
        hdf.append('popstats',pandas.DataFrame(i))

hdf.close()
