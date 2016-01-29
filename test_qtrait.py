import fwdpy as fp
import fwdpy.qtrait as qt
import numpy as np
import pandas
rng = fp.GSLrng(101)
nlist = np.array([1000]*10000,dtype=np.uint32)

hdf = pandas.HDFStore("foo2.h5",'w')
hdf.open()
k=0
for i in range(2):
    print i
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
                            track=True)
    #Evolve another 10N gens after shift optimum to 0.25
    qt.evolve_qtrait_more(rng,pops,
                            nlist[0:],
                            0,
                            0.01,
                            0.001,
                            [],
                            [fp.GaussianS(0,1,1,0.25)],
                            [fp.Region(0,1,1)],
                            0.1,
                            optimum=0.25,
                            track=True)
    #Get the mutation frequency trajectories, etc.
    #for all mutations that ever happened
    #tri is a list of pandas.DataFrame
    tri=[fp.trajectories(p,10) for p in pops]
    #Add a 'replicate ID' to each element in tri
    for j in range(len(tri)):
        tri[j]['rep']=[k]*len(tri[j].index)
        k=k+1
        #write each data frame to h5 file
        hdf.append('tr',tri[j])

hdf.close()
