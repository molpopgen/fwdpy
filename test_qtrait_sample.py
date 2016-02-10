import fwdpy as fp
import fwdpy.qtrait as qt
import numpy as np
import pandas,pickle,gzip

rng = fp.GSLrng(101)

mu=0.000625
r=0.5
sigE=0.1
#hdf = pandas.HDFStore("qtrait_test_out.h5",'w',complevel=6,complib="zlib")
#hdf.open()
k=0
#2,000 replicates in 4-thread chunks
PIK="qtrait_pickle.dat"
NB=500
NPIK=4*NB
f=gzip.open(PIK,"wb")
mun=0.01
littler=0.01
rest=r-littler
ratio=rest/r

neutmutregions=[fp.Region(0,1,1)]
selmutregions=[fp.GaussianS(-ratio,1+ratio,1,0.25)]
recregions= [fp.Region(-ratio,1+ratio,1)]

pickle.dump(NPIK,f)
    
for i in range(NB):
    nlist = np.array([1000]*10000,dtype=np.uint32)
    #Evolve to equilibrium
    pops = qt.evolve_qtrait(rng,
                            4,
                            1000,
                            nlist[0:],
                            mun,
                            mu,
                            r,
                            neutmutregions,
                            selmutregions,
                            recregions,
                            sigE,
                            0.,
                            track=0,trackStats=0)
    #Evolve another N gens after shift optimum to 0.5
    nlist = np.array([1000]*3000,dtype=np.uint32)
    samples = qt.evolve_qtrait_sample(rng,pops,
                                      nlist[0:],
                                      mun,
                                      mu,
                                      r,
                                      neutmutregions,
                                      selmutregions,
                                      recregions,
                                      sigE,
                                      trackSamples=10,nsam=20,
                                      optimum=0.5)
    #with open(PIK,"ab") as f:
    for j in samples:
        pickle.dump(j,f)
#hdf.close()
