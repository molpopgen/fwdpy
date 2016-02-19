import fwdpy as fp
import fwdpy.qtrait as qt
import numpy as np
import pandas,cPickle as pickle,gzip

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
NCORES=4
NPIK=NCORES*NB
f=gzip.open(PIK,"wb")
N=1000
scaled_sigmaMU=250.0
sigMU=scaled_sigmaMU/N
theta_n = 100.0
rho_n = 100.0
mun=theta_n/(4*N)
littler=rho_n/(4*N)
rest=r-littler
ratio=rest/r
sample_interval=0.01 #In units of N generations
print sigMU," ",mun," ",int(sample_interval*N)
neutmutregions=[fp.Region(0,1,1)]
selmutregions=[fp.GaussianS(-ratio,1+ratio,1,sigMU)]
recregions= [fp.Region(-ratio,1+ratio,1)]

pickle.dump(NPIK,f)
    
for i in range(NB):
    nlist = np.array([N]*(10*N),dtype=np.uint32)
    #Evolve to equilibrium
    pops = qt.evolve_qtrait(rng,
                            NCORES,
                            N,
                            nlist[0:],
                            mun,
                            mu,
                            r,
                            neutmutregions,
                            selmutregions,
                            recregions,
                            sigE,
                            0.)
    #Evolve another 3*N gens after shift optimum to 0.5
    nlist = np.array([N]*(3*N),dtype=np.uint32)
    samples = qt.evolve_qtrait_sample(rng,pops,
                                      nlist[0:],
                                      mun,
                                      mu,
                                      r,
                                      neutmutregions,
                                      selmutregions,
                                      recregions,
                                      sigE,
                                      trackSamples=int(sample_interval*N),nsam=20,
                                      optimum=0.5)
    #with open(PIK,"ab") as f:
    for j in samples:
        pickle.dump(j,f)
#hdf.close()
