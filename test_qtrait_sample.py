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
NB=50
NCORES=40
NPIK=NCORES*NB
f=gzip.open(PIK,"wb")
N=1000
scaled_sigmaMU=250.0
sigMU=scaled_sigmaMU/N
theta_n = 100.0
rho_n = 100.0
mun=0.0#theta_n/(4*N)
littler=rho_n/(4*N)
rest=r-littler
ratio=rest/r
sample_interval=0.01 #In units of N generations
print sigMU," ",mun," ",int(sample_interval*N)
neutmutregions=[fp.Region(0,0.1,1)]
#selmutregions=[fp.GaussianS(-ratio,1+ratio,1,sigMU)]
#recregions= [fp.Region(-ratio,1+ratio,1)]
selmutregions=[fp.GaussianS(0,0.1,1,sigMU)]
recregions= [fp.Region(0,0.1,1)]

pickle.dump(NPIK,f)
REP=0
for i in range(NB):
    nlist = np.array([N]*(10*N),dtype=np.uint32)
    #Evolve to equilibrium
    pops=fp.popvec(NCORES,N)
    samples = qt.evolve_qtrait_track(rng,
                                     pops,
                                     nlist[0:],
                                     mun,
                                     mu,
                                     r,
                                     neutmutregions,
                                     selmutregions,
                                     recregions,
                                     sigE,
                                     track=1,
                                     optimum=0.)
    #Evolve another 3*N gens after shift optimum to 0.5
    nlist = np.array([N]*(3*N),dtype=np.uint32)
    samples2 = qt.evolve_qtrait_track(rng,pops,
                                      nlist[0:],
                                      mun,
                                      mu,
                                      r,
                                      neutmutregions,
                                      selmutregions,
                                      recregions,
                                      sigE,
                                      track=1,
                                      optimum=1.5)
    for j in range(len(samples)):
        df = pandas.concat([pandas.DataFrame(samples[j]),pandas.DataFrame(samples2[j])])
        for name,group in df.groupby(['pos','esize']):
            if group.freq.max() < 1:
                if group.generation.max()-group.generation.min()>1:
                    print REP," ",name[0]," ",name[1]," ",group.generation.min()," ",group.generation.max(),' ',group.freq.max(),' ',len(group.freq),' ',group.freq.iloc[-1]
        REP+=1
#hdf.close()
