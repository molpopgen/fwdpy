import fwdpy as fp
import numpy as np
rng = fp.GSLrng(101)
N=1000
nlist=np.array([N]*10*N,dtype=np.uint32)
p=fp.evolve_regions(rng,
                    4,N,
                    nlist,
                    0.25,
                    0.0,
                    0.25,
                    [fp.Region(0,1,1)],
                    [],
                    [fp.Region(0,1,1)])
#p.append(p2)
#print len(p)
#v = fp.view_diploids_pd(p,range(0,N,1))
#v=[fp.view_diploids(i,range(0,N,1)) for i in p]
#b = fp.view_diploids_pd(p,range(0,N,1))
#print b
#for i in v:
#    print i
#g = [i.gen() for i in p]
#print g

#d = fp.view_diploids(p,[0,1])

#print d
     
