#See http://docs.cython.org/src/userguide/memoryviews.html
from cython.view cimport array as cvarray
import numpy as np

def evolve_regions(GSLrng rng,
                    int npops,
                    int N,
                    unsigned[:] nlist,
                    double mu_neutral,
                    double mu_selected,
                    double recrate,
                    list nregions,
                    list sregions,
                    list recregions,
                    double f = 0,
                    const char * fitness = "multiplicative"):
    pops = popvec(npops,N)
    print len(pops)
    nreg = process_regions(nregions)
    sreg = process_regions(sregions)
    recreg = process_regions(recregions)
    print("here")
    print(nreg)
    print(nreg['beg'])
    print(nreg['beg'].tolist())
    print(nreg['end'].tolist())
    print(nreg['weight'].tolist())
    print(fitness)
    v = shwrappervec()
    process_sregion_callbacks(v,sregions)
    evolve_regions_t(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,nreg['beg'].tolist(),nreg['end'].tolist(),nreg['weight'].tolist(),
                    sreg['beg'].tolist(),sreg['end'].tolist(),sreg['weight'].tolist(),&v.vec,
                    recreg['beg'].tolist(),recreg['end'].tolist(),recreg['weight'].tolist(),
                    fitness)
    return pops
    
                
def evolve_regions_more(GSLrng rng,
                        popvec pops,
                        unsigned[:] nlist,
                        double mu_neutral,
                        double mu_selected,
                        double recrate,
                        list nregions,
                        list sregions,
                        list recregions,
                        double f = 0,
                        const char * fitness = "multiplicative"):
    nreg = process_regions(nregions)
    sreg = process_regions(sregions)
    recreg = process_regions(recregions)
    v = shwrappervec()
    process_sregion_callbacks(v,sregions)
    evolve_regions_t(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,nreg['beg'],nreg['end'],nreg['weight'],
                    sreg[0]['beg'],sreg[0]['end'],sreg[0]['weight'],&v.vec,recreg['beg'],recreg['end'],recreg['weight'],
                    fitness)
