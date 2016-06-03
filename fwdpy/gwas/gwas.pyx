from fwdpy.fwdpy cimport popvec,singlepop,singlepop_t
from fwdpy.samplers cimport additive_variance_wrapper

def genotype_matrices(popvec pops):
    rv=[]
    for i in range(len(pops)):
        temp=make_geno_matrix(pops.pops[i].get(),True)
        rv.append(temp)
    return rv

def VAVG(object p):
    if isinstance(p,singlepop):
        #Return element 0, as there is no need to return a vector, as not running sampler over time...
        return additive_variance_wrapper[singlepop_t]((<singlepop>p).pop.get())[0]
    elif isinstance(p,popvec):
        return [additive_variance_wrapper[singlepop_t]((<singlepop>i).pop.get())[0] for i in popvec]
