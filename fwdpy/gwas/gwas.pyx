from fwdpy.fwdpy cimport popvec,singlepop,singlepop_t
from fwdpy.fwdpy cimport additive_variance_wrapper

def genotype_matrices(popvec pops):
    rv=[]
    for i in range(len(pops)):
        temp=make_geno_matrix(pops.pops[i].get(),True)
        rv.append(temp)
    return rv

def VAVG(object p):
    if isinstance(p,singlepop):
        return additive_variance_wrapper[singlepop_t]((<singlepop>p).pop.get())
    elif isinstance(p,popvec):
        return [additive_variance_wrapper[singlepop_t]((<singlepop>i).pop.get()) for i in (<popvec>(p))]
