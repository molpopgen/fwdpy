from fwdpy.fwdpy cimport SpopVec,Spop,singlepop_t
from fwdpy.fwdpy cimport additive_variance_wrapper

def genotype_matrices(SpopVec pops):
    rv=[]
    for i in range(len(pops)):
        temp=make_geno_matrix(pops.pops[i].get(),True)
        rv.append(temp)
    return rv

def VAVG(object p):
    if isinstance(p,Spop):
        return additive_variance_wrapper[singlepop_t]((<Spop>p).pop.get())
    elif isinstance(p,SpopVec):
        return [additive_variance_wrapper[singlepop_t]((<Spop>i).pop.get()) for i in (<SpopVec>(p))]
