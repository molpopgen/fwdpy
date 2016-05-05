from fwdpy.fwdpy cimport popvec,singlepop,singlepop_t

def genotype_matrices(popvec pops):
    rv=[]
    for i in range(len(pops)):
        temp=make_geno_matrix(pops.pops[i].get(),True)
        rv.append(temp)
    return rv
