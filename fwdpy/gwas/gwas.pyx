from fwdpy.fwdpy cimport SpopVec,Spop,singlepop_t

def genotype_matrices(SpopVec pops):
    rv=[]
    for i in range(len(pops)):
        temp=make_geno_matrix(pops.pops[i].get(),True)
        rv.append(temp)
    return rv
