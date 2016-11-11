from fwdpy.fwdpy cimport SpopVec,Spop,singlepop_t

def genotype_matrices(SpopVec pops,list individuals = None):
    rv=[]
    for i in range(len(pops)):
        if individuals is None:
            temp=make_geno_matrix(pops.pops[i].get(),True)
            rv.append(temp)
        else:
            temp=make_geno_matrix(pops.pops[i].get(),individuals,True)
            rv.append(temp)
    return rv
