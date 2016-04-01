from fwdpy cimport singlepop,metapop,popvec,mpopvec

##cdef function act directly on C++ types

cdef copy_pop_details(metapop mpop,int deme):
    copy_deme(mpop.mpop.get(),deme)

cdef remove_pop_details(metapop mpop,int deme):
    remove_deme(mpop.mpop.get(),deme)

cdef merge_pops_details(metapop mpop,int i,int j):
    merge_demes(mpop.mpop.get(),i,j)



##def function are callable from Python:

def make_mpopvec(popvec pops):
    """
    Initialize :class:`fwdpy.mpopvec` from :class:`fwdpy.popvec`

    Example:

    >>> import fwdpy as fp
    >>> import fwdpy.demography as demog
    >>> s=fp.popvec(40,100)
    >>> m=demog.make_mpopvec(s)
    """
    rv=mpopvec(len(pops),[0])
    for i in range(len(pops)):
        rv[i].from_singlepop(pops[i])
    return rv;
    
def copy_pop(mpopvec mpops, int deme):
    for i in mpops:
        copy_pop_details(i,deme)

def merge_pops(mpopvec mpops, int deme_i, int deme_j):
    for i in mpops:
        merge_pops_details(i,deme_i,deme_j)

def remove_pop(mpopvec mpops, int deme):
    for i in mpops:
        remove_pop_details(i,deme)
