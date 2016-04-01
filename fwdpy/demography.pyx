from fwdpp cimport *
from fwdpy cimport singlepop,metapop,popvec,mpopvec

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
    
