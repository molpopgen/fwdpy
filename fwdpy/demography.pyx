from fwdpy cimport singlepop,metapop,popvec,mpopvec,GSLrng

##cdef function act directly on C++ types

cdef copy_pop_details(metapop mpop,size_t deme):
    copy_deme(mpop.mpop.get(),deme)

cdef remove_pop_details(metapop mpop,size_t deme):
    remove_deme(mpop.mpop.get(),deme)

cdef merge_pops_details(metapop mpop,size_t i,size_t j):
    merge_demes(mpop.mpop.get(),i,j)

cdef split_deme_details(GSLrng r,metapop mpop,size_t i, unsigned N_new, bint replacement):
    split_deme(r.thisptr.get(),mpop.mpop.get(),i,N_new,replacement)

cdef admix_demes_details(GSLrng r,metapop mpop,size_t i,size_t j,double p,bint replacement):
    admix_demes(r.thisptr.get(),mpop.mpop.get(),i,j,p,replacement)

cdef swap_demes_details(metapop mpop, size_t i, size_t j):
    swap_demes(mpop.mpop.get(),i,j)
    
##def function are callable from Python:

def make_mpopvec(popvec pops):
    """
    Initialize :class:`fwdpy.fwdpy.mpopvec` from :class:`fwdpy.fwdpy.popvec`

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
    
def copy_pop(mpopvec mpops, size_t deme_index):
    """
    Make an exact copy of an existing deme.

    :param mpops: A :class:`fwdpy.fwdpy.mpopvec`
    :param deme_index: A value in the range :math:`0 \leq x \leq len(mpops)`.

    :return: nothing

    :raise: IndexError if deme_index out of range
    """
    for i in mpops:
        copy_pop_details(i,deme_index)

def merge_pops(mpopvec mpops, size_t deme_i, size_t deme_j):
    """
    Merge population max(deme_i,deme_j) into min(deme_i,deme_j) and remove max(deme_i,deme_j).

    :param mpops: A :class:`fwdpy.fwdpy.mpopvec`
    :param deme_i: A value in the range :math:`0 \leq x \leq len(mpops)`.
    :param deme_j: A value in the range :math:`0 \leq x \leq len(mpops)`.

    :return: nothing

    :raise: IndexError if deme_index out of range or RuntimeError if deme_i == deme_j
    """
    for i in mpops:
        merge_pops_details(i,deme_i,deme_j)

def remove_pop(mpopvec mpops, size_t deme_index):
    """
    Remove a deme from a metapopulation.

    :param mpops: A :class:`fwdpy.fwdpy.mpopvec`
    :param deme_index: A value in the range :math:`0 \leq x \leq len(mpops)`.

    :return: nothing

    :raise: IndexError if deme_index out of range
    """
    for i in mpops:
        remove_pop_details(i,deme_index)

def split_pops(GSLrng rng,mpopvec mpops,size_t deme_index, unsigned N_new, bint replacement):
    """
    Split an existing deme into two.

    :param rng: A :class:`fwdpy.fwdpy.GSLrng`
    :param mpops: A :class:`fwdpy.fwdpy.mpopvec`
    :param deme_index: A value in the range :math:`0 \leq x \leq len(mpops)`.
    :param N_new: The size of the new deme.
    :bint replacement: Whether to sample individuals into new deme with replacement

    :return: nothing

    :raise: IndexError if deme_index out of range or RuntimeError if N_new is larger than parental deme size.

    ..note:: The parental population's size is reduced by a value of N_new.  If splitting **without** replacement,
    then the daughter deme is generated from a list of unique diploids who are remvoed from the parental deme.  **With**
    replacement, both demes are populated by random samples from the parental deme.  The new deme is added to the end of 
    the metapopulation.
    """
    for i in mpops:
        split_deme_details(rng,i,deme_index,N_new,replacement)

def admix_pops(GSLrng rng, mpopvec mpops, size_t deme1, size_t deme2, double admix_proportion,bint replacement):
    """
    Admix two demes

    :param rng: A :class:`fwdpy.fwdpy.GSLrng`
    :param mpops: A :class:`fwdpy.fwdpy.mpopvec`
    :param deme1: A value in the range :math:`0 \leq x \leq len(mpops)`.
    :param deme2: A value in the range :math:`0 \leq x \leq len(mpops)`.
    :param admix_proportion: The probability that an individual in the admixed deme comes from deme1.
    :bint replacement: Whether to sample individuals into new deme with replacement

    ..note:: The new deme is added to the end of the metapopulation.
    """
    for i in mpops:
        admix_demes_details(rng,i,deme1,deme2,admix_proportion,replacement)

def swap_pops(mpopvec mpops, size_t deme_i, size_t deme_j):
    """
    Swap demes.

    :param mpops: A :class:`fwdpy.fwdpy.mpopvec`
    :param deme_i: A value in the range :math:`0 \leq x \leq len(mpops).`
    :param deme_j: A value in the range :math:`0 \leq x \leq len(mpops)`.

    :return: nothing

    :raise: IndexError if deme_i or deme_j is out of range

    ..note:: This function simply swaps the diploids in two demes, and is probably only useful if you
    want the deme order to be a certain way after performing other demographic operations.
    """
    for i in mpops:
        swap_demes_details(i,deme_i,deme_j)
