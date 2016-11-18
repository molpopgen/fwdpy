from fwdpy.fwdpy cimport Spop,MetaPop,SpopVec,MetaPopVec,GSLrng
import math

##cdef function act directly on C++ types

cdef copy_pop_details(MetaPop mpop,size_t deme):
    copy_deme(mpop.mpop.get(),deme)

cdef remove_pop_details(MetaPop mpop,size_t deme):
    remove_deme(mpop.mpop.get(),deme)

cdef merge_pops_details(MetaPop mpop,size_t i,size_t j):
    merge_demes(mpop.mpop.get(),i,j)

cdef split_deme_details(GSLrng r,MetaPop mpop,size_t i, unsigned N_new, bint replacement):
    split_deme(r.thisptr.get(),mpop.mpop.get(),i,N_new,replacement)

cdef admix_demes_details(GSLrng r,MetaPop mpop,size_t i,size_t j,double p,unsigned N_new, bint replacement):
    admix_demes(r.thisptr.get(),mpop.mpop.get(),i,j,p,N_new,replacement)

cdef swap_demes_details(MetaPop mpop, size_t i, size_t j):
    swap_demes(mpop.mpop.get(),i,j)
    
##def function are callable from Python:
from cython.operator cimport dereference as deref

def make_MetaPopVec(SpopVec pops):
    """
    Initialize :class:`fwdpy.fwdpy.MetaPopVec` from :class:`fwdpy.fwdpy.SpopVec`

    Example:

    >>> import fwdpy as fp
    >>> import fwdpy.demography as demog
    >>> s=fp.SpopVec(40,100)
    >>> m=demog.make_MetaPopVec(s)
    """
    rv=MetaPopVec(len(pops),[0])
    for i in range(len(pops)):
        rv.mpops[i].reset(new metapop_t(deref(pops.pops[i].get())))
    return rv
    
def copy_pop(MetaPopVec mpops, size_t deme_index):
    """
    Make an exact copy of an existing deme.

    :param mpops: A :class:`fwdpy.fwdpy.MetaPopVec`
    :param deme_index: A value in the range :math:`0 \leq x \leq len(mpops)`.

    :return: nothing

    :raise: IndexError if deme_index out of range

    Example:

    >>> import fwdpy as fp
    >>> import fwdpy.demography as demog
    >>> #Init 64 metapops with 1 deme of N=1,000
    >>> m = fp.MetaPopVec(64,[1000])
    >>> len(m)
    64
    >>> m[0].popsizes()
    [1000]
    >>> demog.copy_pop(m,0)
    >>> m[0].popsizes()
    [1000, 1000]
    """
    for i in mpops:
        copy_pop_details(i,deme_index)

def merge_pops(MetaPopVec mpops, size_t deme_i, size_t deme_j):
    """
    Merge population max(deme_i,deme_j) into min(deme_i,deme_j) and remove max(deme_i,deme_j).

    :param mpops: A :class:`fwdpy.fwdpy.MetaPopVec`
    :param deme_i: A value in the range :math:`0 \leq x \leq len(mpops)`.
    :param deme_j: A value in the range :math:`0 \leq x \leq len(mpops)`.

    :return: nothing

    :raise: IndexError if deme_index out of range or RuntimeError if deme_i == deme_j

    Example:

    >>> import fwdpy as fp
    >>> import fwdpy.demography as demog
    >>> #Init 64 metapops with 3 demes of N=1,000
    >>> m = fp.MetaPopVec(64,[1000,1000,1000])
    >>> m[0].popsizes()
    [1000, 1000, 1000]
    >>> demog.merge_pops(m,1,2)
    >>> m[0].popsizes()
    [1000, 2000]
    >>> demog.merge_pops(m,0,1)
    >>> m[0].popsizes()
    [3000]
    """
    for i in mpops:
        merge_pops_details(i,deme_i,deme_j)

def remove_pop(MetaPopVec mpops, size_t deme_index):
    """
    Remove a deme from a metapopulation.

    :param mpops: A :class:`fwdpy.fwdpy.MetaPopVec`
    :param deme_index: A value in the range :math:`0 \leq x \leq len(mpops)`.

    :return: nothing

    :raise: IndexError if deme_index out of range

    Example:

    >>> import fwdpy as fp
    >>> import fwdpy.demography as demog
    >>> #Init 64 metapops with 3 demes of N=1,000
    >>> m = fp.MetaPopVec(64,[1000,1000,1000])
    >>> demog.merge_pops(m,1,2)
    >>> m[0].popsizes()
    [1000, 2000]
    >>> demog.remove_pop(m,0)
    >>> m[0].popsizes()
    [2000]
    """
    for i in mpops:
        remove_pop_details(i,deme_index)

def split_pops(GSLrng rng,MetaPopVec mpops,size_t deme_index, unsigned N_new, bint replacement = False):
    """
    Split an existing deme into two.

    :param rng: A :class:`fwdpy.fwdpy.GSLrng`
    :param mpops: A :class:`fwdpy.fwdpy.MetaPopVec`
    :param deme_index: A value in the range :math:`0 \leq x \leq len(mpops)`.
    :param N_new: The size of the new deme.
    :bint replacement: Whether to sample individuals into new deme with replacement

    :return: nothing

    :raise: IndexError if deme_index out of range or RuntimeError if N_new is larger than parental deme size.

    ..note:: The parental population's size is reduced by a value of N_new.  If splitting **without** replacement,
    then the daughter deme is generated from a list of unique diploids who are remvoed from the parental deme.  **With**
    replacement, both demes are populated by random samples from the parental deme.  The new deme is added to the end of 
    the metapopulation.

    Example:

    >>> import fwdpy as fp
    >>> import fwdpy.demography as demog
    >>> #Init 64 metapops with 1 deme of N=1,000
    >>> m = fp.MetaPopVec(64,[1000])
    >>> m[0].popsizes()
    [1000]
    >>> #Split pops into sizes 250 and 750:
    >>> rng = fp.GSLrng(100)
    >>> demog.split_pops(rng,m,0,750)
    >>> m[0].popsizes()
    [250, 750]
    """
    for i in mpops:
        split_deme_details(rng,i,deme_index,N_new,replacement)

def admix_pops(GSLrng rng, MetaPopVec mpops, size_t deme1, size_t deme2, double admix_proportion,unsigned N_new,bint replacement = False):
    """
    Admix two demes

    :param rng: A :class:`fwdpy.fwdpy.GSLrng`
    :param mpops: A :class:`fwdpy.fwdpy.MetaPopVec`
    :param deme1: A value in the range :math:`0 \leq x \leq len(mpops)`.
    :param deme2: A value in the range :math:`0 \leq x \leq len(mpops)`.
    :param admix_proportion: The probability that an individual in the admixed deme comes from deme1.
    :param N_new: population size of admixed population
    :bint replacement: Whether to sample individuals into new deme with replacement

    ..note:: The new deme is added to the end of the metapopulation.

    Example:
    
    >>> import fwdpy as fp
    >>> import fwdpy.demography as demog
    >>> #Init 64 metapops with 2 demes of N=1,000
    >>> m = fp.MetaPopVec(64,[1000,1000])
    >>> #Create new admixed deme with N=500, 25% individuals from deme 0:
    >>> rng=fp.GSLrng(100)
    >>> demog.admix_pops(rng,m,0,1,0.25,500)
    >>> m[0].popsizes()
    [1000, 1000, 500]
    """
    for i in mpops:
        admix_demes_details(rng,i,deme1,deme2,admix_proportion,N_new,replacement)

def swap_pops(MetaPopVec mpops, size_t deme_i, size_t deme_j):
    """
    Swap demes.

    :param mpops: A :class:`fwdpy.fwdpy.MetaPopVec`
    :param deme_i: A value in the range :math:`0 \leq x \leq len(mpops).`
    :param deme_j: A value in the range :math:`0 \leq x \leq len(mpops)`.

    :return: nothing

    :raise: IndexError if deme_i or deme_j is out of range

    ..note:: This function simply swaps the diploids in two demes, and is probably only useful if you
    want the deme order to be a certain way after performing other demographic operations.

    Example:

    >>> import fwdpy as fp
    >>> import fwdpy.demography as demog
    >>> #Init 64 metapops with 2 demes of different size
    >>> m = fp.MetaPopVec(64,[1000,500])
    >>> m[0].popsizes()
    [1000, 500]
    >>> demog.swap_pops(m,0,1)
    >>> m[0].popsizes()
    [500, 1000]
    """
    for i in mpops:
        swap_demes_details(i,deme_i,deme_j)

def exponential_size_change(Nstart,Nstop,time):
    """
    Generate a list of population sizes according to exponential size_change model

    :param Nstart: population size at onset of size change
    :param Nstop: Population size to reach at end of size change
    :param time: Time (in generations) to get from Nstart to Nstop

    :return: A list of integers representing population size over time.
    """
    if time < 1:
        raise RuntimeError("time must be >= 1")
    if Nstart < 1 or Nstop < 1:
        raise RuntimeError("Nstart and Nstop must both be >= 1")
    G = math.exp((math.log(Nstop) - math.log(Nstart))/time)
    rv=[]
    for i in range(time):
        rv.append(round(Nstart*pow(G,i+1)))
    return rv
