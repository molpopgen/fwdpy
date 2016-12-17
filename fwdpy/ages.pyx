from cython.parallel import parallel, prange
from cython.operator cimport dereference as deref
from libcpp.limits cimport numeric_limits
import multiprocessing as mp

def allele_ages( list trajectories, double minfreq = 0.0, unsigned minsojourn = 1 ):
    """
    Calculate allele age information from mutation frequency trajectories.

    The return value is a list of dicts that include the generation when the mutation arose, its effect size, 
    maximum frequency, and the number of times a frequency was recorded for that mutation.
    """
    return  allele_ages_details(trajectories,minfreq,minsojourn)

def merge_trajectories(list trajectories1, list trajectories2):
                       
    """
    Take two sets of mutation trajectories and merge them.

    The intended use case is that trajectories2 are from a later time point of the same
    simulations used to generate trajectories 1.  For example, trajectories 1 may represent 
    mutations arising during evolution to equilibrium, and trajectories2 represent what happend
    during a bottleneck.
    """
    if len(trajectories1) != len(trajectories2):
        raise RuntimeError("the two input lists must be the same length")
    return merge_trajectories_details(trajectories1,trajectories2)

cdef merge_and_copy_dict(dict x,dict y):
    z=x.copy()
    z.update(y)
    return z

cdef tidy_trajectories_details(freqTraj & trajectories,unsigned min_sojourn, double min_freq,
        unsigned remove_arose_after,unsigned remove_gone_before):
    rv=[]
    for ti in trajectories:
        if ti.second.empty() is False:
            if ti.first.origin <= remove_arose_after and ti.second.back().first > remove_gone_before: 
                if ti.second.size()>=min_sojourn or ti.second.back().second==1.0:
                    passes_minfreq = False
                    for x in ti.second:
                        if x.second >= min_freq:
                            passes_minfreq = True
                            break
                    if passes_minfreq is True:
                        rv.extend([merge_and_copy_dict(ti.first,{'generation':i.first,'freq':i.second}) for i in ti.second])
    return rv

#def tidy_trajectories_work(ti):
#    assert(isinstance(ti,tuple))
#    return [merge_and_copy_dict(ti[0],{'generation':i[0],'freq':i[1]}) for i in ti[1]]
#
#cdef tidy_trajectories_details(freqTraj & trajectories,unsigned min_sojourn, double min_freq,
#        unsigned remove_arose_after,unsigned remove_gone_before):
#    temp=[]
#    for ti in trajectories:
#        if ti.second.empty() is False:
#            if ti.first.origin <= remove_arose_after and ti.second.back().first > remove_gone_before: 
#                if ti.second.size()>=min_sojourn or ti.second.back().second==1.0:
#                    passes_minfreq = False
#                    for x in ti.second:
#                        if x.second >= min_freq:
#                            passes_minfreq = True
#                            break
#                    if passes_minfreq is True:
#                        temp.append(ti)
#    rv=[]
#    P=mp.Pool(mp.cpu_count())
#    for i in P.imap_unordered(tidy_trajectories_work,temp,chunksize=100):
#        rv.extend(i)
#    temp=[]
#    P.close()
#    del P
#    del temp
#    return rv

def tidy_trajectories(trajectories, unsigned min_sojourn = 0, double min_freq = 0.0,
        remove_arose_after = None,remove_gone_before = None):
    """
    Take a set of allele frequency trajectories and 'tidy' them for easier coercion into
    a pandas.DataFrame.

    :param trajectories: A list of mutation frequency trajectories from :class:`fwdpy.fwdpy.FreqSampler`.
    :param min_sojourn: Exclude mutations that segregate for fewer generations than this value.
    :param min_freq: Exclude mutations that never reach a frequency :math:`\\geq` this value.
    :param remove_arose_after: (None) Do not include mutations whose origin times are > this value.
    :param remove_gone_before: (None) Do not include mutations whose last recorded generation are <= this value.

    .. note:: The sojourn time filter is not applied to fixations.  I'm assuming you are always interested in those.

    The frequency data for each generation for each mutation are represented as a dict.

    :rtype: list of dicts or generator to such a list
    """
    cdef numeric_limits[unsigned] ul
    cdef unsigned raf = ul.max()
    cdef unsigned rgb = 0
    if remove_arose_after is not None:
        raf=remove_arose_after
    if remove_gone_before is not None:
        rgb=remove_gone_before
    return tidy_trajectories_details(trajectories,min_sojourn,min_freq,raf,rgb)

