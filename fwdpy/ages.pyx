from cython.parallel import parallel, prange
from cython.operator cimport dereference as deref
from libcpp.limits cimport numeric_limits
from cython.operator cimport dereference as deref

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

cdef mergedict(dict x,dict y):
    x.update(y)
    return x

def tidy_trajectories_details(freqTraj & trajectories,tuple args):
    rv=[]
    cdef vector[pair[selected_mut_data,vector[genfreqPair]]].iterator ti = trajectories.begin()
    tie = trajectories.end()
    cdef vector[pair[uint,double]].iterator sb,se
    while ti != tie:
        if deref(ti).second.empty() is False:
            add_trajectory = True
            if len(args) > 0:
                tests = [ai(deref(ti)) for ai in args]
                if tests.count(True) != len(args):
                    add_trajectory = False
            if add_trajectory:
                sb = deref(ti).second.begin()
                se = deref(ti).second.end()
                while sb<se:
                    rv.append(mergedict(deref(ti).first,{'generation':deref(sb).first,'freq':deref(sb).second}))
                    sb+=1
        ti+=1
    return rv

def tidy_trajectories(trajectories, *args):
    """
    Take a set of allele frequency trajectories and 'tidy' them for easier coercion into
    a pandas.DataFrame.

    :param trajectories: A list of mutation frequency trajectories from :class:`fwdpy.fwdpy.FreqSampler`.
    :param args: Callable objects (functions/lambdas) that return True or False.  See note below.

    .. note:: The sojourn time filter is not applied to fixations.  I'm assuming you are always interested in those.

    The frequency data for each generation for each mutation are represented as a dict.

    :rtype: list of dicts or generator to such a list
    """
    return tidy_trajectories_details(trajectories,*args)

