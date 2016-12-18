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

cdef copy_merge_dict(dict x,dict y):
    z=x.copy()
    z.update(y)
    return z

def tidy_trajectories(trajectories, *args):
    """
    Take a set of allele frequency trajectories and 'tidy' them for easier coercion into
    a pandas.DataFrame.

    :param trajectories: A list of mutation frequency trajectories from :class:`fwdpy.fwdpy.FreqSampler`.
    :param args: Callable objects (functions/lambdas) that return True or False.  See note below.

    .. note:: 
        
        The optional arguments must be callable functions (or lambdas) that expect a tuple
        with two elements as an argument.  The first element in the tuple is a dict
        with the following keys: pos, origin, esize, and label representing mutation
        position, generation when it first appearted, effect size, and label, respectively.
        Note that label is not currently used by many types of simualation.  The second
        element in the tuple will be a list of tuples of length two.  The first element
        is the generation and the second is the frequency in that generation.  This list
        of tuples describes the frequency trajectory of the mutation over time. The list
        is also sorted by generation. Your 
        callable should do something with these data and return True if the result is acceptable.
        Returning False means that the mutation will not be returned in the "tidied" data.

    The frequency data for each generation for each mutation are represented as a dict.

    :rtype: list of dicts

    .. versionchanged:: 0.0.4-rc2

        Added ability to filter results with custom functions.

    """
    rv=[]
    for ti in trajectories:
        if len(ti[1]) > 0:
            add_trajectory = True
            if len(args)>0:
                tests=[i(ti) for i in args]
                if tests.count(True) != len(args):
                    add_trajectory=False
            if add_trajectory is True:
                for genfreq in ti[1]:
                    rv.append(copy_merge_dict(ti[0],{'generation':genfreq[0],'freq':genfreq[1]}))
    return rv
