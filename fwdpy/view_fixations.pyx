from libcpp.vector cimport vector
from libcpp.utility cimport pair

cdef vector[pair[uint,popgen_mut_data]] view_fixations_details( const mcont_t & fixations,
                                                                const vector[uint] & fixation_times,
                                                                const unsigned N) nogil:
    cdef vector[pair[uint,popgen_mut_data]] rv
    cdef size_t i=0,j=fixations.size()
    while i!=j:
        rv.push_back(pair[uint,popgen_mut_data](fixation_times[i],get_mutation(fixations[i],N)))
        i+=1
    return rv

def view_fixations_popvec(popvec p):
    cdef vector[vector[pair[uint,popgen_mut_data]]] rv
    cdef size_t npops=p.pops.size()
    cdef int i
    for i in prange(npops,schedule='static',nogil=True,chunksize=1):
        rv[i]=view_fixations_details(p.pops[i].get().fixations,p.pops[i].get().fixation_times,p.pops[i].get().N)
    return rv

def view_fixations_mpopvec(mpopvec p):
    cdef vector[vector[pair[uint,popgen_mut_data]]] rv
    cdef size_t npops=p.mpops.size()
    cdef int i
    ##HACK ALERT
    cdef vector[unsigned] Ns
    for i in range(npops):
        Ns.push_back(sum(p.mpops[i].get().Ns))
    for i in prange(npops,schedule='static',nogil=True,chunksize=1):
        rv[i]=view_fixations_details(p.mpops[i].get().fixations,p.mpops[i].get().fixation_times,Ns[i])
    return rv

def view_fixations(object p):
    """
    Return information on fixed variants

    :param p: a :class:`fwdpy.fwdpy.poptype` or a :class:`fwdpy.fwdpy.popvec`

    :return: A list of tuples. The first element is fixation time, and the second is a dict containing data about the mutation.

    .. note:: You may need to call :func:`fwdpy.fwdpy.view_mutations` to view all types of fixations, depending on the type of simulation you are running.
    """

    #Streamline using casts:
    if isinstance(p,Spop):
        return view_fixations_details((<Spop>p).pop.get().fixations,(<Spop>p).pop.get().fixation_times,(<Spop>p).pop.get().N)
    if isinstance(p,MetaPop):
        return view_fixations_details((<MetaPop>p).mpop.get().fixations,(<MetaPop>p).mpop.get().fixation_times,sum((<MetaPop>p).mpop.get().Ns))
    
    elif isinstance(p,popvec):
        return view_fixations_popvec(p)
    elif isinstance(p,mpopvec):
        return view_fixations_mpopvec(p)
    else:
        raise ValueError("unsupported type")    
