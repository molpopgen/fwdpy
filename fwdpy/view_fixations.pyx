from libcpp.vector cimport vector

cdef popgen_mut_data get_fixed_mutation(const popgenmut & m,
                                        const unsigned ftime,
                                        const unsigned N) nogil:
    cdef popgen_mut_data rv
    rv.pos=m.pos
    rv.n=2*N
    rv.g=m.g
    rv.ftime=ftime
    rv.s=m.s
    rv.h=m.h
    rv.neutral=m.neutral
    rv.label=m.xtra
    return rv

cdef vector[popgen_mut_data] view_fixations_details( const mcont_t & fixations,
                                                     const vector[uint] & fixation_times,
                                                     const unsigned N) nogil:
    cdef vector[popgen_mut_data] rv
    cdef size_t i=0,j=fixations.size()
    while i!=j:
        rv.push_back(get_fixed_mutation(fixations[i],fixation_times[i],N))
        i+=1
    return rv

def view_fixations_popvec(SpopVec p):
    cdef vector[vector[popgen_mut_data]] rv
    cdef size_t npops=p.pops.size()
    rv.resize(npops)
    cdef int i
    for i in prange(npops,schedule='static',nogil=True,chunksize=1):
        rv[i]=view_fixations_details(p.pops[i].get().fixations,p.pops[i].get().fixation_times,p.pops[i].get().N)
    return rv

def view_fixations_mlocuspopvec(MlocusPopVec p):
    cdef vector[vector[popgen_mut_data]] rv
    cdef size_t npops=p.pops.size()
    rv.resize(npops)
    cdef int i
    for i in prange(npops,schedule='static',nogil=True,chunksize=1):
        rv[i]=view_fixations_details(p.pops[i].get().fixations,p.pops[i].get().fixation_times,p.pops[i].get().N)
    return rv

def view_fixations(object p):
    """
    Return information on fixed variants

    :param p: a :class:`fwdpy.fwdpy.PopType` or a :class:`fwdpy.fwdpy.PopVec`

    :return: A list of tuples. The first element is fixation time, and the second is a dict containing data about the mutation.

    .. note:: You may need to call :func:`fwdpy.fwdpy.view_mutations` to view all types of fixations, depending on the type of simulation you are running.
    """

    #Streamline using casts:
    if isinstance(p,Spop):
        return view_fixations_details((<Spop>p).pop.get().fixations,(<Spop>p).pop.get().fixation_times,(<Spop>p).pop.get().N)
    if isinstance(p,MetaPop):
        return view_fixations_details((<MetaPop>p).mpop.get().fixations,(<MetaPop>p).mpop.get().fixation_times,sum((<MetaPop>p).mpop.get().Ns))
    elif isinstance(p,SpopVec):
        return view_fixations_popvec(p)
    elif isinstance(p,MlocusPopVec):
        return view_fixations_mlocuspopvec(p)
    else:
        raise ValueError("unsupported type")    
    
