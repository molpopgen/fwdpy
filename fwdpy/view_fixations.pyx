from libcpp.vector cimport vector

cdef MutationView get_fixed_mutation(const popgenmut & m,
                                        const unsigned ftime,
                                        const unsigned N):
    return MutationView(m.pos,2*N,m.g,ftime,m.s,m.h,m.neutral,m.xtra)

cdef list view_fixations_details( const mcont_t & fixations,
                                 const vector[uint] & fixation_times,
                                 const unsigned N):
    rv=[]
    if fixation_times.empty():
        rv.append(empty_MutationView())
    else:
        for i in range(fixation_times.size()):
            rv.append(get_fixed_mutation(fixations[i],fixation_times[i],N))
    return rv

def view_fixations(object p):
    """
    Return information on fixed variants

    :param p: a :class:`fwdpy.fwdpy.PopType` or a :class:`fwdpy.fwdpy.PopVec`

    :return: A list of tuples. The first element is fixation time, and the second is a dict containing data about the mutation.

    .. note:: You may need to call :func:`fwdpy.views.view_mutations` to view all types of fixations, depending on the type of simulation you are running.
    """
    if isinstance(p,Spop):
        return view_fixations_details((<Spop>p).pop.get().fixations,(<Spop>p).pop.get().fixation_times,(<Spop>p).pop.get().N)
    if isinstance(p,MetaPop):
        return view_fixations_details((<MetaPop>p).mpop.get().fixations,(<MetaPop>p).mpop.get().fixation_times,sum((<MetaPop>p).mpop.get().Ns))
    elif isinstance(p,SpopVec):
        return [view_fixations_details((<Spop>i).pop.get().fixations,(<Spop>i).pop.get().fixation_times,(<Spop>i).pop.get().N) for i in p]
    elif isinstance(p,MlocusPopVec):
        return [view_fixations_details((<MlocusPop>i).pop.get().fixations,(<MlocusPop>i).pop.get().fixation_times,(<MlocusPop>i).pop.get().N) for i in p]
    else:
        raise ValueError("unsupported type")    
    
