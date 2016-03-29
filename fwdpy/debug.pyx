def check_popdata_singlepop(singlepop p):
    cdef bint csum = check_sum[gcont_t](p.pop.get().gametes,2*(<unsigned>p.pop.get().diploids.size()))
    cdef bint pds = popdata_sane[dipvector_t,gcont_t,mcont_t](p.pop.get().diploids,p.pop.get().gametes,p.pop.get().mutations,p.pop.get().mcounts)
    return {'check_sum':csum,'popdata_sane':pds}

def check_popdata_metapop(metapop p):
    cdef bint csum = check_sum[gcont_t](p.mpop.get().gametes,2*(<unsigned>p.mpop.get().diploids.size()))
    lpds=[]
    cdef bint pds
    cdef size_t i = 0
    cdef size_t ndemes = p.mpop.get().diploids.size()
    while i < ndemes:
        pds=popdata_sane[dipvector_t,gcont_t,mcont_t](p.mpop.get().diploids[i],p.mpop.get().gametes,p.mpop.get().mutations,p.mpop.get().mcounts)
        lpds.append(pds)
        i+=1
    return {'check_sum':csum,'popdata_sane':lpds}

def check_popdata_popvec(popvec p):
    return [check_popdata_singlepop(i) for i in p]

def check_popdata_mpopvec(mpopvec p):
    return [check_popdata_metapop(i) for i in p]

def check_popdata(object p):
    """
    Apply fwdpp's debugging functions to population containers.
    
    :param p: A object of type :class:`fwdpy.fwdpy.poptype` or :class:`fwdpy.fwdpy.poptype`

    :rtype: Dictionary with return values (True or False). Any false values reflect a critical data inconsistency, and mean an exception should be raised.
    """
    if isinstance(p,popvec):
        return check_popdata_popvec(p)
    elif isinstance(p,mpopvec):
        return check_popdata_mpopvec(p)
    elif isinstance(p,singlepop):
        return check_popdata_singlepop(p)
    elif isinstance(p,metapop):
        return check_popdata-metapop(p)
    else:
        raise RuntimeError("object type not understood")
