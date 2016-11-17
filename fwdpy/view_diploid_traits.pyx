cdef diploid_traits_singlepop(Spop p):
    rv=[]
    for i in p.pop.get().diploids:
        rv.append({'g':i.g,'e':i.e,'w':i.w})
    return rv

cdef diploid_traits_singlepop_mloc(MlocusPop p):
    rv=[]
    for i in range(p.pop.get().diploids.size()):
        rv.append({'g':p.pop.get().diploids[i][0].g,
                   'e':p.pop.get().diploids[i][0].e,
                   'w':p.pop.get().diploids[i][0].w})
    return rv

cdef diploid_traits_mpop(MetaPop m, deme):
    if deme > m.mpop.get().diploids.size():
        raise RuntimeError("deme value out of range")
    rv=[]
    for i in range(m.mpop.get().diploids[deme].size()):
        rv.append({'g':m.mpop.get().diploids[deme][i].g,
                   'e':m.mpop.get().diploids[deme][i].e,
                   'w':m.mpop.get().diploids[deme][i].w})

def diploid_traits(object p, deme = None):
    """
    Return genetic value (g), environmental value (e), and fitness (w) for all diploids.
    """
    if isinstance(p,Spop):
        return diploid_traits_singlepop(<Spop>p)
    elif isinstance(p,MlocusPop):
        return diploid_traits_singlepop_mloc(<MlocusPop>p)
    elif isinstance(p,MetaPop):
        if deme is None:
            raise RuntimeError("deme cannot be None")
        return diploid_traits_mpop(p,deme)
    if isinstance(p,SpopVec):
        return [diploid_traits_singlepop(<Spop>i) for i in <SpopVec>p]
    if isinstance(p,MlocusPopVec):
        return [diploid_traits_singlepop_mloc(<MlocusPop>i) for i in <MlocusPopVec>i]
    elif isinstance(p,MetaPopVec):
        if deme is None:
            raise RuntimeError("deme cannot be None")
        return [diploid_traits_mpop(<MetaPop>i,deme) for i in <MetaPopVec>p]
    else:
        raise ValueError("unsupported type")    
