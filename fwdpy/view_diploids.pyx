
cdef DiploidView get_diploid(const diploid_t & dip,
                              const gcont_t & gametes,
                              const mcont_t & mutations,
                              const mcounts_cont_t & mcounts) :
    g1=get_gamete(gametes[dip.first],mutations,mcounts)
    g2=get_gamete(gametes[dip.second],mutations,mcounts)
    return DiploidView(g1,g2,dip.g,dip.w,dip.w)

cdef MultiLocusDiploidView get_diploid_mloc (const dipvector_t & dip,
                                         const gcont_t & gametes,
                                         const mcont_t & mutations,
                                         const mcounts_cont_t & mcounts) :
    loci1=[]
    loci2=[]
    for j in range(<int>dip.size()):
        loci1.append(get_gamete(gametes[dip[j].first],mutations,mcounts))
        loci2.append(get_gamete(gametes[dip[j].second],mutations,mcounts))
    return MultiLocusDiploidView(loci1,loci2,dip[0].g,dip[0].e,dip[0].w)


cdef list view_diploids_details(const dipvector_t & diploids,
                                                const gcont_t & gametes,
                                                const mcont_t & mutations,
                                                const mcounts_cont_t & mcounts,
                                                const vector[unsigned] & indlist) :
    rv=[]
    for i in range(indlist.size()):
        if i >= diploids.size():
            raise IndexError("index greater than population size")
        rv.append(get_diploid(diploids[indlist[i]],gametes,mutations,mcounts))
    return rv


cdef list view_diploids_details_mloc(const vector[dipvector_t] & diploids,
                                                          const gcont_t & gametes,
                                                          const mcont_t & mutations,
                                                          const mcounts_cont_t & mcounts,
                                                          const vector[unsigned] & indlist) :
    rv=[]
    for i in range(indlist.size()):
        if i >= diploids.size():
            raise IndexError("index greater than population size")
        rv.append(get_diploid_mloc(diploids[indlist[i]],gametes,mutations,mcounts))
    return rv

def view_diploids_singlepop(Spop p, list indlist):
    return view_diploids_details(p.pop.get().diploids,p.pop.get().gametes,p.pop.get().mutations,p.pop.get().mcounts,indlist)

def view_diploids_singlepop_mloc(MlocusPop p, list indlist):
    return view_diploids_details_mloc(p.pop.get().diploids,p.pop.get().gametes,p.pop.get().mutations,p.pop.get().mcounts,indlist)

def view_diploids_metapop(MetaPop p, list indlist, unsigned deme):
    psizes = p.popsizes()
    for i in indlist:
        for ps in psizes:
            if i >= ps:
                raise IndexError("index greater than deme size")
    if deme >= len(p.popsizes()):
        raise IndexError("view_diploids: deme index out of range")
    return view_diploids_details(p.mpop.get().diploids[deme],p.mpop.get().gametes,p.mpop.get().mutations,p.mpop.get().mcounts,indlist)
    
def view_diploids(object p, list indlist, deme = None):
    """
    Get detailed list of a set of diploids in the population

    :param p: a :class:`fwdpy.fwdpy.PopType` or a :class:`fwdpy.fwdpy.PopVec`
    :param deme: if p is a :class`fwdpy.fwdpy.MetaPop`, deme is the index of the deme to sample
    
    :rtype: a list of dictionaries.  See Note.

    Example:

    >>> import fwdpy
    >>> import fwdpy.views
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> popsizes=np.tile(popsizes,100)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> dips = fwdpy.views.view_diploids(pops,[0,101,201,301])
    >>> dips = [fwdpy.views.view_diploids(i,[0,101,201,301]) for i in pops]

    .. note:: :class:`fwdpy.fwdpy.MetaPopVec` currently not supported
    """
    if isinstance(p,Spop):
        return view_diploids_singlepop(p,indlist)
    elif isinstance(p,MlocusPop):
        return view_diploids_singlepop_mloc(p, indlist)
    elif isinstance(p,MetaPop):
        if deme is None:
            raise RuntimeError("view_diploids: deme index required for metapopulation")
        return view_diploids_metapop(p,indlist,deme)
    elif isinstance(p,SpopVec):
        return [view_diploids_singlepop(<Spop>i,indlist) for i in <SpopVec>p]
    elif isinstance(p,MlocusPopVec):
        return [view_diploids_singlepop_mloc(<MlocusPop>i,indlist) for i in <MlocusPopVec>p]
    else:
        raise RuntimeError("view_diploids: unsupported object type")
