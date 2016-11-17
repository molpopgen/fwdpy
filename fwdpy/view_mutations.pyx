#cdef popgen_mut_data get_mutation(const popgenmut & m,
cdef MutationView get_mutation(const popgenmut & m,size_t n,int index) :
    return MutationView(m.pos,n,m.h,None,m.s,m.h,m.neutral,m.xtra,index)

cdef list view_mutations_details(const mcont_t & mutations,
                                 const mcounts_cont_t & mcounts) :
    rv=[]
    for i in range(mcounts.size()):
        #skip extinct mutation
        if mcounts[i]:
            rv.append(get_mutation(mutations[i],mcounts[i],i))
    return rv

def view_mutations_singlepop(Spop p):
    return view_mutations_details(p.pop.get().mutations,p.pop.get().mcounts)

def view_mutations_singlepop_mloc(MlocusPop p):
    return view_mutations_details(p.pop.get().mutations,p.pop.get().mcounts)

def view_mutations_metapop(MetaPop p,unsigned deme):
    raise RuntimeError("needs refactoring")
    if deme >= len(p.popsizes()):
        raise IndexError("view_mutations: deme index out of range")
    #get the gametes from this population
    gams = view_gametes_metapop(p,deme)
    #extract the mutations from each gamete
    allmuts = []
    umuts = []
    for g in gams:
        for m in g['neutral']:
            allmuts.append(m)
            if umuts.count(m) == 0:
                umuts.append(m)
        for m in g['selected']:
            allmuts.append(m)
            if umuts.count(m) == 0:
                umuts.append(m)

    rv = []
    dummy=0
    for i in umuts:
        rv.append(i)
        rv[dummy]['n'] = allmuts.count(i)
        dummy+=1
    return rv

def view_mutations(object p, deme = None):
    """
    Get detailed list of all mutations in the population

    :param p: a :class:`fwdpy.fwdpy.PopType` or a :class:`fwdpy.fwdpy.PopVec`

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
    >>> #muts[0] will be a list and muts[0][0] will be a dict
    >>> muts = [fwdpy.views.view_mutations(i) for i in pops]

    .. note:: :class:`fwdpy.fwdpy.MetaPopVec` currently not supported
    """
    if isinstance(p,Spop):
        return view_mutations_singlepop(p)
    elif isinstance(p,MlocusPop):
        return view_mutations_singlepop_mloc(p)
    elif isinstance(p,SpopVec):
        return [view_mutations_singlepop(<Spop>i) for i in <SpopVec>p]# view_mutations_popvec(p)
    elif isinstance(p,MlocusPopVec):
        return [view_mutations_singlepop_mloc(<MlocusPop>i) for i in <MlocusPopVec>p]
        #return view_mutations_popvec_mloc(p)
    elif isinstance(p,MetaPop):
        if deme is None:
            raise RuntimeError("view_mutations: deme cannot be none for metapops")
        return view_mutations_metapop(p,deme)
    else:
        raise RuntimeError("view_mutations: unsupported object type")
    
