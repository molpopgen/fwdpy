cdef GameteView get_gamete(const gamete_t & g,
                            const mcont_t & mutations,
                            const mcounts_cont_t & mcounts,int index) :
    n=[]
    for i in range(<int>g.mutations.size()):
        n.append(get_mutation(mutations[g.mutations[<size_t>i]],mcounts[g.mutations[<size_t>i]],g.mutations[<size_t>i]))
    s=[]
    for i in range(<int>g.smutations.size()):
        s.append(get_mutation(mutations[g.smutations[<size_t>i]],mcounts[g.smutations[<size_t>i]],g.smutations[<size_t>i]))
    return GameteView(n,s,g.n,index)

cdef list view_gametes_details(const gcont_t & gametes,const mcont_t & mutations, const vector[uint] & mcounts):
    rv=[]
    for i in range(gametes.size()):
        if gametes[i].n:
            rv.append(get_gamete(gametes[i],mutations,mcounts,i))
    return rv

def view_gametes_metapop(MetaPop p, unsigned deme):
    raise RuntimeError("needs refactoring")
    if deme >= len(p.popsizes()):
        raise IndexError("view_gametes: deme index out of range")
    temp1 = view_diploids(p,list(range(p.mpop.get().diploids[deme].size())),deme)
    #Get unique list of haplotypes
    unique_gams = []
    allgams = []
    for i in temp1:
        if unique_gams.count(i['chrom0'])==0:
            unique_gams.append(i['chrom0'])
        if unique_gams.count(i['chrom1'])==0:
            unique_gams.append(i['chrom1'])
        allgams.append(i['chrom0'])
        allgams.append(i['chrom1'])
    #clear temp1 and fill it with unique gametes + their counts in this deme
    temp1=[]
    dummy=0
    for i in unique_gams:
        temp1.append(i)
        temp1[dummy]['n'] = allgams.count(i)
        dummy+=1
    return temp1

def view_gametes(object p ,deme = None):
    """
    Get detailed list of all gametes in the population

    :param p: a :class:`fwdpy.fwdpy.PopType` or a :class:`fwdpy.fwdpy.PopVec`
    :param deme: If p is a :class:`fwdpy.fwdpy.MetaPop`, deme is the index of the deme to view

    :rtype: a list of dictionaries.  See note.

    Example for a single deme:

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
    >>> dips = [fwdpy.views.view_gametes(i) for i in pops]

    .. note:: :class:`fwdpy.fwdpy.MetaPopVec` currently not supported
    """
    if isinstance(p,Spop):
        return view_gametes_details((<Spop>p).pop.get().gametes,(<Spop>p).pop.get().mutations,(<Spop>p).pop.get().mcounts)
    elif isinstance(p,MlocusPop):
        return view_gametes_details((<MlocusPop>p).pop.get().gametes,(<MlocusPop>p).pop.get().mutations,(<MlocusPop>p).pop.get().mcounts)
    elif isinstance(p,SpopVec):
        return [view_gametes_details((<Spop>i).pop.get().gametes,(<Spop>i).pop.get().mutations,(<Spop>i).pop.get().mcounts) for i in <SpopVec>p]
    elif isinstance(p,MlocusPopVec):
        return [view_gametes_details((<MlocusPop>i).pop.get().gametes,(<MlocusPop>i).pop.get().mutations,(<MlocusPop>i).pop.get().mcounts) for i in <MlocusPopVec>p]
    elif isinstance(p,MetaPop):
        if deme is None:
            raise RuntimeError("view_gametes: deme cannot be None when p is a MetaPop")
        return view_gametes_metapop(p,deme)
    else:
        raise RuntimeError("view_gametes: unsupported object type")
