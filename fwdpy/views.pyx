from cython.operator import dereference as deref,postincrement as inc

cdef get_mutation( const cpplist[popgenmut].iterator & itr):
    return {'pos':deref(itr).pos,'n':deref(itr).n,'g':deref(itr).g,'s':deref(itr).s,'h':deref(itr).h}

cdef get_gamete( const cpplist[gamete_t].iterator & itr ):
    cdef vector[cpplist[popgenmut].iterator].iterator beg = deref(itr).mutations.begin()
    cdef vector[cpplist[popgenmut].iterator].iterator end = deref(itr).mutations.end()
    neutral = []
    selected = []
    
    while beg != end:
        neutral.append(get_mutation(deref(beg)))
        inc(beg)

    beg = deref(itr).smutations.begin()
    end = deref(itr).smutations.end()
    while beg != end:
        selected.append(get_mutation(deref(beg)))
        inc(beg)
    return {'n':deref(itr).n,'neutral':neutral,'selected':selected}

cdef get_diploid( const vector[diploid_t].iterator itr ):
    return {'chrom0':get_gamete(deref(itr).first),
            'chrom1':get_gamete(deref(itr).second)}

cdef view_mutations_details(cpplist[popgenmut].iterator beg,cpplist[popgenmut].iterator end):
    rv=[]
    while beg != end:
        rv.append(get_mutation(beg))
        inc(beg)
    return rv

cdef view_gametes_details( cpplist[gamete_t].iterator beg,
                           cpplist[gamete_t].iterator end ):
    rv=[]
    while beg != end:
        rv.append(get_gamete(beg))
        inc(beg)
    return rv

##This really should be const...
cdef view_diploids_details( vector[diploid_t] & diploids,
                            const vector[unsigned] indlist ):
    cdef vector[diploid_t].iterator itr = diploids.begin()
    rv=[]
    for i in range(indlist.size()):
        if i >= diploids.size():
            raise IndexError("view_diploids: index out of range")
        rv.append(get_diploid(itr+indlist[i]))
    return rv

def view_mutations_singlepop(singlepop p):
    cdef cpplist[popgenmut].iterator beg = p.pop.get().mutations.begin()
    cdef cpplist[popgenmut].iterator end = p.pop.get().mutations.end()
    return view_mutations_details(beg,end)

def view_mutations_metapop(metapop p):
    cdef cpplist[popgenmut].iterator beg = p.mpop.get().mutations.begin()
    cdef cpplist[popgenmut].iterator end = p.mpop.get().mutations.end()
    return view_mutations_details(beg,end)

def view_mutations( poptype p ):
    """
    Get detailed list of all mutations in the population

    :param p: a :class:`fwdpy.fwdpy.poptype`

    :rtype: a list of dictionaries.  See Note.

    Example:

    >>> import fwdpy
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> popsizes=np.tile(popsizes,10000)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> muts = [fwdpy.view_mutations(i) for i in pops]
    >>> type(muts[0])
    <type 'list'>
    >>> type(muts[0][0])
    <type 'dict'>
    """
    if isinstance(p,singlepop):
        return view_mutations_singlepop(p)
    elif isinstance(p,metapop):
        return view_mutations_metapop(p)
    else:
        raise RuntimeError("view_mutations: unsupported poptype")
    
def view_gametes_singlepop( singlepop p ):
    cdef cpplist[gamete_t].iterator beg = p.pop.get().gametes.begin()
    cdef cpplist[gamete_t].iterator end = p.pop.get().gametes.end()
    return view_gametes_details(beg,end)

def view_gametes_metapop( metapop p, deme):
    if deme >= len(p.popsizes()):
        raise IndexError("view_gametes: deme index out of ramge")
    temp1 = view_diploids(p,range(p.mpop.get().diploids[deme].size()),deme)
    #Get unique list of haplotypes
    temp2 = []
    temp3 = []
    for i in temp1:
        if temp2.count(i['chrom0'])==0:
            temp2.append(i['chrom0'])
        if temp2.count(i['chrom1'])==0:
            temp2.append(i['chrom1'])
        temp3.append(i['chrom0'])
        temp3.append(i['chrom1'])
    #clear temp1 and fill it with unique gametes + their counts in this deme
    temp1=[]
    dummy=0
    for i in temp2:
        temp1.append(i)
        temp1[dummy]['n'] = temp3.count(i)
        dummy+=1
    return temp1

def view_gametes( poptype p ,deme = None):
    """
    Get detailed list of all gametes in the population

    :param p: a :class:`fwdpy.fwdpy.poptype`
    :param deme: If p is a :class:`fwdpy.fwdpy.metapop`, deme is the index of the deme to view

    :rtype: a list of dictionaries.  See note.

    Example for a single deme:

    >>> import fwdpy
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> popsizes=np.tile(popsizes,10000)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> dips = [fwdpy.view_gametes(i) for i in pops]

    Example for a metapopulation:

    >>> from __future__ import print_function
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.001,0.0),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> # Evolve for 5N generations initially
    >>> popsizes=np.tile(popsizes,10000)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> #Now, "bud" off a daughter population of same size, and evolve both for another 100 generations
    >>> mpops = fwdpy.evolve_regions_split(rng,pops,popsizes[0:100],popsizes[0:100],0.001,0.0001,0.001,nregions,sregions,rregions,[0]*2)
    >>> gams = [fwdpy.view_gametes(i,0) for i in mpops]
    >>> #The sum of the gamete counts must be 2*(deme size):
    >>> n=0
    >>> for i in gams[0]: n += i['n']
    >>> n
    2000
    >>> #Exceptions will be thrown if deme is out of range:
    >>> gams = [fwdpy.view_gametes(i,2) for i in mpops]
    Traceback (most recent call last):
     ...
    IndexError: view_gametes: deme index out of ramge
    """
    if isinstance(p,singlepop):
        return view_gametes_singlepop(p)
    elif isinstance(p,metapop):
        return view_gametes_metapop(p,deme)
    else:
        raise RuntimeError("view_gametes: unsupported poptype")

def view_diploids_singlepop( singlepop p, list indlist ):
    return view_diploids_details(p.pop.get().diploids,indlist)
    
def view_diploids_metapop( metapop p, list indlist, unsigned deme ):
    if deme >= len(p.popsizes()):
        raise RuntimeError("view_diploids: deme index out of range")
    return view_diploids_details(p.mpop.get().diploids[deme],indlist)
    
def view_diploids( poptype p, list indlist, deme = None ):
    """
    Get detailed list of a set of diploids in the population

    :param p: a :class:`fwdpy.fwdpy.poptype`

    :rtype: a list of dictionaries.  See Note.

    Example:

    >>> import fwdpy
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> popsizes=np.tile(popsizes,10000)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> dips = [fwdpy.view_diploids(i,[0,101,201,301]) for i in pops]
    """
    if isinstance(p,singlepop):
        return view_diploids_singlepop(p,indlist)
    elif isinstance(p,metapop):
        if deme is None:
            raise RuntimeError("view_diploids: deme index required for metapopulation")
        return view_diploids_metapop(p,indlist,deme)
    else:
        raise RuntimeError("view_diploids: unsupported poptype")
