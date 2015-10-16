from cython.operator import dereference as deref,postincrement as inc

cdef extern from "<algorithm>" namespace "std":
    OUTPUT move[INPUT,OUTPUT](INPUT,INPUT,OUTPUT)
cdef extern from "<iterator>" namespace "std":
    cdef cppclass iterator[Category,T,Distance,Pointer,Reference]:
        pass
    cdef cppclass output_iterator_tag:
        pass
    cdef cppclass back_insert_iterator[T](iterator[output_iterator_tag,void,void,void,void]):
        pass
    back_insert_iterator[CONTAINER] back_inserter[CONTAINER](CONTAINER &)

cdef get_mutation( const mlist_t_itr & itr):
    return {'pos':deref(itr).pos,'n':deref(itr).n,'g':deref(itr).g,'s':deref(itr).s,'h':deref(itr).h,'neutral':deref(itr).neutral}

cdef get_gamete( const cpplist[gamete_t].iterator & itr ):
    cdef vector[mlist_t_itr].iterator beg = deref(itr).mutations.begin()
    cdef vector[mlist_t_itr].iterator end = deref(itr).mutations.end()
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

cdef get_diploid( const dipvector_t_itr & itr ):
    return {'chrom0':get_gamete(deref(itr).first),
            'chrom1':get_gamete(deref(itr).second)}

cdef view_mutations_details(mlist_t_itr beg,mlist_t_itr end):
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
    cdef dipvector_t_itr itr = diploids.begin()
    rv=[]
    for i in range(indlist.size()):
        if indlist[i] >= diploids.size():
            raise IndexError("view_diploids: index out of range")
        rv.append(get_diploid(itr+indlist[i]))
    return rv

def view_mutations_singlepop(singlepop p):
    cdef mlist_t_itr beg = p.pop.get().mutations.begin()
    cdef mlist_t_itr end = p.pop.get().mutations.end()
    return sorted(view_mutations_details(beg,end),key = lambda x:x['pos'])

def view_mutations_metapop(metapop p,deme):
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
    return sorted(rv, key = lambda x : x['pos'])

def view_mutations( poptype p, deme = None ):
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
    >>> popsizes=np.tile(popsizes,100)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> muts = [fwdpy.view_mutations(i) for i in pops]
    >>> type(muts[0])
    <type 'list'>
    >>> type(muts[0][0])
    <type 'dict'>

    Metapopulation example:

    >>> mpops = fwdpy.evolve_regions_split(rng,pops,popsizes[0:100],popsizes[0:100],0.001,0.0001,0.001,nregions,sregions,rregions,[0]*2)
    >>> muts_deme_0 = [fwdpy.view_mutations(i,0) for i in mpops]
    >>> muts_deme_1 = [fwdpy.view_mutations(i,1) for i in mpops]
    
    """
    if isinstance(p,singlepop):
        return view_mutations_singlepop(p)
    elif isinstance(p,metapop):
        if deme is None:
            raise RuntimeError("view_mutations: deme cannot be none for metapops")
        return view_mutations_metapop(p,deme)
    else:
        raise RuntimeError("view_mutations: unsupported poptype")
    
def view_gametes_singlepop( singlepop p ):
    cdef glist_t_itr beg = p.pop.get().gametes.begin()
    cdef glist_t_itr end = p.pop.get().gametes.end()
    return sorted(view_gametes_details(beg,end),key=lambda x:x['n'],reverse=True)

def view_gametes_metapop( metapop p, deme):
    if deme >= len(p.popsizes()):
        raise IndexError("view_gametes: deme index out of ramge")
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
    return sorted(temp1,key=lambda x:x['n'],reverse=True)

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
    >>> popsizes=np.tile(popsizes,100)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> dips = [fwdpy.view_gametes(i) for i in pops]

    Example for a metapopulation:

    >>> mpops = fwdpy.evolve_regions_split(rng,pops,popsizes[0:100],popsizes[0:100],0.001,0.0001,0.001,nregions,sregions,rregions,[0]*2)
    >>> gams = [fwdpy.view_gametes(i,0) for i in mpops]
    """
    if isinstance(p,singlepop):
        return view_gametes_singlepop(p)
    elif isinstance(p,metapop):
        if deme is None:
            raise RuntimeError("view_gametes: deme cannot be None when p is a metapop")
        return view_gametes_metapop(p,deme)
    else:
        raise RuntimeError("view_gametes: unsupported poptype")

def view_diploids_singlepop( singlepop p, list indlist ):
    return view_diploids_details(p.pop.get().diploids,indlist)
    
def view_diploids_metapop( metapop p, list indlist, unsigned deme ):
    if deme >= len(p.popsizes()):
        raise IndexError("view_diploids: deme index out of range")
    return view_diploids_details(p.mpop.get().diploids[deme],indlist)
    
def view_diploids( poptype p, list indlist, deme = None ):
    """
    Get detailed list of a set of diploids in the population

    :param p: a :class:`fwdpy.fwdpy.poptype`
    :param deme: if p is a :class`fwdpy.fwdpy.metapop`, deme is the index of the deme to sample
    
    :rtype: a list of dictionaries.  See Note.

    Example:

    >>> import fwdpy
    >>> import numpy as np
    >>> nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
    >>> sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
    >>> rregions = [fwdpy.Region(0,3,1)]
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes = np.array([1000],dtype=np.uint32)
    >>> popsizes=np.tile(popsizes,100)
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> dips = [fwdpy.view_diploids(i,[0,101,201,301]) for i in pops]

    And now viewing from a metapop:
    
    >>> pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions)
    >>> #Now, "bud" off a daughter population of same size, and evolve both for another 100 generations
    >>> mpops = fwdpy.evolve_regions_split(rng,pops,popsizes[0:100],popsizes[0:100],0.001,0.0001,0.001,nregions,sregions,rregions,[0]*2)
    >>> dips = [fwdpy.view_diploids(i,[0,101,201,301],0) for i in mpops]
    """
    if isinstance(p,singlepop):
        return view_diploids_singlepop(p,indlist)
    elif isinstance(p,metapop):
        if deme is None:
            raise RuntimeError("view_diploids: deme index required for metapopulation")
        return view_diploids_metapop(p,indlist,deme)
    else:
        raise RuntimeError("view_diploids: unsupported poptype")

cdef diploid_view_to_sample_init_containers( list mutations, map[double,string] * rmap, list info, const size_t ttl_nsam ):
    cdef map[double,string].iterator map_itr
    cdef char ancestral = '0'
    rmap_info = []
    for mut in mutations:
        map_itr = rmap.find(mut['pos'])
        if map_itr == rmap.end():
            rmap.insert( pair[double,string](mut['pos'],string(ttl_nsam,ancestral)) )
            info.append(mut)
    return info

cdef diploid_view_to_sample_fill_containers_details( list mutations, map[double,string] * data, const unsigned offset ):
    cdef char derived = '1'
    cdef map[double,string].iterator map_itr
    for mut in mutations:
        map_itr = data.find(mut['pos'])
        if map_itr == data.end():
            raise RuntimeError("diploid_view_to_sample_fill_containers_details: fatal error due to mutation position not being found")
        if offset >= deref(map_itr).second.size():
            raise IndexError("diploid_view_to_sample_fill_containers_details: offset out of range")
        deref(map_itr).second[offset]=derived

cdef diploid_view_to_sample_fill_containers( list view, map[double,string] * neutral, map[double,string] * selected ): 
    cdef unsigned I = 0
    cdef size_t ttl_nsam = 2*len(view)
    for dip in view:
        diploid_view_to_sample_fill_containers_details(dip['chrom0']['neutral'],neutral,2*I)
        diploid_view_to_sample_fill_containers_details(dip['chrom1']['neutral'],neutral,2*I+1)
        diploid_view_to_sample_fill_containers_details(dip['chrom0']['selected'],selected,2*I)
        diploid_view_to_sample_fill_containers_details(dip['chrom1']['selected'],selected,2*I+1)
        I+=1
    if I != len(view):
        raise RuntimeError("diploid_view_to_sample_fill_containers: indexing incorrect")
                        
ctypedef vector[pair[double,string]].iterator vec_itr
ctypedef vector[pair[double,string]] fwdpy_sample_t
ctypedef back_insert_iterator[fwdpy_sample_t] back_insert_itr

#This doesn't appear to actually move, but it is cool that it compiles!
cdef copy_map( map[double,string] & m, vector[pair[double,string]] & v):
    move[map[double,string].iterator,back_insert_itr](m.begin(),m.end(),back_inserter[fwdpy_sample_t](v))

def diploid_view_to_sample(list view):
    cdef size_t ttl_nsam = 2*len(view)
    if ttl_nsam == 0:
        return []

    cdef map[double,string] neutral
    cdef map[double,string] selected
    cdef map[double,string].iterator map_itr
    neutral_info = []
    selected_info = []
    cdef char ancestral = '0'
    cdef char derived = '0'
    for dip in view:
        neutral_info = diploid_view_to_sample_init_containers(dip['chrom0']['neutral'],&neutral,neutral_info,ttl_nsam)
        neutral_info = diploid_view_to_sample_init_containers(dip['chrom1']['neutral'],&neutral,neutral_info,ttl_nsam)
        selected_info = diploid_view_to_sample_init_containers(dip['chrom0']['selected'],&selected,selected_info,ttl_nsam)
        selected_info = diploid_view_to_sample_init_containers(dip['chrom1']['selected'],&selected,selected_info,ttl_nsam)
    if neutral.size() != len(neutral_info):
        raise RuntimeError("diploid_view: unequal container sizes for neutral mutations")
    if selected.size() != len(selected_info):
        raise RuntimeError("diploid_view: unequal container sizes for selected mutations")
    diploid_view_to_sample_fill_containers(view,&neutral,&selected)

    ##Now, we convert from maps to vectors.
    ##We attempt to do this with C++11 move, which may or may not happen...
    cdef vector[pair[double,string]] vneutral
    cdef vector[pair[double,string]] vselected
    copy_map(neutral,vneutral)
    copy_map(selected,vselected)

    return {'neutral':vneutral,'selected':vselected,'neutral_info':neutral_info,'selected_info':selected_info}
