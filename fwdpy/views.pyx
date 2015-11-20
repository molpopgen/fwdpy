from cython.operator import dereference as deref,postincrement as inc
from cython.parallel import parallel, prange
import pandas as pd

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

    
cdef popgen_mut_data get_mutation( const mlist_t_itr & itr) nogil:
    cdef popgen_mut_data rv
    rv.pos=deref(itr).pos
    rv.n=deref(itr).n
    rv.g=deref(itr).g
    rv.s=deref(itr).s
    rv.h=deref(itr).h
    rv.neutral==deref(itr).neutral
    return rv;

cdef gamete_data get_gamete( const cpplist[gamete_t].iterator & itr ) nogil:
    cdef vector[mlist_t_itr].iterator beg = deref(itr).mutations.begin()
    cdef vector[mlist_t_itr].iterator end = deref(itr).mutations.end()
    #neutral = []
    #selected = []
    cdef gamete_data rv;
    while beg != end:
        rv.neutral.push_back(get_mutation(deref(beg)))
        inc(beg)

    beg = deref(itr).smutations.begin()
    end = deref(itr).smutations.end()
    while beg != end:
        rv.selected.push_back(get_mutation(deref(beg)))
        inc(beg)
    rv.n=deref(itr).n
    return rv
    #return {'n':deref(itr).n,'neutral':neutral,'selected':selected}

cdef diploid_data get_diploid( const dipvector_t_itr & itr ) nogil:
   cdef diploid_data rv
   rv.g=deref(itr).g
   rv.e=deref(itr).e
   rv.w=deref(itr).w
   rv.chrom0=get_gamete(deref(itr).first)
   rv.chrom1=get_gamete(deref(itr).second)
   return rv

cdef vector[popgen_mut_data] view_mutations_details(mlist_t_itr beg,mlist_t_itr end) nogil:
    cdef vector[popgen_mut_data] rv
    while beg != end:
        rv.push_back(get_mutation(beg))
        inc(beg)
    return rv

cdef vector[gamete_data] view_gametes_details( cpplist[gamete_t].iterator beg,
                           cpplist[gamete_t].iterator end ) nogil:
    cdef vector[gamete_data] rv
    while beg != end:
        rv.push_back(get_gamete(beg))
        inc(beg)
    return rv

##This really should be const...
cdef vector[diploid_data] view_diploids_details( vector[diploid_t] & diploids,
                                                 const vector[unsigned] indlist ) nogil:
    cdef dipvector_t_itr itr = diploids.begin()
    cdef vector[diploid_data] rv
    for i in range(indlist.size()):
        #if indlist[i] >= diploids.size():
        #    raise IndexError("view_diploids: index out of range")
        rv.push_back(get_diploid(itr+indlist[i]))
    return rv

def view_mutations_singlepop(singlepop p):
    cdef mlist_t_itr beg = p.pop.get().mutations.begin()
    cdef mlist_t_itr end = p.pop.get().mutations.end()
    cdef vector[popgen_mut_data] rv;
    with nogil:
        rv = view_mutations_details(beg,end)        
    return sorted(view_mutations_details(beg,end),key = lambda x:x['pos'])

def view_mutations_popvec(popvec p):
    cdef mlist_t_itr beg,end
    cdef vector[vector[popgen_mut_data]] rv;
    cdef int npops = p.pops.size(),i
    rv.resize(npops)
    for i in prange(npops,schedule='guided',nogil=True):
        rv[i] = view_mutations_details(p.pops[i].get().mutations.begin(),p.pops[i].get().mutations.end())

    return rv

def view_mutations_metapop(metapop p,unsigned deme):
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

def view_mutations( object p, deme = None ):
    """
    Get detailed list of all mutations in the population

    :param p: a :class:`fwdpy.fwdpy.poptype` or a :class:`fwdpy.fwdpy.popvec`

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
    elif isinstance(p,popvec):
        return view_mutations_popvec(p)
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

def view_gametes_popvec(popvec p):
    cdef glist_t_itr beg,end
    cdef glist_t_itr 
    cdef int npops = p.pops.size(),i
    cdef vector[vector[gamete_data]] temp
    temp.resize(npops)
    for i in prange(npops,schedule='guided',nogil=True):
        temp[i]=view_gametes_details(p.pops[i].get().gametes.begin(),p.pops[i].get().gametes.end())
    rv=[]
    for i in range(npops):
        rv.append(temp[i])
        rv[i]=sorted(rv[i],key=lambda x:x['n'],reverse=True)

    return rv
        
    
    

def view_gametes_metapop( metapop p, unsigned deme ):
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

def view_gametes( object p ,deme = None):
    """
    Get detailed list of all gametes in the population

    :param p: a :class:`fwdpy.fwdpy.poptype` or a :class:`fwdpy.fwdpy.popvec`
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
    elif isinstance(p,popvec):
        return view_gametes_popvec(p)
    elif isinstance(p,metapop):
        if deme is None:
            raise RuntimeError("view_gametes: deme cannot be None when p is a metapop")
        return view_gametes_metapop(p,deme)
    else:
        raise RuntimeError("view_gametes: unsupported poptype")

def view_diploids_singlepop( singlepop p, list indlist ):
    for i in indlist:
        if i >= p.popsize():
            raise IndexError("index greater than population size")
    return view_diploids_details(p.pop.get().diploids,indlist)

def view_diploids_popvec( popvec p, list indlist ):
    cdef int npops = len(p),i
    cdef vector[vector[diploid_data]] rv
    rv.resize(npops);
    cdef vector[unsigned] il
    for i in indlist:
        il.push_back(i)
    for i in prange(npops,schedule='guided',nogil=True):
        rv[i] = view_diploids_details(p.pops[i].get().diploids,il)

    return rv
        
def view_diploids_metapop( metapop p, list indlist, unsigned deme ):
    psizes = p.popsizes()
    for i in indlist:
        for ps in psizes:
            if i >= ps:
                raise IndexError("index greater than deme size")
    if deme >= len(p.popsizes()):
        raise IndexError("view_diploids: deme index out of range")
    return view_diploids_details(p.mpop.get().diploids[deme],indlist)
    
def view_diploids( object p, list indlist, deme = None ):
    """
    Get detailed list of a set of diploids in the population

    :param p: a :class:`fwdpy.fwdpy.poptype` or a :class:`fwdpy.fwdpy.popvec`
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
    elif isinstance(p,popvec):
        return view_diploids_popvec(p,indlist)
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
    """
    Convert a "view" of diploids into the data types returned from :func:`fwdpy.fwdpy.get_samples`
    
    :param view: An object returned by :func:`fwdpy.fwdpy.view_diploids`

    :return: A dict.  See Note

    :rtype: dict

    .. note:: The return value is a dict containing two lists of tuples ('neutral', and 'selected')  These lists are the same format 
       as objects returned by :func:`fwdpy.fwdpy.get_samples`.  The dict also contains 'neutral_info' and 'selected_info', which are
       pandas.DataFrame objects with the information about the mutations present in the other two lists.

    """
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

    neutral_df = pd.DataFrame()
    if len(neutral_info)>0:
        neutral_df = pd.concat([pd.DataFrame(i,index=[0]) for i in neutral_info])
        neutral_df['index']=list(range(len(neutral_df.index)))
        neutral_df = neutral_df.set_index('index')

    selected_df = pd.DataFrame()
    if len(selected_info)>0:
        selected_df = pd.concat([pd.DataFrame(i,index=[0]) for i in selected_info])
        selected_df['index']=list(range(len(selected_df.index)))
        selected_df = selected_df.set_index('index')

    return {'neutral':vneutral,'selected':vselected,'neutral_info':neutral_df,'selected_info':selected_df}
