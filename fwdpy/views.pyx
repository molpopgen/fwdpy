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

    
cdef popgen_mut_data get_mutation( const popgenmut & m, size_t n) nogil:
    cdef popgen_mut_data rv
    rv.pos=m.pos
    rv.n=<unsigned>n
    rv.g=m.g
    rv.s=m.s
    rv.h=m.h
    rv.neutral=m.neutral
    return rv;

cdef gamete_data get_gamete( const gamete_t & g, const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil:
    cdef gamete_data rv;
    cdef size_t i=0,j=g.mutations.size()
    while i<j:
        rv.neutral.push_back(get_mutation(mutations[i],mcounts[i]))
        i+=1

    i=0
    j=g.smutations.size()
    while i<j:
        rv.selected.push_back(get_mutation(mutations[i],mcounts[i]))
        i+=1
    rv.n=g.n
    return rv

cdef diploid_data get_diploid( const diploid_t & dip, const gcont_t & gametes, const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil:
   cdef diploid_data rv
   rv.g=dip.g
   rv.e=dip.e
   rv.w=dip.w
   rv.chrom0=get_gamete(gametes[dip.first],mutations,mcounts)
   rv.chrom1=get_gamete(gametes[dip.second],mutations,mcounts)
   rv.n0 = <unsigned>rv.chrom0.selected.size()
   rv.n1 = <unsigned>rv.chrom1.selected.size()
   cdef size_t i = 0
   while i < rv.chrom0.selected.size():
       rv.sh0+=(rv.chrom0.selected[i].s*rv.chrom0.selected[i].h)
       i+=1
   i=0
   while i < rv.chrom1.selected.size():
       rv.sh1+=(rv.chrom1.selected[i].s*rv.chrom1.selected[i].h)
       i+=1
   return rv

cdef vector[popgen_mut_data] view_mutations_details(const mcont_t & mutations, const mcounts_cont_t & mcounts) nogil:
    cdef vector[popgen_mut_data] rv
    cdef size_t i=0,j=mutations.size();
    while i!=j:
        #skip extinct mutation
        if mcounts[i]:
            rv.push_back(get_mutation(mutations[i],mcounts[i]))
        i+=1
    return rv

cdef vector[gamete_data] view_gametes_details( const singlepop_t * pop ) nogil:
    cdef vector[gamete_data] rv
    cdef size_t i = 0, j = pop.gametes.size()
    while i!=j:
        #skip extinct gamets
        if pop.gametes[i].n:
            rv.push_back(get_gamete(pop.gametes[i],pop.mutations,pop.mcounts));
        i+=1
    return rv

##This really should be const...
cdef vector[diploid_data] view_diploids_details( const dipvector_t & diploids,
                                                 const gcont_t & gametes,
                                                 const mcont_t & mutations,
                                                 const mcounts_cont_t & mcounts,
                                                 const vector[unsigned] & indlist ) nogil:
    cdef vector[diploid_data] rv
    for i in range(indlist.size()):
        rv.push_back(get_diploid(diploids[indlist[i]],gametes,mutations,mcounts))
    return rv

def view_mutations_singlepop(singlepop p):
    cdef mcont_t_itr beg = p.pop.get().mutations.begin()
    cdef mcont_t_itr end = p.pop.get().mutations.end()
    cdef vector[popgen_mut_data] rv;
    with nogil:
        rv = view_mutations_details(p.pop.get().mutations,p.pop.get().mcounts)
    return rv

def view_mutations_popvec(popvec p):
    cdef mcont_t_itr beg,end
    cdef vector[vector[popgen_mut_data]] rv;
    cdef size_t npops = p.pops.size()
    cdef int i
    rv.resize(npops)
    #for i in range(npops):
    for i in prange(npops,schedule='static',nogil=True,chunksize=1):
        rv[i] = view_mutations_details(p.pops[i].get().mutations,p.pops[i].get().mcounts)

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
    return rv

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
    
    .. note:: :class:`fwdpy.fwdpy.mpopvec` currently not supported
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
        raise RuntimeError("view_mutations: unsupported object type")
    
def view_gametes_singlepop( singlepop p ):
    #cdef gcont_t_itr beg = p.pop.get().gametes.begin()
    #cdef gcont_t_itr end = p.pop.get().gametes.end()
    return view_gametes_details(p.pop.get())

def view_gametes_popvec(popvec p):
    cdef:
        gcont_t_itr beg,end
        size_t npops = p.pops.size()
        int i
        vector[vector[gamete_data]] rv
    rv.resize(npops)
    #for i in range(npops):
    for i in prange(npops,schedule='static',nogil=True,chunksize=1):
        rv[i]=view_gametes_details(p.pops[i].get())
        
def view_gametes_metapop( metapop p, unsigned deme ):
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

    .. note:: :class:`fwdpy.fwdpy.mpopvec` currently not supported
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
        raise RuntimeError("view_gametes: unsupported object type")

def view_diploids_singlepop( singlepop p, list indlist ):
    for i in indlist:
        if i >= p.popsize():
            raise IndexError("index greater than population size")
    return view_diploids_details(p.pop.get().diploids,p.pop.get().gametes,p.pop.get().mutations,p.pop.get().mcounts,indlist)

def view_diploids_popvec( popvec p, list indlist ):
    cdef size_t npops = len(p),i
    cdef vector[vector[diploid_data]] rv
    rv.resize(npops);
    cdef vector[unsigned] il
    for i in indlist:
        il.push_back(<unsigned>(i))
    #for i in range(npops):
    for i in prange(npops,schedule='static',nogil=True,chunksize=1):
        rv[i] = view_diploids_details(p.pops[i].get().diploids,
                                      p.pops[i].get().gametes,
                                      p.pops[i].get().mutations,
                                      p.pops[i].get().mcounts,il)
    return rv
        
def view_diploids_metapop( metapop p, list indlist, unsigned deme ):
    psizes = p.popsizes()
    for i in indlist:
        for ps in psizes:
            if i >= ps:
                raise IndexError("index greater than deme size")
    if deme >= len(p.popsizes()):
        raise IndexError("view_diploids: deme index out of range")
    return view_diploids_details(p.mpop.get().diploids[deme],p.mpop.get().gametes,p.mpop.get().mutations,p.mpop.get().mcounts,indlist)
    
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

    .. note:: :class:`fwdpy.fwdpy.mpopvec` currently not supported
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
        raise RuntimeError("view_diploids: unsupported object type")

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
    if <size_t>I != len(view):
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
    cdef size_t size_ = len(neutral_info)
    if neutral.size() != size_:
        raise RuntimeError("diploid_view: unequal container sizes for neutral mutations")
    size_ = len(selected_info)
    if selected.size() != size_:
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

cdef struct diploid_view_data:
    vector[double] s,h
    vector[unsigned] n,g,ind,chrom

cdef void resize_dip_view_data( diploid_view_data & dv, size_t nr ) nogil:
    dv.s.resize(nr)
    dv.h.resize(nr)
    dv.n.resize(nr)
    dv.g.resize(nr)
    dv.ind.resize(nr)
    dv.chrom.resize(nr)

cdef size_t fill_dip_view_data( vector[popgen_mut_data].iterator gbeg,
                                  vector[popgen_mut_data].iterator gend,
                                  diploid_view_data & rv,
                                  const size_t ROW,
                                  const size_t IND,
                                  const size_t ch) nogil:
   cdef size_t R=ROW
   while gbeg != gend:
       rv.s[R] = deref(gbeg).s
       rv.h[R] = deref(gbeg).h
       rv.n[R] = deref(gbeg).n
       rv.g[R] = deref(gbeg).g
       rv.ind[R]=<unsigned>IND
       rv.chrom[R]=<unsigned>ch
       R+=1
       inc(gbeg)
   return R

cdef diploid_view_data view_diploids_pd_details(const singlepop_t * pop,
                                                const vector[unsigned] & indlist,
                                                bint selectedOnly) nogil:
    cdef vector[diploid_data] v = view_diploids_details(pop.diploids,
                                                        pop.gametes,
                                                        pop.mutations,
                                                        pop.mcounts,
                                                        indlist)
    cdef vector[diploid_data].iterator beg,end
    beg = v.begin()
    end = v.end()
    cdef size_t nr = 0

    #Determine length of vectors that we'll need
    while beg != end:
        dd = deref(beg)
        if selectedOnly == False:
            nr += deref(beg).chrom0.neutral.size()
            nr += deref(beg).chrom1.neutral.size()
        nr += deref(beg).chrom0.selected.size()
        nr += deref(beg).chrom1.selected.size()
        inc(beg)
    #Construct and size return value
    cdef diploid_view_data rv
    resize_dip_view_data(rv,nr)

    #Fill the return value
    beg = v.begin()
    end = v.end()
    cdef:
        size_t ROW=0
        size_t IND=0
    while beg != end:
        if selectedOnly == False:
            ROW=fill_dip_view_data(deref(beg).chrom0.neutral.begin(),deref(beg).chrom0.neutral.end(),rv,ROW,IND,0)
            ROW=fill_dip_view_data(deref(beg).chrom1.neutral.begin(),deref(beg).chrom1.neutral.end(),rv,ROW,IND,1)
        ROW=fill_dip_view_data(deref(beg).chrom0.selected.begin(),deref(beg).chrom0.selected.end(),rv,ROW,IND,0)
        ROW=fill_dip_view_data(deref(beg).chrom1.selected.begin(),deref(beg).chrom1.selected.end(),rv,ROW,IND,1)
        IND+=1
        inc(beg)

def view_diploids_pd_popvec( popvec p, vector[unsigned] & indlist, bint selectedOnly ):
    cdef size_t npops = p.pops.size()
    cdef int i
    cdef vector[diploid_view_data] rv
    rv.resize(npops)
    for i in prange(npops,schedule='static',nogil=True,chunksize=1):
        rv[i]=view_diploids_pd_details(p.pops[i].get(),indlist,selectedOnly)
    return rv

def view_diploids_pd_singlepop( singlepop p, vector[unsigned] & indlist, bint selectedOnly ):
    return view_diploids_pd_details(p.pop.get(),indlist,selectedOnly)

def view_diploids_pd( object p, list indlist, bint selectedOnly = True ):
    """
    Get detailed list of a set of diploids in the population

    :param p: a :class:`fwdpy.fwdpy.poptype` or a :class:`fwdpy.fwdpy.popvec`
    :param deme: if p is a :class`fwdpy.fwdpy.metapop`, deme is the index of the deme to sample
    
    :rtype: pandas.DataFrame or a list of such objects

    .. note:: :class:`fwdpy.fwdpy.mpopvec` is not yet supported
    """
    if isinstance(p,popvec):
        return [pd.DataFrame(i) for i in view_diploids_pd_popvec(p,indlist,selectedOnly)]
    elif isinstance(p,singlepop):
        return pd.DataFrame( view_diploids_pd_singlepop(p,indlist,selectedOnly) )

cdef diploid_traits_singlepop(singlepop p):
    rv=[]
    for i in range(p.pop.get().diploids.size()):
        rv.append( {'g':p.pop.get().diploids[i].g,
                    'e':p.pop.get().diploids[i].e,
                    'w':p.pop.get().diploids[i].w} )
        
    return rv;

cdef diploid_traits_popvec(popvec p):
    return [diploid_traits_singlepop(i) for i in p]

cdef diploid_traits_mpop(metapop m, deme):
    if deme > m.mpop.get().diploids.size():
        raise RuntimeError("deme value out of range")
    rv=[]
    for i in range(m.mpop.get().diploids[deme].size()):
        rv.append( {'g':m.mpop.get().diploids[deme][i].g,
                    'e':m.mpop.get().diploids[deme][i].e,
                    'w':m.mpop.get().diploids[deme][i].w} )

cdef diploid_traits_mpopvec(mpopvec p,deme):
    return [diploid_traits_mpop(i,deme) for i in p]

def diploid_traits( object p, deme = None ):
    """
    Return genetic value (g), environmental value (e), and fitness (w) for all diploids.

    .. note:: "Standard population genetic" models do not update these values during simulation.
    """
    if isinstance(p,singlepop):
        return diploid_traits_singlepop(p)
    elif isinstance(p,metapop):
        if deme is None:
            raise RuntimeError("deme cannot be None")
        return diploid_traits_mpop(p,deme)
    if isinstance(p,popvec):
        return diploid_traits_popvec(p)
    elif isinstance(p,mpopvec):
        if deme is None:
            raise RuntimeError("deme cannot be None")
        return diploid_traits_mpopvec(p,deme)
    else:
        raise ValueError("unsupported type")    
