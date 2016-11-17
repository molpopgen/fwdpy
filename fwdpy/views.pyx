from cython.operator import dereference as deref,postincrement as inc
from libcpp.limits cimport numeric_limits
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

cdef class MutationView(object):
    def __cinit__(self,float pos,uint32_t n,uint32_t g,ftime, float s,
            float h, bint neutral,uint16_t label,key):
        self.pos=pos
        self.n=n
        self.g=g
        self.ftime=ftime
        self.s=s
        self.h=h
        self.neutral=neutral
        self.label=label
        self.key=key
    def __repr__(self):
        r = b'position: ' + format(self.pos) + b', '
        r += b'count: ' + format(self.n) + b', '
        r += b'origin time: ' + format(self.g) + b', '
        r += b'fixation time: ' + format(self.ftime) + b', '
        r += b's: ' + format(self.s) + b', '
        r += b'h: ' + format(self.h) + b', '
        r += b'neutral: ' + format(self.neutral) + b', '
        r += b'label: ' + format(self.label) + b', '
        r += b'key: ' + format(self.key)
        return r
    def as_dict(self):
        return {'position':self.pos,'count':self.n,'origin':self.g,
                'fixation':self.ftime,'s':self.s,'h':self.h,'neutral':self.neutral,
                'label':self.label,'key':self.key}
cdef class GameteView(object):
    def __cinit__(self,list neutral_mutations,list selected_mutations, int count,key):
        self.neutral=neutral_mutations
        self.selected=selected_mutations
        self.n=count
        self.key=key
    def __repr__(self):
        r = b'Neutral variants:\n' + b'\n'.join([str(i) for i in self.neutral]) + b'\n'
        r += b'Selected variants:\n' + b'\n'.join([str(i) for i in self.selected]) + b'\n'
        return r

cdef class DiploidView(object):
    def __cinit__(self,GameteView a,GameteView b,float genetic_value,float env_value,float fitness,key):
        self.first=a
        self.second=b
        self.g=genetic_value
        self.e=env_value
        self.w=fitness
        self.key=key
cdef class MultiLocusDiploidView(object):
    def __cinit__(self,list a,list b,float genetic_value, float env_value, float fitness,key):
        self.first=a
        self.second=b
        self.g=genetic_value
        self.e=env_value
        self.w=fitness
        self.key=key
        
include "view_mutations.pyx"
include "view_fixations.pyx"
include "view_gametes.pyx"
include "view_diploids.pyx"
include "view_diploid_traits.pyx"

cdef diploid_view_to_sample_init_containers(list mutations, map[double,string] * rmap, list info, const size_t ttl_nsam):
    cdef map[double,string].iterator map_itr
    cdef char ancestral = '0'
    rmap_info = []
    for mut in mutations:
        map_itr = rmap.find(mut['pos'])
        if map_itr == rmap.end():
            rmap.insert(pair[double,string](mut['pos'],string(ttl_nsam,ancestral)))
            info.append(mut)
    return info

cdef diploid_view_to_sample_fill_containers_details(list mutations, map[double,string] * data, const unsigned offset):
    cdef char derived = '1'
    cdef map[double,string].iterator map_itr
    for mut in mutations:
        map_itr = data.find(mut['pos'])
        if map_itr == data.end():
            raise RuntimeError("diploid_view_to_sample_fill_containers_details: fatal error due to mutation position not being found")
        if offset >= deref(map_itr).second.size():
            raise IndexError("diploid_view_to_sample_fill_containers_details: offset out of range")
        deref(map_itr).second[offset]=derived

cdef diploid_view_to_sample_fill_containers(list view, map[double,string] * neutral, map[double,string] * selected): 
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

cdef diploid_mloc_view_to_sample_fill_containers(list view, map[double,string] * neutral, map[double,string] * selected): 
    cdef unsigned I = 0
    cdef size_t ttl_nsam = 2*len(view)
    for dip in view:
        for i in range(len(dip['chrom0'])):
            diploid_view_to_sample_fill_containers_details(dip['chrom0'][i]['neutral'],neutral,2*I)
            diploid_view_to_sample_fill_containers_details(dip['chrom1'][i]['neutral'],neutral,2*I+1)
            diploid_view_to_sample_fill_containers_details(dip['chrom0'][i]['selected'],selected,2*I)
            diploid_view_to_sample_fill_containers_details(dip['chrom1'][i]['selected'],selected,2*I+1)
        I+=1
    if <size_t>I != len(view):
        raise RuntimeError("diploid_mloc_view_to_sample_fill_containers: indexing incorrect")

                  
ctypedef vector[pair[double,string]].iterator vec_itr
ctypedef vector[pair[double,string]] fwdpy_sample_t
ctypedef back_insert_iterator[fwdpy_sample_t] back_insert_itr

#This doesn't appear to actually move, but it is cool that it compiles!
cdef copy_map(map[double,string] & m, vector[pair[double,string]] & v):
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

def diploid_mloc_view_to_sample(list view):
    """
    Convert a "view" of diploids into the data types returned from :func:`fwdpy.fwdpy.get_samples`
    
    :param view: An object returned by :func:`fwdpy.fwdpy.view_diploids` from a multilocus population.

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
        for i in range(len(dip['chrom0'])):
            neutral_info = diploid_view_to_sample_init_containers(dip['chrom0'][i]['neutral'],&neutral,neutral_info,ttl_nsam)
            neutral_info = diploid_view_to_sample_init_containers(dip['chrom1'][i]['neutral'],&neutral,neutral_info,ttl_nsam)
            selected_info = diploid_view_to_sample_init_containers(dip['chrom0'][i]['selected'],&selected,selected_info,ttl_nsam)
            selected_info = diploid_view_to_sample_init_containers(dip['chrom1'][i]['selected'],&selected,selected_info,ttl_nsam)
    cdef size_t size_ = len(neutral_info)
    if neutral.size() != size_:
        raise RuntimeError("diploid_view: unequal container sizes for neutral mutations")
    size_ = len(selected_info)
    if selected.size() != size_:
        raise RuntimeError("diploid_view: unequal container sizes for selected mutations")
    diploid_mloc_view_to_sample_fill_containers(view,&neutral,&selected)

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

cdef void resize_dip_view_data(diploid_view_data & dv, size_t nr) :
    dv.s.resize(nr)
    dv.h.resize(nr)
    dv.n.resize(nr)
    dv.g.resize(nr)
    dv.ind.resize(nr)
    dv.chrom.resize(nr)

cdef size_t fill_dip_view_data(vector[popgen_mut_data].iterator gbeg,
                               vector[popgen_mut_data].iterator gend,
                               diploid_view_data & rv,
                               const size_t ROW,
                               const size_t IND,
                               const size_t ch) :
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
                                                bint selectedOnly) :
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

def view_diploids_pd_popvec(SpopVec p,
                            vector[unsigned] & indlist,
                            bint selectedOnly):
    cdef size_t npops = p.pops.size()
    cdef int i
    cdef vector[diploid_view_data] rv
    rv.resize(npops)
    for i in range(npops):
        rv[i]=view_diploids_pd_details(p.pops[i].get(),indlist,selectedOnly)
    return rv

def view_diploids_pd_singlepop(Spop p,
                               vector[unsigned] & indlist,
                               bint selectedOnly):
    return view_diploids_pd_details(p.pop.get(),indlist,selectedOnly)

def view_diploids_pd(object p,
                     list indlist,
                     bint selectedOnly = True):
    """
    Get detailed list of a set of diploids in the population

    :param p: a :class:`fwdpy.fwdpy.PopType` or a :class:`fwdpy.fwdpy.PopVec`
    :param deme: if p is a :class`fwdpy.fwdpy.MetaPop`, deme is the index of the deme to sample
    
    :rtype: pandas.DataFrame or a list of such objects

    .. note:: :class:`fwdpy.fwdpy.MetaPopVec` is not yet supported
    """
    if isinstance(p,SpopVec):
        return [pd.DataFrame(i) for i in view_diploids_pd_popvec(p,indlist,selectedOnly)]
    elif isinstance(p,Spop):
        return pd.DataFrame(view_diploids_pd_singlepop(p,indlist,selectedOnly))

