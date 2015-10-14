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
    return {'neutral':neutral,'selected':selected}

cdef get_diploid( const vector[diploid_t].iterator itr ):
    return {'chrom0':get_gamete(deref(itr).first),
            'chrom1':get_gamete(deref(itr).second)}

cdef view_mutations_details(cpplist[popgenmut].iterator beg,cpplist[popgenmut].iterator end):
    rv=[]
    while beg != end:
        rv.append(get_mutation(beg))
        inc(beg)
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
    if isinstance(p,singlepop):
        return view_mutations_singlepop(p)
    elif isinstance(p,metapop):
        return view_mutations_metapop(p)
    else:
        raise RuntimeError("view_mutations: unsupported poptype")
    
cdef view_gametes_details( cpplist[gamete_t].iterator beg,
                           cpplist[gamete_t].iterator end ):
    rv=[]
    while beg != end:
        print "here!"
        rv.append(get_gamete(beg))
        inc(beg)
    return rv

def view_gametes_singlepop( singlepop p ):
    cdef cpplist[gamete_t].iterator beg = p.pop.get().gametes.begin()
    cdef cpplist[gamete_t].iterator end = p.pop.get().gametes.end()
    return view_gametes_details(beg,end)

def view_gametes_metapop( metapop p ):
    cdef cpplist[gamete_t].iterator beg = p.mpop.get().gametes.begin()
    cdef cpplist[gamete_t].iterator end = p.mpop.get().gametes.end()
    return view_gametes_details(beg,end)

def view_gametes( poptype p ):
    if isinstance(p,singlepop):
        return view_gametes_singlepop(p)
    elif isinstance(p,metapop):
        return view_gametes_metapop(p)
    else:
        raise RuntimeError("view_gametes: unsupported poptype")

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

def view_diploids_singlepop( singlepop p, list indlist ):
    return view_diploids_details(p.pop.get().diploids,indlist)
    
def view_diploids_metapop( metapop p, list indlist, unsigned deme ):
    if deme >= len(p.popsizes()):
        raise RuntimeError("view_diploids: deme index out of range")
    return view_diploids_details(p.mpop.get().diploids[deme],indlist)
    
def view_diploids( poptype p, list indlist, deme = None ):
    if isinstance(p,singlepop):
        return view_diploids_singlepop(p,indlist)
    elif isinstance(p,metapop):
        if deme is None:
            raise RuntimeError("view_diploids: deme index required for metapopulation")
        return view_diploids_metapop(p,indlist,deme)
    else:
        raise RuntimeError("view_diploids: unsupported poptype")
