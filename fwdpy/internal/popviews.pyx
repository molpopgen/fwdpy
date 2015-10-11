from cython.operator cimport dereference as deref,preincrement as inc
from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.list cimport list as cpplist
from fwdpy.fwdpy cimport popgenmut,dipvector_t

cdef fill_diploid_view_details( map[string,vector[double]] & rv,
                                const vector[cpplist[popgenmut].iterator] & muts,
                                const size_t ind,
                                const unsigned twoN,
                                const int remove_fixed,
                                const int chrom ):
    cdef vector[cpplist[popgenmut].iterator].const_iterator beg = muts.const_begin()
    cdef vector[cpplist[popgenmut].iterator].const_iterator end = muts.const_end()
    cdef cpplist[popgenmut].iterator mitr
    while beg != end:
        mitr = deref(beg)
        p = float(deref(mitr).n)/float(twoN)
        if p < 1. or remove_fixed == 0:
            rv["ind"].push_back(ind);
            rv["pos"].push_back(deref(mitr).pos);
            rv["esize"].push_back(deref(mitr).s);
            rv["h"].push_back(deref(mitr).h);
            rv["p"].push_back(p);
            rv["hap"].push_back(float(chrom));
            rv["origin"].push_back(float(deref(mitr).g));
	  
cdef diploid_view_details(const dipvector_t & diploids,
                          const size_t ind,
                          const int remove_fixed):
    if ind > diploids.size():
        raise RuntimeError("diploid_view_details: individual index out of range");
    twoN = 2*diploids.size()
    rvlabels = ["ind","pos","esize","h","p","hap","origin"]
    cdef map[string,vector[double]] rv
    for i in rvlabels:
        rv[i] = vector[double]()
    if deref(diploids[ind].first).mutations.empty():
        if deref(diploids[ind].first).smutations.empty():
            if deref(diploids[ind].second).mutations.empty():
                if deref(diploids[ind].second).smutations.empty():
                    for i in rvlabels:
                        if i == 'ind':
                            rv[i].push_back(ind)
                        else:
                        #WARNING: there's a big assumption here!
                            rv[i].push_back(np.nan)
                    return rv
    
def diploid_view_singlepop(singlepop pop, int ind, bint removeFixed):
    return diploid_view_details(pop.pop.get().diploids,ind,removeFixed)

def diploid_view_metapop(metapop mpop, int ind, bint removeFixed, int deme):
    return diploid_view_details(mpop.mpop.get().diploids[deme],ind,removeFixed)
