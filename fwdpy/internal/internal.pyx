# distutils: language = c++
# distutils: sources = fwdpy/internal/callbacks.cc

include "sregionCallbacks.pyx"
include "slim.pyx"

cdef class region_manager_wrapper:
    def __cinit__(self):
        self.thisptr = new region_manager()
    def __dealloc__(self):
        del self.thisptr
        
##Returns a region_manager_wrapper
def make_region_manager(list nregions,
                        list sregions,
                        list recregions):
    nreg = process_regions(nregions)
    sreg = process_regions(sregions)
    recreg = process_regions(recregions)
    v = shwrappervec()
    process_sregion_callbacks(v,sregions)
    rv = region_manager_wrapper()
    rv.thisptr.callbacks = v.vec
    rv.thisptr.nb = nreg['beg'].tolist()
    rv.thisptr.ne = nreg['end'].tolist()
    rv.thisptr.nw = nreg['weight'].tolist()
    rv.thisptr.sb = sreg['beg'].tolist()
    rv.thisptr.se = sreg['end'].tolist()
    rv.thisptr.sw = sreg['weight'].tolist()
    rv.thisptr.rb = recreg['beg'].tolist()
    rv.thisptr.re = recreg['end'].tolist()
    rv.thisptr.rw = recreg['weight'].tolist()
    return rv;
