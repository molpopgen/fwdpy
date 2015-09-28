# distutils: language = c++
# distutils: sources = fwdpy/internal/callbacks.cc

include "sregionCallbacks.pyx"
include "slim.pyx"

cdef class region_manager_wrapper:
    def __cinit__(self):
        self.rm = new region_manager()
    def __dealloc__(self):
        del self.rm
        
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
    rv.rm.callbacks = v.vec
    rv.rm.nb = nreg['beg'].tolist()
    rv.rm.ne = nreg['end'].tolist()
    rv.rm.nw = nreg['weight'].tolist()
    rv.rm.sb = sreg['beg'].tolist()
    rv.rm.se = sreg['end'].tolist()
    rv.rm.sw = sreg['weight'].tolist()
    rv.rm.rb = sreg['beg'].tolist()
    rv.rm.re = sreg['end'].tolist()
    rv.rm.rw = sreg['weight'].tolist()
    return rv;
