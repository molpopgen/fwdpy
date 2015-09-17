# distutils: language = c++
# distutils: sources = fwdpy/qtrait/qtrait_impl.cc
from libcpp.vector cimport vector
from fwdpy.fwdpy cimport *
from fwdpy.internal.internal cimport shwrappervec
import fwdpy.internal as internal

cdef extern from "evolve_qtraits.hpp" namespace "fwdpy::qtrait":
    void evolve_qtraits_t( GSLrng_t * rng, vector[shared_ptr[singlepop_t] ] * pops,
        const unsigned * Nvector,
        const size_t Nvector_length,
        const double mu_neutral,
        const double mu_selected,
        const double littler,
        const double f,
        const double sigmaE,
        const double optimum,
        const int track,
        const vector[double] & nbegs,
        const vector[double] & nends,
        const vector[double] & nweights,
        const vector[double] & sbegs,
        const vector[double] & sends,
        const vector[double] & sweights,
        const vector[shmodel] * callbacks,
        const vector[double] & rbeg,
        const vector[double] & rend,
        const vector[double] & rweight,
        const char * fitness)

def evolve_qtrait(GSLrng rng,
                    int npops,
                    int N,
                    unsigned[:] nlist,
                    double mu_neutral,
                    double mu_selected,
                    double recrate,
                    list nregions,
                    list sregions,
                    list recregions,
                    double sigmaE,
                    double optimum = 0.,
                    bint track = False,
                    double f = 0.,
                    const char * fitness = "multiplicative"):
    pops = popvec(npops,N)
    nreg = internal.process_regions(nregions)
    sreg = internal.process_regions(sregions)
    recreg = internal.process_regions(recregions)
    v = shwrappervec()
    internal.process_sregion_callbacks(v,sregions)
    evolve_qtraits_t(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,optimum,track,
                    nreg['beg'].tolist(),nreg['end'].tolist(),nreg['weight'].tolist(),
                    sreg['beg'].tolist(),sreg['end'].tolist(),sreg['weight'].tolist(),&v.vec,
                    recreg['beg'].tolist(),recreg['end'].tolist(),recreg['weight'].tolist(),
                    fitness)
    return pops
