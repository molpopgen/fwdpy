# distutils: language = c++
# distutils: sources = fwdpy/qtrait/qtrait_impl.cc fwdpy/qtrait/ew2010.cc
from libcpp.vector cimport vector
from libcpp.map cimport map
from fwdpy.fwdpy cimport *
from fwdpy.internal.internal cimport shwrappervec
import fwdpy.internal as internal
import pandas

cdef extern from "qtraits.hpp" namespace "fwdpy::qtrait":
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
        const vector[double] & rweight)
    map[string,double] qtrait_pop_props( const singlepop_t * pop );
    map[string,vector[double]] get_qtrait_traj(const singlepop_t *pop,const unsigned minsojourn,const double minfreq)
    map[string,vector[double]] qtrait_esize_freq(const singlepop_t * pop)
    vector[double] ew2010_traits_cpp(GSLrng_t * rng, const singlepop_t * pop, const double tau, const double sigma)
    
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
                    double f = 0.):
    pops = popvec(npops,N)
    nreg = internal.process_regions(nregions)
    sreg = internal.process_regions(sregions)
    recreg = internal.process_regions(recregions)
    v = shwrappervec()
    internal.process_sregion_callbacks(v,sregions)
    evolve_qtraits_t(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,optimum,track,
                    nreg['beg'].tolist(),nreg['end'].tolist(),nreg['weight'].tolist(),
                    sreg['beg'].tolist(),sreg['end'].tolist(),sreg['weight'].tolist(),&v.vec,
                    recreg['beg'].tolist(),recreg['end'].tolist(),recreg['weight'].tolist())
    return pops

def evolve_qtrait_more(GSLrng rng,
                    popvec pops,
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
                    double f = 0.):
    nreg = internal.process_regions(nregions)
    sreg = internal.process_regions(sregions)
    recreg = internal.process_regions(recregions)
    v = shwrappervec()
    internal.process_sregion_callbacks(v,sregions)
    evolve_qtraits_t(rng.thisptr,&pops.pops,&nlist[0],len(nlist),mu_neutral,mu_selected,recrate,f,sigmaE,optimum,track,
                    nreg['beg'].tolist(),nreg['end'].tolist(),nreg['weight'].tolist(),
                    sreg['beg'].tolist(),sreg['end'].tolist(),sreg['weight'].tolist(),&v.vec,
                    recreg['beg'].tolist(),recreg['end'].tolist(),recreg['weight'].tolist())

def popstats( singlepop pop ):
    return pandas.DataFrame(qtrait_pop_props(pop.pop.get()).items(),columns=['stat','value'])

def trajectories( singlepop pop, int minsojourn = 0, double minfreq = 0.):
    return pandas.DataFrame.from_dict(get_qtrait_traj(pop.pop.get(),minsojourn,minfreq))

def esize_freq(singlepop pop):
    return pandas.DataFrame.from_dict(qtrait_esize_freq(pop.pop.get()))

def ew2010_traits(GSLrng rng, singlepop pop, double tau, double sigma):
    return ew2010_traits_cpp(rng.thisptr,pop.pop.get(),tau,sigma)
