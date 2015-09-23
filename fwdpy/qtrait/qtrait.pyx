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
    map[double,double] ew2010_assign_effects(GSLrng_t * rng, const singlepop_t * pop, const double tau, const double sigma) except +
    vector[double] ew2010_traits_cpp(const singlepop_t * pop, const map[double,double] & effects) except +
    
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
    """
    Evolve a quantitative trait with variable mutation, fitness effects, and recombination rates.

    :param rng: a :class:`GSLrng`
    :param npops: The number of populations to simulate.  This is equal to the number of threads that will be used!
    :param N: The diploid population size to simulate
    :param nlist: An array view of a NumPy array.  This represents the population sizes over time.  The length of this view is the length of the simulation in generations. The view must be of an array of 32 bit, unsigned integers (see example).
    :param mu_neutral: The mutation rate to variants not affecting fitness ("neutral" mutations).  The unit is per gamete, per generation.
    :param mu_selected: The mutation rate to variants affecting fitness ("selected" mutations).  The unit is per gamete, per generation.
    :param recrate: The recombination rate in the regions (per diploid, per generation)
    :param nregions: A list specifying where neutral mutations occur
    :param sregions: A list specifying where selected mutations occur
    :param recregions: A list specifying how the genetic map varies along the region
    :param sigmaE: The standard deviation in random variation to add to trait value
    :param optimum: The optimum trait value.
    :param track: whether or not to record the frequency trajectories of mutations.  True = simulations are much slower!
    :param f: The selfing probabilty
    """
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
    """
    Continue to evolve a quantitative trait with variable mutation, fitness effects, and recombination rates.

    :param rng: a :class:`GSLrng`
    :param pops: A list of populations simulated by :func:`fwdpy.qtrait.qtrait.evolve_qtrait`
    :param N: The diploid population size to simulate
    :param nlist: An array view of a NumPy array.  This represents the population sizes over time.  The length of this view is the length of the simulation in generations. The view must be of an array of 32 bit, unsigned integers (see example).
    :param mu_neutral: The mutation rate to variants not affecting fitness ("neutral" mutations).  The unit is per gamete, per generation.
    :param mu_selected: The mutation rate to variants affecting fitness ("selected" mutations).  The unit is per gamete, per generation.
    :param recrate: The recombination rate in the regions (per diploid, per generation)
    :param nregions: A list specifying where neutral mutations occur
    :param sregions: A list specifying where selected mutations occur
    :param recregions: A list specifying how the genetic map varies along the region
    :param sigmaE: The standard deviation in random variation to add to trait value
    :param optimum: The optimum trait value.
    :param track: whether or not to record the frequency trajectories of mutations.  True = simulations are much slower!
    :param f: The selfing probabilty
    """
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
    """
    Get statistics like VG, etc. from a population

    :return: A pandas.DataFrame.

    :rtype: pandas.DataFrame
    """
    return pandas.DataFrame(qtrait_pop_props(pop.pop.get()).items(),columns=['stat','value'])

def trajectories( singlepop pop, int minsojourn = 0, double minfreq = 0.):
    """
    Get the frequency trajactories of mutations

    :param pop: A :class:`fwdpy.fwdpy.singlepop` simulated using :func:`fwdpy.qtrait.qtrait.evolve_qtrait` and/or :func:`fwdp.qtrait.qtrait.evolve_qtrait_more`
    :param minsojounrn: Exclude all mutations that did not exist in the population for :math:`\geq\ \mathrm{minsojourn}` generations
    :param minfreq: Exclude all mutations that did not read a minimum frequency of :math:`\geq\ \mathrm{minfreq}`

    :return: A pandas.DataFrame

    :rtype: pandas.DataFrame
    """
    return pandas.DataFrame.from_dict(get_qtrait_traj(pop.pop.get(),minsojourn,minfreq))

def esize_freq(singlepop pop):
    """
    Returns effect size vs frequency

    :param pop: A :class:`fwdpy.fwdpy.singlepop` simulated using :func:`fwdpy.qtrait.qtrait.evolve_qtrait` and/or :func:`fwdp.qtrait.qtrait.evolve_qtrait_more`

    :rtype: pandas.DataFrame
    """
    return pandas.DataFrame.from_dict(qtrait_esize_freq(pop.pop.get()))

def ew2010_effects(GSLrng rng, singlepop pop, double tau, double sigma):
    """
    :param pop: A :class:`fwdpy.fwdpy.singlepop` simulated using :func:`fwdpy.qtrait.qtrait.evolve_qtrait` and/or :func:`fwdp.qtrait.qtrait.evolve_qtrait_more`
    :param tau: The coupling of trait value to fitness effect of mutation
    :param sigma: The standard deviation for Gaussian noise applied to trait value.  Generates the :math:`\epsilon` term in E-W's paper
    """
    if tau < 0.:
        raise RuntimeError("tau cannot be < 0.")
    if sigma < 0.:
        raise RuntimeError("sigma cannot be < 0.")
    return ew2010_assign_effects(rng.thisptr,pop.pop.get(),tau,sigma)
    
def ew2010_traits(singlepop pop,list effects):
    """
    Implement model of Eyre-Walker 2010
    
    .. note:: The citation is www.pnas.org/cgi/doi/10.1073/pnas.0906182107.
       We implement the simple additive case here.
    """
    return ew2010_traits_cpp(pop.pop.get(),effects)
