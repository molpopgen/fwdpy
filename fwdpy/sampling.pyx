from cython.operator cimport dereference as deref
from fwdpy.fwdpp cimport sep_sample_t,sample_t,gsl_rng
import numpy as np
import pandas as pd

##Undocumented fxns are wrappers to enable run-time polymorphism within the Py environs.
##These fxns make calls to the C++ layer

#cdef sample_t ms_sample_single_deme(GSLrng rng, singlepop pop, int nsam, bint removeFixed) nogil:
#    return sample_single[singlepop_t](rng.thisptr.get(),deref(pop.pop.get()),nsam, int(removeFixed))

#cdef sep_sample_t ms_sample_single_deme_sep(GSLrng rng, singlepop pop, int nsam, bint removeFixed) nogil:
#    return sample_sep_single[singlepop_t](rng.thisptr.get(),deref(pop.pop.get()),nsam,removeFixed)

#cdef sep_sample_t ms_sample_metapop_sep(GSLrng rng, metapop pop, int nsam, bint removeFixed,int deme) nogil:
#    return sample_separate[metapop_t](rng.thisptr.get(),deref(pop.mpop.get()),deme,nsam,removeFixed)

#cdef vector[sep_sample_t] ms_sample_singlepop_mloc(GSLrng rng, singlepop_mloc pop, int nsam, bint removeFixed) nogil:
#    return sample_sep_single_mloc[multilocus_t](rng.thisptr.get(),deref(pop.pop.get()),nsam,removeFixed)

#cdef get_sh_single(const sample_t & ms_sample,
#                    singlepop pop,
#                    vector[double] * s,
#                    vector[double] * h,
#                    vector[double] * p,
#                    vector[double] * a,
#                    vector[uint16_t] * l):
#    get_sh(ms_sample,pop.pop.get(),s,h,p,a,l)

#cdef get_sh_metapop(const sample_t & ms_sample,
#                    metapop pop,
#                    vector[double] * s,
#                    vector[double] * h,
#                    vector[double] * p,
#                    vector[double] * a,
#                    vector[uint16_t] * l):
#    get_sh(ms_sample,pop.mpop.get(),s,h,p,a,l)

def ms_sample(GSLrng rng, poptype pop, int nsam, bint removeFixed = True):
    """
    Take a sample from a set of simulated populations.

    :param rng: a :class:`GSLrng`
    :param pops: An object inheriting from :class:`poptype`
    :param nsam: List of sample sizes (no. chromosomes) to sample.
    :param removeFixed: if True, only polymorphic sites are retained

    Note: nsam will likely be changed to a list soon, to accomodate multi-deme simulations

    Example:

    >>> import fwdpy,array
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes=array.array('I',[1000]*1000)
    >>> pop = fwdpy.evolve_regions(rng,3,1000,popsizes[0:],0.001,0,0.001,[fwdpy.Region(0,1,1)],[],[fwdpy.Region(0,1,1)])
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    """
    if isinstance(pop,singlepop):
        return sample_single[singlepop_t](rng.thisptr.get(),deref((<singlepop>pop).pop.get()),nsam, int(removeFixed))
    else:
        raise ValueError("ms_sample: unsupported type of popcontainer")

def get_samples(GSLrng rng, poptype pop, int nsam, bint removeFixed = True, deme = None):
    """
    Take a sample from a set of simulated populations.

    :param rng: a :class:`GSLrng`
    :param pop: An object inheriting from :class:`poptype`
    :param nsam: The sample size to take.
    :param removeFixed: if True, only polymorphic sites are retained
    :param deme: Optional.  If 'pop' is a :class:`metapop`, deme is required and represents the sub-population to sample.

    :return: A list. Element 0 is neutral mutations, and element 1 is selected mutations.  Within each list is a tuple of size 2.  The first element is the mutation position.  The second element is the genotype for each of the 'nsam' chromosomes.  Genotypes are coded as 0 = the ancestral state and 1 = the derived state.  For each site, each pair of genotypes constitutes a single diploid.  In other words, for nsam = 50, the data will represent the complete haplotypes of 25 diploids.

    :raise: IndexError if 'deme' is out of range and pop is a :class:`fwdpy.fwdpy.metapop`

    Please note that if you desire an odd 'nsam', you should input nsam+2 and randomly remove one haplotype to obtain your desired sample size.  This is due to an issue with how we are sampling chromosomes from the population.

    Example:

    >>> import fwdpy,array
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes=array.array('I',[1000]*1000)
    >>> pop = fwdpy.evolve_regions(rng,3,1000,popsizes[0:],0.001,0,0.001,[fwdpy.Region(0,1,1)],[],[fwdpy.Region(0,1,1)])
    >>> s = [fwdpy.get_samples(rng,i,10) for i in pop]
    """
    if isinstance(pop,singlepop):
        return sample_sep_single[singlepop_t](rng.thisptr.get(),deref((<singlepop>pop).pop.get()),nsam, int(removeFixed))
        #return ms_sample_single_deme_sep(rng,pop,nsam,removeFixed)
    elif isinstance(pop,metapop):
        if deme is None:
            raise RuntimeError("deme may not be set to None when sampling from a meta-population")
        if deme >= len(pop):
            raise RuntimeError("value for deme out of range. len(pop) = "+str(len(pop))+", but deme = "+str(deme))
        return sample_separate[metapop_t](rng.thisptr.get(),deref((<metapop>pop).mpop.get()),deme,nsam,removeFixed)
        #return ms_sample_metapop_sep(rng,pop,nsam,removeFixed,deme)
    elif isinstance(pop,singlepop_mloc):
        return sample_sep_single_mloc[multilocus_t](rng.thisptr.get(),deref((<singlepop_mloc>pop).pop.get()),nsam,removeFixed)
        #return ms_sample_singlepop_mloc(rng,pop,nsam,removeFixed)
    else:
        raise ValueError("ms_sample: unsupported type of popcontainer")

def get_sample_details( sample_t ms_sample, poptype pop ):
    """
    Get additional details for population samples

    :param ms_samples: A list returned by :func:`ms_sample`
    :param pops: A :class:`poptype`

    :return: A pandas.DataFrame containing the selection coefficient (s), dominance (h), populations frequency (p), and age (a) for each mutation.

    :rtype: pandas.DataFrame

    Example:

    >>> import fwdpy,array
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes=array.array('I',[1000]*1000)
    >>> pop = fwdpy.evolve_regions(rng,3,1000,popsizes[0:],0.001,0,0.001,[fwdpy.Region(0,1,1)],[],[fwdpy.Region(0,1,1)])
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    >>> details = [fwdpy.get_sample_details(i,j) for i,j in zip(s,pop)]
    """
    cdef vector[double] h
    cdef vector[double] s
    cdef vector[double] p
    cdef vector[double] a
    cdef vector[uint16_t] l
    if isinstance(pop,singlepop):
        get_sh(ms_sample,
               (<singlepop>pop).pop.get().mutations,
               (<singlepop>pop).pop.get().mcounts,
               (<singlepop>pop).pop.get().N,
               (<singlepop>pop).pop.get().generation,
               &s,&h,&p,&a,&l)
    elif isinstance(pop,singlepop_mloc):
        get_sh(ms_sample,
               (<singlepop_mloc>pop).pop.get().mutations,
               (<singlepop_mloc>pop).pop.get().mcounts,
               (<singlepop_mloc>pop).pop.get().N,
               (<singlepop_mloc>pop).pop.get().generation,
               &s,&h,&p,&a,&l)
    elif isinstance(pop,metapop):
        get_sh(ms_sample,
               (<metapop>pop).mpop.get().mutations,
               (<metapop>pop).mpop.get().mcounts,
               sum((<metapop>pop).mpop.get().mcounts.Ns),
               (<metapop>pop).mpop.get().generation,
               &s,&h,&p,&a,&l)
    return pandas.DataFrame({'s':s,'h':h,'p':p,'a':a,'label':l})


###### Functions for manipulating samples.
def nderived_site(tuple site):
    """
    Get the number of derived mutations at a site.

    :param site: A tuple.  See example

    .. note:: In general, it will be more convenient to call :func:`fwdpy.fwdpy.nderived` on a list of tuples.

    Example:

    >>> import fwdpy
    >>> #Create a site at position 0.1 with the
    >>> #genotypes as the second element. 0/1 = ancestral/derived
    >>> x = (0.1,'01100111')
    >>> fwdpy.nderived_site(x)
    5
    >>> #Process simulation results:
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes=array.array('I',[1000]*1000)
    >>> pop = fwdpy.evolve_regions(rng,3,1000,popsizes[0:],0.001,0,0.001,[fwdpy.Region(0,1,1)],[],[fwdpy.Region(0,1,1)])
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    >>> for i in s[0]: ndi = fwdpy.nderived_site(i)
    """
    return str(site[1]).count('1')

def nderived( list sample ):
    """
    Convenience wrapper around :func:`fwdpy.fwdpy.nderived`

    :param sample: a sample from a population.  For example, the return value of :func:`fwdpy.fwdpy.ms_sample` or :func:`fwdpy.fwdpy.get_samples`

    Example:

    >>> import fwdpy,array
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes=array.array('I',[1000]*1000)
    >>> pop = fwdpy.evolve_regions(rng,3,1000,popsizes[0:],0.001,0,0.001,[fwdpy.Region(0,1,1)],[],[fwdpy.Region(0,1,1)])
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    >>> nd = [fwdpy.nderived(i) for i in s]
    """
    return [nderived_site(i) for i in sample]

def getfreq(tuple site,bint derived = True):
    """
    Get mutation frequencies

    :param site: A tuple. See example
    :param derived:  If True, report derived allele frequency (DAF).  If False, return minor allele freqency (MAF)

    .. note:: Do **not** use this function to calculate :math:`\pi` (a.k.a. :math:`\\hat\\theta_\pi`, a.k.a. "sum of site heterozygosity").
       :math:`\pi` for a **sample** is not :math:`2\sum_ip_iq_i`. because the sample is *finite*.
       In general, it will be more convenient to call :func:`fwdpy.fwdpy.getfreqs` on a list of tuples.

    Example:

    >>> import fwdpy
    >>> #DAF = 7/10
    >>> #MAF = 3/10
    >>> x = (0.1,'0111111100')
    >>> round(fwdpy.getfreq(x,True),3)
    0.7
    >>> round(fwdpy.getfreq(x,False),3)
    0.3
    """
    o = nderived_site(site)
    dfreq = float(o)/float(len(site[1]))
    if derived is True:
        return dfreq
    return min(dfreq,1.0-dfreq)

def getfreqs(list sample,bint derived = True):
    """
    Convenience wrapper around :func:`fwdpy.fwdpy.getfreq`

    :param sample: a sample from a population.  For example, the return value of :func:`fwdpy.fwdpy.ms_sample` or :func:`fwdpy.fwdpy.get_samples`
    :param derived: If True, report derived allele frequency (DAF).  If False, return minor allele freqency (MAF).

    .. note:: Do **not** use this function to calculate :math:`\pi` (a.k.a. :math:`\\hat\\theta_\\pi`, a.k.a. "sum of site heterozygosity").
       :math:`\pi` for a **sample** is not :math:`2\sum_ip_iq_i`. because the sample is *finite*.

    Example:

    >>> import fwdpy,array
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes=array.array('I',[1000]*1000)
    >>> pop = fwdpy.evolve_regions(rng,3,1000,popsizes[0:],0.001,0,0.001,[fwdpy.Region(0,1,1)],[],[fwdpy.Region(0,1,1)])
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    >>> freqs = [fwdpy.getfreqs(i) for i in s]
    """
    return [getfreq(i,derived) for i in sample]

def freqfilter( list sample,
                float minfreq,
                bint derived = True ):
    """
    Remove low-frequency variants from a sample.

    :param sample: a sample from a population.  For example, the return value of :func:`fwdpy.fwdpy.ms_sample` or :func:`fwdpy.fwdpy.get_samples`
    :param minfreq: Remove all sites with frequency < minfreq
    :param derived: if True, filter on derived allele frequency.  If False, filter on minor allele frequency.

    Example:

    >>> import fwdpy,array
    >>> rng = fwdpy.GSLrng(100)
    >>> popsizes=array.array('I',[1000]*1000)
    >>> pop = fwdpy.evolve_regions(rng,3,1000,popsizes[0:],0.001,0,0.001,[fwdpy.Region(0,1,1)],[],[fwdpy.Region(0,1,1)])
    >>> s = [fwdpy.ms_sample(rng,i,10) for i in pop]
    >>> s2 = [fwdpy.freqfilter(i,0.2) for i in s]
    """
    rv=list()
    for i in sample:
        if type(i) is not tuple:
            raise RuntimeError("values is sample must be tuples, not"+str(type(i)))
        if getfreq(i,derived)>=minfreq:
            rv.append(i)
    return rv

def hapmatrix(p,const vector[size_t] & diploids, deme = None):
    """
    Obtain compact representations of a "haplotype matrix" from a population object

    :param p: A population object or vector of such objects
    :param diploids: Indexes of individuals to include in matrix
    :param deme: Deme index (only needed if p represents a meta-population)

    :return: A dict that will allow construction of a 0/1 matrix represting ancestral and derived genotypes,
    respectively, corresponding the the individuals in 'diploids'.

    .. note:: Examples will be provided elsewhere.
    """
    if isinstance(p,singlepop):
        return make_haplotype_matrix((<singlepop>p).pop.get(),diploids)
    elif isinstance(p,singlepop_mloc):
        return make_haplotype_matrix((<singlepop_mloc>p).pop.get(),diploids)
    elif isinstance(p,popvec):
        return [make_haplotype_matrix((<singlepop>i).pop.get(),diploids) for i in <popvec>p]
    elif isinstance(p,popvec_mloc):
        return [make_haplotype_matrix((<singlepop_mloc>i).pop.get(),diploids) for i in <popvec_mloc>p]
    elif isinstance(p,metapop):
        if deme is None:
            raise RuntimeError("deme cannot be None")
        return make_haplotype_matrix((<metapop>p).mpop.get(),diploids,deme)
    elif isinstance(p,mpopvec):
        if deme is None:
            raise RuntimeError("deme cannot be None")
        return [make_haplotype_matrix((<metapop>i).mpop.get(),diploids,deme) for i in <mpopvec>p]
    else:
        raise RuntimeError("object type not understood")

cdef row_col_indexes(list indexes,
                     size_t ncol):
    rv=[]
    for i in indexes:
        t=float(i+1)/float(ncol)
        row=np.floor(t)
        col=(t-row)*float(ncol)-1.0
        rv.append((int(row),int(col)))
    return rv

def genomatrix(dict m):
    """
    Generate a "genotype matrix" from a haplotype matrix

    :param m: The return value from :func:`fwdpy.fwdpy.hapmatrix`

    :return: A dict containing the indexes of Aa and aa genotypes
    for neutral (n) and selected (s) markers in a 1d array of length m['nrow']*m['ncol_n']/2

    ..note :: The return value is most useful when used to create a 2d numpy array representing 
    genotypes.
    """

    #Step 1: get the row/column index of each mutation in the haplotype matrix
    n=row_col_indexes(m['n'],m['ncol_n'])
    s=row_col_indexes(m['s'],m['ncol_s'])

    #Step 2: For each pair of rows, collect Aa and aa genotypes
    nAa=[]
    naa=[]
    sAa=[]
    saa=[]
    newrows=m['nrow']/2
    newrow=0
    for row in range(0,m['nrow'],2):
        #Get data for this and the next row in the hap matrix
        x=[i[1] for i in n if (i[0]==row or i[0]==(row+1))]
        xu=set(x)
        for i in xu:
            ci=x.count(i)
            if ci == 1:
                nAa.append(newrow*m['ncol_n']+i)
            elif ci == 2:
                naa.append(newrow*m['ncol_n']+i)

        #do same for selected variants
        x=[i[1] for i in s if (i[0]==row or i[0]==(row+1))]
        xu=set(x)
        for i in xu:
            ci=x.count(i)
            if ci == 1:
                sAa.append(newrow*m['ncol_s']+i)
            elif ci == 2:
                saa.append(newrow*m['ncol_s']+i)

        newrow+=1

    return {'nAa':nAa,'naa':naa,'sAa':sAa,'saa':saa}
