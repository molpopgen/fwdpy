
Extending fwdpy with Cython
===========================

.. code:: python

    %load_ext Cython

.. code:: python

    import fwdpy as fp
    import numpy as np

.. code:: python

    fwdpy_includes = fp.get_includes()
    fwdpp_includes = fp.get_fwdpp_includes()

.. code:: python

    %%cython --cplus --compile-args=-std=c++11 -I $fwdpy_includes -I $fwdpp_includes -l sequence -l gsl -l gslcblas
    from fwdpy.fwdpy cimport *
    cdef vector[unsigned] sfs_cpp(const singlepop_t * pop):
        cdef vector[unsigned] rv
        rv.resize(2*pop.N,0)
        cdef size_t i = 0
        for i in range(pop.mcounts.size()):
            if pop.mcounts[i]>0:
                rv[pop.mcounts[i]-1]+=1
        return rv
    
    def sfs(Spop pop):
        return sfs_cpp(pop.pop.get())

.. code:: python

    N=1000
    theta=100.
    nlist=np.array([N]*(10*N),dtype=np.uint32)
    rng = fp.GSLrng(135123)
    nregions=[fp.Region(0,1,1)]
    sregions=[]
    recregions=nregions
    pops = fp.evolve_regions(rng,1,N,nlist,theta/(4.*float(N)),0.,theta/(4.*float(N)),nregions,sregions,recregions)

.. code:: python

    sfs_pop=sfs(pops[0])
    print(sfs_pop[0:10])


.. parsed-literal::

    [106, 49, 40, 37, 14, 16, 11, 16, 2, 4]


