
Extending fwdpy with Cython
===========================

Cython makes it easy to `share
definitions <http://docs.cython.org/en/latest/src/userguide/sharing_declarations.html>`__
between packages, making it easy to use fwdpy's types to write custom
code. Further, as fwdpy depends on and installs
`fwdpp <https://molpopgen.github.io/fwdpp>`__, you get access to many of
that library's features. Even better, you can write your extensions and
ignore a lot of gory details regarding compiling and linking--Cython
handles that for you!

This document serves as a rapid-fire tutorial both to the C++ types that
underly fwdpy and how to use Cython to write your own extensions.

Plugins vs. packages
--------------------

For the most part, Cython code is written in files with the extension
.pyx. If you can write all of your extension in Cython, you may simply
compile and import your module into a script using
`pyximport <http://cython.readthedocs.io/en/latest/src/reference/compilation.html>`__.
For user's familiar with `Rcpp <http://rcpp.org>`__, think of pyimport
as the analog to their ``sourceCpp`` function.

When using pyxbld, the function ``fwdpy.make_pyxbld`` will help you out
a lot--please see its documentation in the reference manual.

If you start writing a lot of extensions or your extensions require
C++11 features that Cython cannot handle, then you may want to consider
writing a full-blown package for your extensions. There are lots of
examples online, from the Cython documentation to how the `fwdpy source
code <http://github.com/molpopgen/fwdpy>`__ is organized.

Linux vs OS X
-------------

An important note
-----------------

Many of the example functions below actually end up replicating things
that are already doable in fwdpy. In other words, you don't need any of
the stuff below to do what is below. These are examples for the point of
documenting the C++/Cython API that you have access to.

.. code:: python

    #This is a 'magic' command allowing us to 
    #use Cython in a Jupyter notebook, which is
    #what we use to write this document.
    %load_ext Cython

.. code:: python

    import fwdpy as fp
    import numpy as np

Finding the headers
===================

fwdpy provides functions that reveal the locations of both the fwdpy C++
header files and the fwdpp C++ header files that are installed along
with fwdpy. You need to know these locations!

.. code:: python

    fwdpy_includes = fp.get_includes()
    fwdpp_includes = fp.get_fwdpp_includes()

Example 1: the site-frequency spectrum of all mutations
=======================================================

The first function that we will write will calculate the
site-frequency-spectrum (SFS) of the entire population. We impose the
following constraints to keep things simple:

-  We will only process single-deme objects (type fwdpy.Spop).

Cython 'magic' lines
--------------------

Every Cython code block in this document begins with a line starting
"%%cython". That's another 'magic' command for the Jupyter notebooks. It
contains info needed to compile each code block. You can basically
ignore that.

On to our code for the SFS

.. code:: python

    %%cython --cplus --compile-args=-std=c++11 -I $fwdpy_includes -I $fwdpp_includes -l sequence -l gsl -l gslcblas
    #Import all Cython symbols defined
    #in fwdpy's main module
    from fwdpy.fwdpy cimport *
    import numpy as np
    #Now, we define a C++ function that:
    #1. Takes the C++ representation as an argument
    #2. Returns a C++ vector of unsigned integers
    cdef vector[unsigned] sfs_cpp(const singlepop_t * pop):
        #declare our return value.
        #This is a standard C++ vector.
        #The C++ vector is imported as a 
        #side-effect of cimporting fwdpp's
        #Cython API
        cdef vector[unsigned] rv
        #For a population of N diploids,
        #there are N bins in the SFS 
        #(including fixations, which
        #we don't deal with here).
        #So we initialize the return
        #value to 2N zeroes
        rv.resize(2*pop.N,0)
        
        #i is a dummy variable
        cdef size_t i = 0
        #A population contains a 
        #vector[unsigned] that represents
        #the count (no. occurrences) of
        #every mutation.  Warning: it also
        #conatains mutations with a count of
        #0 (zero) because fwdpp internally
        #puts new variants in those spaces...
        for i in range(pop.mcounts.size()):
            #...so we check that
            #a mutation's count
            #is nonzero...
            if pop.mcounts[i]>0:
                #...and increment our return value
                #accordingly.
                rv[pop.mcounts[i]-1]+=1
        #Return the SFS to Python
        return rv
    
    def sfs(Spop pop):
        """
        This is the Python function that will return the 
        SFS for a fwdpy.Spop object.
        
        Note that we can specify the argument type in the
        "def" line.  
        
        This docstring can be processed by Sphinx, and so
        we use Sphinx grammar for documenting the params,
        and we make sure to provide a link to the documentation
        of the parameter's expected type:
        
        :param pop: A :class:`fwdpy.fwdpy.Spop`
        
        :return: The site-frequency spectrum for pop
        
        :rtype: numpy.array with dtype numpy.uint32
        """
        #Here, we call our Cython function.
        #The fwdpy.Spop type contains a
        #std::unique_ptr[singlepop_t] object
        #called "pop".  So, we send the raw pointer
        #to our Cython function:
        return np.array(sfs_cpp(pop.pop.get()),dtype=np.uint32)

.. code:: python

    N=1000
    theta=100.
    nlist=np.array([N]*(10*N),dtype=np.uint32)
    rng = fp.GSLrng(135123)
    nregions=[fp.Region(0,1,1)]
    sregions=[]
    recregions=nregions
    pops = fp.evolve_regions(rng,10,N,nlist,theta/(4.*float(N)),0.,theta/(4.*float(N)),nregions,sregions,recregions)

.. code:: python

    sfs_pop=sfs(pops[0])
    print(sfs_pop[0:10])


.. parsed-literal::

    [106  49  40  37  14  16  11  16   2   4]


.. code:: python

    print(type(sfs_pop))


.. parsed-literal::

    <type 'numpy.ndarray'>


