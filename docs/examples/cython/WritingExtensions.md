
# Extending fwdpy with Cython

Cython makes it easy to [share definitions](http://docs.cython.org/en/latest/src/userguide/sharing_declarations.html) between packages, making it easy to use fwdpy's types to write custom code.  Further, as fwdpy depends on and installs [fwdpp](https://molpopgen.github.io/fwdpp), you get access to many of that library's features.  Even better, you can write your extensions and ignore a lot of gory details regarding compiling and linking--Cython handles that for you!

This document serves as a rapid-fire tutorial both to the C++ types that underly fwdpy and how to use Cython to write your own extensions.

## Plugins vs. packages

For the most part, Cython code is written in files with the extension .pyx.  If you can write all of your extension in Cython, you may simply compile and import your module into a script using [pyximport](http://cython.readthedocs.io/en/latest/src/reference/compilation.html).  For user's familiar with [Rcpp](http://rcpp.org), think of pyimport as the analog to their `sourceCpp` function.

When using pyximport, the function `fwdpy.make_pyxbld` will help you out a lot--please see its documentation in the reference manual.

If you start writing a lot of extensions or your extensions require C++11 features that Cython cannot handle, then you may want to consider writing a full-blown package for your extensions. There are lots of examples online, from the Cython documentation to how the [fwdpy source code](http://github.com/molpopgen/fwdpy) is organized.

## Linux vs OS X

Due to issues with compiler support on OS X, Linux is the intended platform for using fwdpy. It is possible to install the package if you use GCC, which you can install via Anaconda.  

When compiling extensions, Python's distutils attempts to force use of the same compiler used to build Python.  On OS X, that means clang, but fwdpy requires GCC on OS X.  Thus, you need to force the use of GCC via the CC/CXX environment variables.

## An important note

Many of the example functions below actually end up replicating things that are already doable in fwdpy.  In other words, you don't need any of the stuff below to do what is below.  These are examples for the point of documenting the C++/Cython API that you have access to.

## Cython 'magic' lines

Every Cython code block in this document begins with a line starting "%%cython".  That's another 'magic' command for the Jupyter notebooks.  It contains info needed to compile each code block.  You can basically ignore that.


```python
#This is a 'magic' command allowing us to 
#use Cython in a Jupyter notebook, which is
#what we use to write this document.
%load_ext Cython
```


```python
import fwdpy as fp
import numpy as np
```

# Finding the headers
fwdpy provides functions that reveal the locations of both the fwdpy C++ header files and the fwdpp C++ header files that are installed along with fwdpy.  You need to know these locations!


```python
fwdpy_includes = fp.get_includes()
fwdpp_includes = fp.get_fwdpp_includes()
```

# Example 1: the site-frequency spectrum of all mutations

The first function that we will write will calculate the site-frequency-spectrum (SFS) of the entire population.  We impose the following constraints to keep things simple:

* We will only process single-deme objects (type fwdpy.Spop).

On to our code for the SFS:


```python
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
```


```python
N=1000
theta=100.
nlist=np.array([N]*(10*N),dtype=np.uint32)
rng = fp.GSLrng(135123)
nregions=[fp.Region(0,1,1)]
sregions=[]
recregions=nregions
#Simulate 10 populations
pops = fp.evolve_regions(rng,10,N,nlist,theta/(4.*float(N)),0.,theta/(4.*float(N)),nregions,sregions,recregions)
```


```python
sfs_pop=sfs(pops[0])
print(sfs_pop[0:10])
print(type(sfs_pop))
```

    [106  49  40  37  14  16  11  16   2   4]
    <type 'numpy.ndarray'>


Get the mean SFS for our 10 replicates:


```python
mean_sfs = np.sum([sfs(i) for i in pops],axis=0)/10.
mean_sfs
```




    array([ 109. ,   50.1,   36.5, ...,    0. ,    0. ,    0. ])



## Pythonic or not?

The `sfs_cpp` function takes a const pointer for an argument.  If we relax that constraint, we can write some of the details in a more relaxed, Pythonic manner:


```python
%%cython --cplus --compile-args=-std=c++11 -I $fwdpy_includes -I $fwdpp_includes -l sequence -l gsl -l gslcblas
from fwdpy.fwdpy cimport *
import numpy as np
#A non-const pointer now:
cdef vector[unsigned] sfs_cpp_pythonic(singlepop_t * pop):
    cdef vector[unsigned] rv
    rv.resize(2*pop.N,0)
    cdef size_t i = 0
    #When operating in a non-const
    #context, you can use 
    #Python-like syntax
    #to iterate over C++
    #containers:
    for i in pop.mcounts:
        if i>0:
            rv[i-1]+=1
    return rv

def sfs_pythonic(Spop pop):
    """
    This is another Python function that will return the 
    SFS for a fwdpy.Spop object.
    
    :param pop: A :class:`fwdpy.fwdpy.Spop`
    
    :return: The site-frequency spectrum for pop
    
    :rtype: numpy.array with dtype numpy.uint32
    """
    return np.array(sfs_cpp_pythonic(pop.pop.get()),dtype=np.uint32)
```

We get the same results:


```python
mean_sfs = np.sum([sfs_pythonic(i) for i in pops],axis=0)/10.
mean_sfs
```




    array([ 109. ,   50.1,   36.5, ...,    0. ,    0. ,    0. ])



Why would you use the more complex first method?  From a C++ purist's perspective, the latter function protoype (with the non-const pointer argument) is annoying.  While the function does not modify the input value, but you cannot know that without reading its implementation in detail.  Personally, I like having the function fail to compile if I accidentally try to modify a constant object.

## Getting the SFS from fwdpy

Remember, the above code replicates existing fwdpy features.  To get the SFS, use "views" of the mutations in your simulation:


```python
mean_sfs_views = np.array([0.]*2*N)
for v in fp.view_mutations(pops):
    for m in v:
        mean_sfs_views[m['n']-1]+=1
mean_sfs_views /= 10.
mean_sfs_views
```




    array([ 109. ,   50.1,   36.5, ...,    0. ,    0. ,    0. ])



# Separating the neutral and selected SFS


```python
#Now simulated selected variants
sregions=[fp.GammaS(0,1,0.9,-0.043,0.23,1),
         fp.ExpS(0,1,0.1,0.01,1)]
theta_selected = 0.1*theta
#Re-simulate 10 populations
pops = fp.evolve_regions(rng,10,N,nlist,theta/(4.*float(N)),theta_selected/(4.*float(N)),theta/(4.*float(N)),nregions,sregions,recregions)
```


```python
%%cython --cplus --compile-args=-std=c++11 -I $fwdpy_includes -I $fwdpp_includes -l sequence -l gsl -l gslcblas
#Import all Cython symbols defined
#in fwdpy's main module
from fwdpy.fwdpy cimport *
from libcpp.utility cimport pair
import numpy as np

ctypedef vector[unsigned] vu

cdef pair[vu,vu] sfs_sep_cpp(const singlepop_t * pop):
    cdef pair[vu,vu] rv
    rv.first.resize(2*pop.N,0)
    rv.second.resize(2*pop.N,0)
    cdef size_t i = 0
    for i in range(pop.mcounts.size()):
        if pop.mcounts[i]>0:
            #Populations store their mutations
            #in a vector. A mutation
            #contains a boolean recording its
            #"neutrality":
            if pop.mutations[i].neutral is True:
                #The first element will be the
                #neutral SFS
                rv.first[pop.mcounts[i]-1]+=1
            else:
                #The second will be the selected
                #SFS
                rv.second[pop.mcounts[i]-1]+=1
    #Return the SFS to Python.
    #Cython auto-converts the
    #pair of vectors to a 
    #tuple of lists
    return rv

def sfs_sep(Spop pop):
    """
    This is the Python function that will return the 
    SFS for a fwdpy.Spop object.  The sfs will be 
    separate for neutral variants
    
    :param pop: A :class:`fwdpy.fwdpy.Spop`
    
    :return: The site-frequency spectrum for pop, separating
    neutral and selected variants
    
    :rtype: tuple of numpy.array with dtype numpy.uint32
    """
    return np.array(sfs_sep_cpp(pop.pop.get()),dtype=np.uint32)
```

Let's apply our new function and get the mean normalized SFS for neutral and selected variants.


```python
pop_sfs_sep = [sfs_sep(i) for i in pops]
#Note that we need to cast one array from uint32 to float,
#so that numpy promotes the calculation to floating-point.
mean_norm_sfs_neut = np.sum([i[0].astype(np.float)/np.sum(i[0]) for i in pop_sfs_sep],axis=0) / float(len(pops))
mean_norm_sfs_sel = np.sum([i[1].astype(np.float)/np.sum(i[1]) for i in pop_sfs_sep],axis=0) / float(len(pops))
print(mean_norm_sfs_neut)
print(mean_norm_sfs_sel)
```

    [ 0.18459531  0.07918477  0.06201459 ...,  0.          0.          0.        ]
    [ 0.22176807  0.08613146  0.07475673 ...,  0.          0.          0.        ]

