Modeling regional variation in mutation and recombination
======================================================================

Several of the simulation routines allow the details of the mutation and recombination models to vary along a "sequence" or "region".  A user is able to specify the details of such variation by passing _lists_ to package functions.  For example, you are able to:

* Vary the neutral mutation rate along a sequence.
* Vary the distribution of selection coefficients (and the dominance associated with selected mutations) along a sequence.
* Vary the recombination rate along a sequence.

The implementation of such variation along a region is *discrete*.  A region is specified by a beginning, and end, and a weight, plus any additional data required to specify selection coefficients, dominance, etc.

Background
--------------------------------------------------
The models are parameterized through Python's "new-style" class system.

Mutation rates, recombination rates, and a weighting system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A simulation will typically have a mutation rate, :math:`\\mu`, which represents the mean of a Poisson number of mutations per gamete per generation), and a recombination rate, :math:`r`, which again is the mean of Poisson number of crossover events (per diploid, per generation).  These parameters are the _total_ rates across an entire simulated region.  Variation in these parameters along the region are affected by a set of positions coupled with "weights", which the user specifies using S4 classes.

The base class: :class:`fwdpy.fwdpy.Region`

A :class:`fwdpy.fwdpy.Region` is an Python class with the following members:

* :math:`b`, which is the beginning/start of the region. The type is "float". 
* :math:`e`, which is the end/stop of the region. The type is "float".
* :math:`w`, which is a weighting factor associated with the region. The type is "float".

The members are used to inform the C++ code about the relative abundance of new mutations or recombination events will occur in what region.  Briefly, the number of events that occur in region :math:`i` are proportional to :math:`w_i/\sum_i w`, *i.e*, the weight assigned to region :math:`i` divided by the sum of weights assigned to all regions.  The weights for mutation events and for recombination events are considered separately.  Thus, in order to model a correlation between mutational processes and recombination, it is up to the user to generate regions whose weights are correlated.

fwdpy allows the :math:`w` slot to be interpreted in one of two ways:

* It is *not*  affected by the length of region.  Interally, the weight assigned is simply :math:`w`. 
* It is affected by the length of a region :math:`(e - b)`.

These two options are determined by arguments to class constructors, which we will see in examples below.  The latter is the default.

These two approaches allow for considerable modeling flexibility.  For example, the latter approach allows :math:`w` to be interpreted as a "per base-pair" rate.  Imagine that you wanted to simulate variation in recombination along discrete 100 kilobase chunks, and the rate of crossing-over *per base pair* increases in each chunk, and includes an initial chunk with no recombination:

1. start=1,stop= :math:`10^5`, :math:`r_{bp}=0`
2. start= :math:`10^5`,stop= :math:`2 \times 10^5`, :math:`r_{bp}=10^{-8}`
3. start= :math:`2 \times 10^5`,stop= :math:`3 \times 10^5`, :math:`r_{bp}=10^{-7}`  


This model boils down to the relative number of crossing overs per region occuring in the ratio :math:`0 : 10^{-8} : 10^{-7}`.  This is easily represented using fwdpy's classes:

"``>>>``"
   >>> import fwdpy 
   >>> recRegions = [fwdpy.Region(1,1e5,0),fwdpy.Region(1e5,2e5,1e-8),fwdpy.Region(2e5,3e5,1e-7)]

For this hypothetical example, the region lengths are all even, and thus an equivalent specification would be this:

"``>>>``"
   >>> import fwdpy 
   >>> recRegions = [fwdpy.Region(1,1e5,0,False),fwdpy.Region(1e5,2e5,1e-8,False),fwdpy.Region(2e5,3e5,1e-7,False)]
   
Specific examples
-------------------

Mutations not affecting fitness ("neutral" mutations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You specify regions where neutral mutations arise via the class :class:`fwdpy.fwdpy.Region`.  A region has a beginning, end, and a weight Thus, the following list would specify that 100% of neutral mutations occur on the continuous interval [0,1):

.. code-block:: python
		
		neutralRegions = [fwdpy.Region(0,1,1)]

The beginning and end positions can be whatever you like:

.. code-block:: python

		#With a weight of 1, we're just rescaling the position here.
		neutralRegions = [fwdpy.Region(0,100,1)]

To specify variation in the netural mutation process along a sequence, combine multiple regions in your list:

"``>>>``"

>>> #If coupled=False for the second region, the effect would be that region2's mutation rate per base pair is 10x less than region 1!!
>>> neutralRegions = [fwdpy.Region(beg=0,end=1,weight=1),fwdpy.Region(beg=2,end=12,weight=1,coupled=True)]


Internally, the total "mutational weight" of the first region will be a function of its length, which is 1(1-0)=1.  The second region's total weight will be 1*(12-2)=10, and it will have 10xas many new mutations arising as the first region.

Mutations affecting fitness
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Type types of mutations affecting fitness that we consider will have two parameters associated with them:

* :math:`s`, the selection coefficient
* :math:`h`, the effect of the mutation in a heterozygote (a.k.a. the "dominance" of the mutation).

In a simulation, we may place a distribution on either :math:`s` itself or on the scaled selection parameter :math:`\alpha = 2Ns`.  These two methods are represented by the class :class:`fwdpy.fwdpy.Sregion`.  These classes contain/extend the :class:`Region` class described above, and thus inherit their members.  :class:`Sregion` adds :math:`h`, which is the dominance of a mutation, and then classes extending :class:`Sregion` add details about the distribution of fitness effects.  These classes are:

* :class:`fwdpy.fwdpy.ConstantS`
* :class:`fwdpy.fwdpy.UniformS`
* :class:`fwdpy.fwdpy.GammaS`
* :class:`fwdpy.fwdpy.GaussianS`
  
Crossover rate variation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Just like neutral mutations, intervals with different crossover rates are specified by different :class:`Region` objects.  Let's set up the following concrete example:

* A region where crossovers occur between positions [0,1)
* Positions [0,0.45) and [0.55,1) have uniform recombintion rates
* Positions [0.45,0.55) are a recombination hotspot with 100x the background intensity (per "base pair").

The above model can be represented as:

"``>>>``"

>>> import fwdpy
>>> #recrate[2] is the hotspot:
>>> recrates = [fwdpy.Region(0,0.45,1),fwdpy.Region(0.55,1,1),fwdpy.Region(0.45,0.55,100)]

Internally, this is what will happen to the above input:

* The total weight on the first region will be :math:`w = w \times (e-b) = 1\times(0.45-0) = 0.45`
* The weight on the second region will be :math:`1\times(1-0.55) = 0.45`
* The weight on the hotspot will be :math:`100\times(0.55-0.45) = 10`

This gives us what we want: the hotspot is 100x hotter "per base", and is 10% of the total region in length.  We therefore expect 10x as many crossovers in that region as in the flanking regions.

How to set up a model
---------------------------------

When setting up a model, it is important that you think in terms of conditional probabilities.  In other words, if the total rate to neutral variants is :math:`\mu_n`, then the weights passed along to a function have the interpretations "Given that a neutral mutation occurs, the probability that it occurs in a certain interval is :math:`x`, where :math:`x` is determined by the relative weight assigned to an interval.

The 'weights' that you assign are *relative* and need not sum to 1.  Each weight must be :math:`\geq 0`, though.

Example
~~~~~~~~~~~

"``>>>``"

>>> import fwdpy
>>> import numpy as np
>>> rng = fwdpy.GSLrng(100)
>>> ##Some basic parameters
>>> N=1000
>>> theta=100.0
>>> rho=100.0
>>> ##All neutral muts are [0,1)
>>> nregions = [ fwdpy.Region(0,1,1) ]
>>> #Selected mutations.  All are additive, to keep this example simple.
>>> ##Strongly deleterious mutations to the "left"
>>> ##Weaker mutations (2Ns = 10 on average) to the "right"
>>> ## 1% of selected mutations will be positively selected
>>> ## and uniform throughout the region.  The distribution
>>> ## of s will be exponential with mean 1e-3
>>> smodels = [fwdpy.ConstantS(-1,0,0.99/2,-0.1),fwdpy.ExpS(1,2,0.99/2,-10),fwdpy.ExpS(-1,2,0.01,0.001)]
>>> ##Recombination models--10x hotspot in the middl
>>> rregions = [fwdpy.Region(-1,1,1),fwdpy.Region(0.45,0.55,10)]
>>> #set up list of population sizes,
>>> #which are NumPy arrays of ints
>>> popsizes = np.array([N],dtype=np.uint32) 
>>> popsizes = np.tile(popsizes,10*N)
>>> pops = fwdpy.evolve_regions(rng,1,N,popsizes[0:],theta/(4*N),0.1*theta/(4*N),rho/(4*N),nregions,smodels,rregions)
>>> #Take a sample of size n = 10 from the population via list comprehension
>>> popSample = [fwdpy.ms_sample(rng,i,10) for i in pops]

A simple example of "exons"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's consider the following model of a "two-exon gene":

* Exon 1 spans positions 0.2 to 0.4.
* Exon 2 spans positions 0.6 to 1.
* There will be selection against amino acid replacements.

For the sake of simplicity, we will ignore:

* mutations in the introns and UTR regions

We can define the following region lengths:

* :math:`l_1 = 0.2` is the length of Exon 1.
* :math:`l_2 = 0.4` is the length of Exon 2.

Thus the relative weights assigned to Exon 1 and 2 must satisfy :math:`w_2/w_1 = l_2/l_1 = 2`.

Further, within an exon, :math:`\approx 3/4` of new mutations will be amino acid replacements.  Thus, 3/4 of the total mutational weight in a region will be on selected mutations.

Our model looks like this:

"``>>>``"

		>>> nregions = [fwdpy.Region(0.2,0.4,0.25),fwdpy.Region(0.6,1,0.5)]
		>>> sregions = [fwdpy.ConstantS(0.2,0.4,0.75,-0.01),fwdpy.ConstantS(0.6,1,1.5,-0.01) ]


If we pass *the same neutral and selected mutation rates to evolve.regions*, then the above model satisfies:

* The total number of mutations occurring in Exon 2 is 2x the number occuring in Exon 1.
* Within an exon, 3/4 of all new mutations are deleterious.
