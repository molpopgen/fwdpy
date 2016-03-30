Processing SLiM input files
=================================================

**fwdpy** has the ability to read in the same input files as Phillip Messer's SLiM_ simulation software, which is accomplished by the function :func:`fwdpy.fwdpy.readslim`.

This function does **not** read in all of the information in an input file.  Instead it reads in the following sections:

* #MUTATION RATE
* #MUTATION TYPES
* #GENOMIC ELEMENT TYPES
* #CHROMOSOME ORGANIZATION
* #RECOMBINATION RATE

These blocks are parsed and converted into data types that may serve as the input parameters for functions like :func:`fwdpy.fwdpy.evolve_regions`.

Differences between SLiM and fwdpy
------------------------------------------------

There are several important differences between SLiM_ and this package.

Scaling of fitnesses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For diploid genotypes :math:`AA`, :math:`Aa`, and :math:`aa`, SLiM_ calculates fitness as :math:`1`, :math:`1+sh`, and :math:`1+s`, respectively, where :math:`s` is the selection coefficient and :math:`h` the dominance associated with the :math:`a` allele.  In **fwdpy**, the fitnesses of the three genotypes are :math:`1`, :math:`1+sh`, and :math:`1+2s`, respectively.   To account for this difference, :func:`fwdpy.fwdpy.readslim` halves the selection coefficient (or the mean :math:`s` for distributions on selection coefficients) and doubles the dominance found in the #MUTATION TYPES block of a SLiM_ input file.

**IMPORTANT:** SLiM_ input files should be generated according to that software's manual!!!  In other words, co-dominance should be coded as :math:`h=0.5` and you should allow :func:`fwdpy.fwdpy.readslim` to convert it to one automatically!  Further, the fitnesses of the three genotypes should be treated as  :math:`1`, :math:`1+sh`, and :math:`1+s`, respectively. The intent here is to allow pre-existing SLiM_ input files to be used in **fwdpy** with no additional modification.

Discrete vs continuous positions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SLiM_'s mutation model involves *discrete* positions along the genome.  Therefore, the output from that software contains mutations with postitions whose values are integers.  Likewise, recombination events occur between integer-valued positions in SLiM_.

Internally, this package uses fwdpp_ to model mutation and recombination positions as *continuous* random variables along genomic intervals in a fashion similar to the coalescent simulation ms_.  (Techincally, ms_ treats mutation as continuous and recombination breakpoints as discrete, but let's ignore that point for now.)

An additional difference is that SLiM_ is paramterized in terms of "per base pair" mutation and recombination rates and this package uses "per region" rates.

:func:`fwdpy.fwdpy.evolve_regions` adjusts the SLiM_ input file to account for these differences.  Let's consider the following setup for recombination rates:

     | #RECOMBINATION RATE     
     | 10000 1e-7
     | 20000 2e-7
     | 30000 3e-7

The block above means that genomic positions 1 through 10,000 have a recombination rate of :math:`10^{-7}` *per base pair* (bp).  Positions 10,001 through 20,000 have a recombination rate of :math:`2 \times 10^{-7}/\mathrm{bp}`, and positions 20,001 through 30,000 have a recombination rate of :math:`3 \times 10^{-7}/\mathrm{bp}`.  The total/cumulative recombination rate is :math:`10^4\times(10^{-7} + 2\times 10^{-7} + 3 \times 10^{-7}) = 0.006`
The function :func:`fwdpy.fwdpy.readslim` will  convert the above block to three half-open, continuous intervals:

     | :math:`[0,10000)`, with *total* recombination rate :math:`10^{-7}\times(10^4-0) = 10^{-3}`
     | :math:`[10000,20000)` with *total* recombination rate :math:`2 \times 10^{-7}\times(2\times 10^4-10^4) = 2\times10^{-3}`
     | :math:`[20000,30000)` with *total* recombination rate :math:`3 \times 10^{-7}\times(3\times 10^4-2\times 10^4) = 3\times10^{-3}`

You can see that the total recombination rate is the same (:math:`0.006`).  The difference from SLiM_ is that crossover positions will be drawn from continuous uniform distributions.

The intervals where mutations occur are treated similarly.  Further, for a given genomic element, the mutation rate to a particular mutation type is the mutation rate specified in the #MUTATION RATE block times the weight assigned to that mutation type divided by the sum of all weights for all mutation types occuring in that genomic element.

A word of caution
-----------------------------

SLiM_ allows different mutations to have the same position label.  This software does not.  This difference may be relevant when comparing simulation output to theoretical predictions, which typically come from the infinitely-many sites assumption that this package models.  For a region of length :math:`L`, as :math:`\mu_{bp}` (the mutation rate per base pair) the output of SLiM_ and **fwdpy** should be comparable for large :math:`L` and :math:`\mu_{bp} \to 0`.  When mutation rates are large relative to :math:`L`, SLiM_'s mutation scheme becomes a quasi-finite-sites model.

Limitations
-------------------------

* Gene conversion is not currently supported
     
.. _SLiM: http://messerlab.org/software/
.. _fwdpp: http://github.com/molpopgen/fwdpp
.. _ms: http://home.uchicago.edu/~rhudson1
