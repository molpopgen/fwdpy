Complex demographic operations
======================================================

In addition to arbitrary changes in population size over time *fwdpy* supports complex demographic operations on metapopulations.  These operations act on :class:`fwdpy.fwdpy.mpopvec` objects and are carried out by the demography module:

>>> import fwdpy.demography as demog

The relevant functions are:

1. :func:`fwdpy.demography.make_mpopvec`
2. :func:`fwdpy.demography.copy_pop`
3. :func:`fwdpy.demography.merge_pops`
4. :func:`fwdpy.demography.remove_pop`
5. :func:`fwdpy.demography.split_pops`
6. :func:`fwdpy.demography.admix_pops`
7. :func:`fwdpy.demography.swap_pops`
   
