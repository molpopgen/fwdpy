##Code for adding a mutation into a population

def add_mutation(PopType p,
                 const vector[size_t] & indlist,
                 const vector[short] & clist,
                 double pos,
                 double s,
                 double h):
    """
    Add a mutation to a population at a specific position.

    :param p: A :class:`fwdpy.fwdpy.PopType`.
    :param indlist: A list of individuals to whom the mutation will be added.  Must be unsigned indexes [0,N).
    :param clist: A list of chromosomes to which to add the new mutation in each individual. 0 = chrom1, 1 = chrom2, 2 = both.
    :param pos: The position of the new mutation
    :param s: The selection coefficient/effect size of the new mutation
    :param h: The dominance coefficient of the new mutation

    .. note:: Only :class:`fwdpy.fwdpy.Spop` is currently supported

    .. note:: RuntimeError will be thrown if pos already exists in p or if any element in clist and/or indlist is out of range.

    Example:

    >>> import fwdpy as fp
    >>> import fwdpy.views 
    >>> #init a monomorphic pop
    >>> p = fp.SpopVec(1,1000)
    >>> x = fp.add_mutation(p[0],[0,10,20],[0,1,2],0.01,0.0,1.0)
    >>> #x will have value 0 b/c it is first mutation entered into pop
    >>> x
    0
    >>> #Note here how the position returned is not exactly
    >>> #the position input.  This occurs because not all 
    >>> #floating-point values are representable exactly.
    >>> fp.views.view_mutations(p[0])
    [position: 0.00999999977648, count: 4, origin time: 1, fixation time: None, s: 0.0, h: 1.0, neutral: True, label: 0.]
    """
    if isinstance(p,Spop):
        return add_mutation_cpp((<Spop>p).pop.get(),indlist,clist,<double>pos,s,h)
    else:
        raise RuntimeError("PopType not supported")

