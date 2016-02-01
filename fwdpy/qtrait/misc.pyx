def popstats( singlepop pop ):
    """
    Get statistics like VG, etc. from a population

    :return: List of dict
    """
    cdef vector[qtrait_stats_cython] temp=convert_qtrait_stats(pop.pop.get())
    rv = []
    for i in range(temp.size()):
        rv.append({'stat':temp[i].stat,'value':temp[i].value,'generation':temp[i].generation});
    return rv

