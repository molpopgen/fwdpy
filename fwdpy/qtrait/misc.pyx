def popstats( singlepop pop ):
    """
    Get statistics like VG, etc. from a population

    :return: A pandas.DataFrame.

    :rtype: pandas.DataFrame
    """
    return pandas.DataFrame(qtrait_pop_props(pop.pop.get()).items(),columns=['stat','value'])

