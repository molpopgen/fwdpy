def popstats( singlepop pop ):
    """
    Get statistics like VG, etc. from a population

    :return: List of dict
    """
    return convert_qtrait_stats(pop.pop.get())


