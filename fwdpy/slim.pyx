import os.path
import internal

#Implementation details found in internal/slim.pyx

def readslim(const char * filename):
    """
    Process a "SLiM" input file.

    This function parses the input files for Phillip Messer's "SLiM" software (http://messerlab.org/software/).

    Thus function only parses the following blocks in a file:
    1. MUTATION RATE
    2. MUTATION TYPES
    3. RECOMBINATION RATE
    4. GENOMIC ELEMENTS
    5. CHROMOSOME ORGANIZATION

    :param filename: An input file for the SLiM simulation software

    :return: A dict containing the following data members:

    'nregions' A list of type :class:`fwdpy.fwdpy.Region`.  This list represents regions where neutral muations occur.
    
    'sregions' A list of types inheriting from :class:`fwdpy.fwdpy.Sregion`.  This list represents regions where mutations affecting fitness occur.
    
    'recregion' A list of type :class:'fwdpy.fwdpy.Region` representing recombination rate variation
    
    'mu_neutral' is the total mutation rate to neutral variants (per gamete, per generation).
    
    'mu_selected' is the total mutation rate to selected variants (per gamete, per generation).
    
    'recrate' is the total recombination rate (per diploid, per generation).
    
    :rtype: A dictionary containing the parsing results.
    """
    if not os.path.isfile(filename):
        raise IOError("fwdpy.readslim: "+filename+" does not appear to be a regular file")
    f = open(filename,'r')
    lines = f.read().splitlines()
    f.close()
    ##First character must be a hash
    if not lines[0][0] == '#':
        raise RuntimeError(filename+" does not appear to be a slim input file. The first line is "+lines[0])

    mtypes = internal.parse_slim_muttypes(lines)
    recdata = internal.parse_slim_recrates(lines)
    gtypes = internal.parse_slim_elements(lines)
    mutrate = internal.parse_slim_mutrate(lines)
    org = internal.parse_slim_organization(lines,mtypes,gtypes,mutrate)

    return {'nregions':org['nregions'],
            'sregions':org['sregions'],
            'mu_neutral':org['mu_neutral'],
            'mu_selected':org['mu_selected'],
            'recregions':recdata['recregions'],
            'recrate':recdata['recrate'],
            }

