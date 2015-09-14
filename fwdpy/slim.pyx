import os.path
import internal

#Implementation details found in internal/slim.pyx

def readslim(const char * filename):
    """
    Process a "SliM" input file
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

