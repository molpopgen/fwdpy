import os.path
import fwdpy

## These functions are intentionally undocumented
def get_next_block_starts( lines, start):
    nbloc = [i for i,j in enumerate(lines[start+1:]) if j[0]=='#']
    end = len(lines)
    if len(nbloc):
        end=nbloc[0]
    return end

def parse_slim_mutrate(lines):
    start = lines.index("#MUTATION RATE")
    return float(lines[start+1])
        
def parse_slim_muttypes(lines):
    start = lines.index("#MUTATION TYPES")
    end = get_next_block_starts(lines,start)

    mtypes = {}
    for i in lines[start+1:start+1+end]:
        t=i.split()
        if not (t[2] == 'e' or t[2] == 'g' or t[2] == 'f'):
            raise RuntimeError("invalid mutation type specified in SLiM input file")
        else:
            if t[2]=='e':
                if len(t[1:]) != 3:
                    raise RuntimeError("invalid specification of exponential DFE: "+i)
            elif t[2]=='g':
                if len(t[1:]) != 4:
                    raise RuntimeError("invalid specification of gamma DFE: "+i)
            elif t[2]=='f':
                if len(t[1:]) != 3:
                    raise RuntimeError("invalid specification of fixed DFE: "+i)
            ##convert all elements except 1 to float
            for j in range(0,len(t[1:]),1):
                if j != 1:
                    t[j+1]=float(t[1:][j])
            mtypes[t[0]] = t[1:]
    return mtypes

def parse_slim_elements(lines):
    start = lines.index("#GENOMIC ELEMENT TYPES")
    end = get_next_block_starts(lines,start)
    gtypes = {}
    for i in lines[start+1:start+1+end]:
        t=i.split()
        inner={}
        for j in range(0,len(t[1:]),2):
            inner[t[1:][j]]=float(t[1:][j+1])
        gtypes[t[0]]=inner

    return gtypes

def parse_slim_organization(lines,mtypes, elements, mutrate):
    start = lines.index("#CHROMOSOME ORGANIZATION")
    end = get_next_block_starts(lines,start)
    #need to get sum of all weights per mutation type
    ttlweights = {}
    for key1 in elements:
        ttlweights[key1]=0.
        for key2 in elements[key1]:
            ttlweights[key1]=ttlweights[key1] + elements[key1][key2]
    nregions=[]
    sregions=[]
    mun=0
    mus=0
    for i in lines[start+1:start+1+end]:
        t=i.split()
        ebeg=float(t[1])
        eend=float(t[2])
        for key in elements[t[0]]:
            mt = mtypes[key]
            weight=elements[t[0]][key]/ttlweights[t[0]]
            ##In this block, we halve s or mean s,
            ##and double h, to convert from SLiM's
            ##fitness model of 1,1+hs,1+s to the
            ##1,1+sh,1+2s used here.
            if mt[1] == 'f':
                if mt[2] == 0.:  #is a neutral mutation
                    mun = mun + mutrate*weight*(eend-ebeg+1.)
                    nregions.append(fwdpy.Region(ebeg-1.,eend,mutrate*weight))
                else:
                    mus = mus + mutrate*weight*(eend-ebeg+1.)
                    sregions.append(fwdpy.ConstantS(ebeg-1.,eend,mutrate*weight,0.5*mt[2],2*mt[0]))
            elif mt[1] == 'e':
                mus = mus + mutrate*weight*(eend-ebeg+1.)
                sregions.append(fwdpy.ExpS(ebeg-1.,eend,mutrate*weight,0.5*mt[2],2*mt[0]))
            elif mt[1] == 'g':
                mus = mus + mutrate*weight*(eend-ebeg+1.)
                sregions.append(fwdpy.GammaS(ebeg-1.,eend,mutrate*weight,0.5*mt[2],mt[3],2*mt[0]))
            else:
                raise RuntimeError("invalid DFE encountered")
    return {'nregions':nregions,'sregions':sregions,'mu_neutral':mun,'mu_selected':mus}

def parse_slim_recrates(lines):
    start = lines.index("#RECOMBINATION RATE")
    end = get_next_block_starts(lines,start)
    regions = []
    recrate = 0.
    laststart=float(0.0)
    for i in lines[start+1:start+1+end]:
        t=i.split()
        stop = float(t[0])
        weight=float(t[1])
        if weight > 0.:
            ##NEED TO DOCUMENT THE + 1
            regions.append(fwdpy.Region(laststart,stop,weight))
            recrate = recrate + weight*(stop-laststart)
        laststart=float(stop)

    return {'recrate':recrate,'recregions':regions}


def readslim(filename):
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

    .. note:: This function is automatically imported into fwdpy.
    """
    if not os.path.isfile(filename):
        raise IOError("fwdpy.readslim: "+filename+" does not appear to be a regular file")
    f = open(filename,'r')
    lines = f.read().splitlines()
    ##Remove all leading and trailing whitespace, just in case the user added some in
    lines = [i.strip() for i in lines]
    f.close()
    ##First character must be a hash
    if not lines[0][0] == '#':
        raise RuntimeError(filename+" does not appear to be a slim input file. The first line is "+lines[0])

    mtypes = parse_slim_muttypes(lines)
    recdata = parse_slim_recrates(lines)
    gtypes = parse_slim_elements(lines)
    mutrate = parse_slim_mutrate(lines)
    org = parse_slim_organization(lines,mtypes,gtypes,mutrate)

    return {'nregions':org['nregions'],
            'sregions':org['sregions'],
            'mu_neutral':org['mu_neutral'],
            'mu_selected':org['mu_selected'],
            'recregions':recdata['recregions'],
            'recrate':recdata['recrate'],
            }

