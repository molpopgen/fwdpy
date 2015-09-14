import fwdpy

## These functions are intentionally undocumented

def parse_slim_mutrate(list lines):
    start = lines.index("#MUTATION RATE")
    return float(lines[start+1])
        
def parse_slim_muttypes(list lines):

    start = lines.index("#MUTATION TYPES")
    nbloc = [i for i,j in enumerate(lines[start+1:]) if j[0]=='#']
    end = len(lines)
    if len(nbloc):
        end=nbloc[0]

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

def parse_slim_elements(list lines):
    start = lines.index("#GENOMIC ELEMENT TYPES")
    nbloc = [i for i,j in enumerate(lines[start+1:]) if j[0]=='#']
    end = len(lines)
    if len(nbloc):
        end=nbloc[0]

    gtypes = {}
    for i in lines[start+1:start+1+end]:
        t=i.split()
        inner={}
        for j in range(0,len(t[1:]),2):
            inner[t[1:][j]]=float(t[1:][j+1])
        gtypes[t[0]]=inner

    return gtypes

def parse_slim_organization(list lines,dict mtypes, dict elements, double mutrate):
    start = lines.index("#CHROMOSOME ORGANIZATION")
    nbloc = [i for i,j in enumerate(lines[start+1:]) if j[0]=='#']
    end = len(lines)
    if len(nbloc):
        end=nbloc[0]

    #need to get sum of all weights
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
            if mt[1] == 'f':
                if mt[2] == 0.:  #is a neutral mutation
                    mun = mun + mutrate*weight*(eend-ebeg)
                    nregions.append(fwdpy.Region(ebeg,eend,mutrate*weight))
                else:
                    mus = mus + mutrate*weight*(eend-ebeg)
                    sregions.append(fwdpy.ConstantS(ebeg,eend,mutrate*weight,mt[2],mt[0]))
            elif mt[1] == 'e':
                mus = mus + mutrate*weight*(eend-ebeg)
                sregions.append(fwdpy.ExpS(ebeg,eend,mutrate*weight,mt[2],mt[0]))
            elif mt[1] == 'g':
                mus = mus + mutrate*weight*(eend-ebeg)
                sregions.append(fwdpy.GammaS(ebeg,eend,mutrate*weight,mt[2],mt[3],mt[0]))
            else:
                raise RuntimeError("invalid DFE encountered")
    return {'nregions':nregions,'sregions':sregions,'mu_neutral':mun,'mu_selected':mus}

def parse_slim_recrates(list lines):
    start = lines.index("#RECOMBINATION RATE")
    nbloc = [i for i,j in enumerate(lines[start+1:]) if j[0]=='#']
    end = len(lines)
    if len(nbloc):
        end=nbloc[0]

    regions = []
    recrate = 0.
    laststart=float(1.0)
    for i in lines[start+1:start+1+end]:
        t=i.split()
        stop = float(t[0])
        weight=float(t[1])
        if weight > 0.:
            regions.append(fwdpy.Region(laststart,stop,weight))
            recrate = recrate + weight*(stop-laststart+1.)
            laststart=float(stop)

    return {'recrate':recrate,'recregions':regions}
