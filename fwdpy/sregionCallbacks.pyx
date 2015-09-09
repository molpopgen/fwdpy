## Testing
class Region(object):
    def __init__(self,beg,end,weight,coupled=False):
        self.b=float(beg)
        self.e=float(end)
        self.w=float(weight)
        self.c=coupled

class Sregion(Region):
    def __init__(self,beg,end,weight,h=1.0,coupled=False):
        self.h=float(h)
        super(Sregion,self).__init__(beg,end,weight,coupled)
        
class GammaS(Sregion):
    def __init__(self,beg,end,weight,mean,shape,h=1.0,coupled=False):
        self.mean=float(mean)
        self.shape=float(shape)
        super(GammaS,self).__init__(beg,end,weight,h,coupled)

class ConstantS(Sregion):
    def __init__(self,beg,end,weight,s,h=1.0,coupled=False):
        self.s=float(s)
        super(ConstantS,self).__init__(beg,end,weight,h,coupled)

class UniformS(Sregion):
    def __init__(self,beg,end,weight,lo,hi,h=1.0,coupled=False):
        self.lo=float(lo)
        self.hi=float(hi)
        super(ConstantS,self).__init__(beg,end,weight,h,coupled)

class ExpS(Sregion):
    def __init__(self,beg,end,weight,mean,h=1.0,coupled=False):
        self.mean=float(mean)
        super(ExpS,self).__init__(beg,end,weight,h,coupled)

class GaussianS(Sregion):
    def __init__(self,beg,end,weight,sd,h=1.0,coupled=False):
        self.sd=float(sd)
        super(GaussianS,self).__init__(beg,end,weight,h,coupled)

def process_regions(list l):
    starts=list()
    stops=list()
    weights=list()
    for i in range(len(l)):
        if isinstance(l[i],Region):
            starts.append(l[i].b)
            stops.append(l[i].e)
            wTemp=l[i].w
            if(l[i].c):
                wTemp = wTemp*(l[i].e-l[i].b)
            weights.append(wTemp)
        else:
            raise ValueError("invalid callback type encountered")
    return pandas.DataFrame({'beg':starts,'end':stops,'weight':weights})

def process_sregion_callbacks( shwrappervec v, list sregions ):
    cdef shmodel temp
    for i in range(len(sregions)):
        if not isinstance(sregions[i],Sregion):
            raise ValueError("fwdpy.process_sregion_callbacks: invalid object type")
        if isinstance(sregions[i],GammaS):
            make_gamma_s(&temp,sregions[i].mean,sregions[i].shape)
            make_constant_h(&temp,sregions[i].h)
            v.vec.push_back(temp)
        elif isinstance(sregions[i],ConstantS):
            make_constant_s(&temp,sregions[i].s)
            make_constant_h(&temp,sregions[i].h)
            v.vec.push_back(temp)
        elif isinstance(sregions[i],ExpS):
            make_exp_s(&temp,sregions[i].mean)
            make_constant_h(&temp,sregions[i].h)
            v.vec.push_back(temp)
        elif isinstance(sregions[i],UniformS):
            make_uniform_s(&temp,sregions[i].lo,sregions[i].hi)
            make_constant_h(&temp,sregions[i].h)
            v.vec.push_back(temp)
        elif isinstance(sregions[i],GaussianS):
            make_gaussian_s(&temp,sregions[i].sd)
            make_constant_h(&temp,sregions[i].h)
            v.vec.push_back(temp)
        else:
            raise ValueError("fwdpy.process_sregion_callbacks: unsupported Sregion type")
            
