## Testing
class Region(object):
    def __init__(self,beg,end,weight,coupled):
        self.b=beg
        self.e=end
        self.w=weight
        self.c=coupled

class Sregion(Region):
    def __init__(self,beg,end,weight,h=1.0,coupled=False):
        self.h=h
        super(Sregion,self).__init__(beg,end,weight,coupled)
        
class GammaS(Sregion):
    def __init__(self,beg,end,weight,mean,shape,h=1.0,coupled=False):
        self.mean=mean
        self.shape=shape
        super(GammaS,self).__init__(beg,end,weight,h,coupled)

class ConstantS(Sregion):
    def __init__(self,beg,end,weight,s,h=1.0,coupled=False):
        self.s=s
        super(ConstantS,self).__init__(beg,end,weight,h,coupled)

class UniformS(Sregion):
    def __init__(self,beg,end,weight,lo,hi,h=1.0,coupled=False):
        self.lo=lo
        self.hi=hi
        super(ConstantS,self).__init__(beg,end,weight,h,coupled)

class ExpS(Sregion):
    def __init__(self,beg,end,weight,mean,h=1.0,coupled=False):
        self.mean=mean
        super(ExpS,self).__init__(beg,end,weight,h,coupled)

class GaussianS(Sregion):
    def __init__(self,beg,end,weight,sd,h=1.0,coupled=False):
        self.sd=sd
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

def process_sregion_callbacks(list l):
    rv = list()
    starts=list()
    stops=list()
    weights=list()
    for i in range(len(l)):
        if isinstance(l[i],GammaS):
            tsh = shwrapper()
            make_gamma_s(&tsh.thisptr,l[i].mean,l[i].shape)
            make_constant_h(&tsh.thisptr,l[i].h)
            rv.append(tsh)
        elif isinstance(l[i],ConstantS):
            tsh=shwrapper()
            make_constant_s(&tsh.thisptr,l[i].s)
            make_constant_h(&tsh.thisptr,l[i].h)
            rv.append(tsh)
        elif isinstance(l[i],ExpS):
            tsh=shwrapper()
            make_exp_s(&tsh.thisptr,l[i].mean)
            make_constant_h(&tsh.thisptr,l[i].h)
            rv.append(tsh)
        elif isinstance(l[i],UniformS):
            tsh=shwrapper()
            make_uniform_s(&tsh.thisptr,l[i].lo,l[i].hi)
            make_constant_h(&tsh.thisptr,l[i].h)
            rv.append(tsh)
        elif isinstance(l[i],GaussianS):
            tsh=shwrapper()
            make_gaussian_s(&tsh.thisptr,l[i].sd)
            make_constant_h(&tsh.thisptr,l[i].h)
            rv.append(tsh)
        else:
            raise ValueError("invalid callback type encountered")
    return [process_regions(l),rv]
