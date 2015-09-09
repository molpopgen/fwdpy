import math

## Testing
class Region(object):
    def __init__(self,float beg,float end,float weight,coupled=False):
        if math.isinf(beg):
            raise ValueError("fwdpy.Region: beg not finite")
        if math.isinf(end):
            raise ValueError("fwdpy.Region: end not finite")
        if math.isinf(weight):
            raise ValueError("fwdpy.Region: weight not finite")
        if math.isnan(beg):
            raise ValueError("fwdpy.Region: beg not a number")
        if math.isnan(end):
            raise ValueError("fwdpy.Region: end not a number")
        if math.isnan(weight):
            raise ValueError("fwdpy.Region: weight not a number")
        if weight <= 0.0:
            raise ValueError("fwdpy.Region: weight <= 0.0")
        self.b=float(beg)
        self.e=float(end)
        self.w=float(weight)
        self.c=coupled
        
class Sregion(Region):
    def __init__(self,float beg,float end,float weight,float h=1.0,coupled=False):
        if math.isinf(h):
            raise ValueError("fwdpy.Sregion: h not finite")
        if math.isnan(h):
            raise ValueError("fwdpy.Segion: h not a number")
        self.h=float(h)
        super(Sregion,self).__init__(beg,end,weight,coupled)
        
class GammaS(Sregion):
    def __init__(self,float beg,float end,float weight,float mean,float shape,float h=1.0,coupled=False):
        if math.isinf(mean):
            raise ValueError("fwdpy.GammaS: mean not finite")
        if math.isinf(shape):
            raise ValueError("fwdpy.GammaS: shape not finite")
        if math.isnan(mean):
            raise ValueError("fwdpy.GammaS: mean not a number")
        if math.isnan(shape):
            raise ValueError("fwdpy.GammaS: shape not a number")
        self.mean=float(mean)
        self.shape=float(shape)
        super(GammaS,self).__init__(beg,end,weight,h,coupled)
        
class ConstantS(Sregion):
    def __init__(self,float beg,float end,float weight,float s,float h=1.0,coupled=False):
        if math.isinf(s):
            raise ValueError("fwdpy.ConstantS: s not finite")
        if math.isnan(s):
            raise ValueError("fwdpy.ConstantS: s not a number")
        self.s=float(s)
        super(ConstantS,self).__init__(beg,end,weight,h,coupled)

class UniformS(Sregion):
    def __init__(self,float beg,float end,float weight,float lo,float hi,float h=1.0,coupled=False):
        if math.isinf(lo):
            raise ValueError("fwdpy.UniformS: lo not finite")
        if math.isinf(hi):
            raise ValueError("fwdpy.UniformS: hi not finite")
        if math.isnan(lo):
            raise ValueError("fwdpy.UniformS: lo not a number")
        if math.isnan(hi):
            raise ValueError("fwdpy.UniformS: hi not a number")
        self.lo=float(lo)
        self.hi=float(hi)
        super(ConstantS,self).__init__(beg,end,weight,h,coupled)

class ExpS(Sregion):
    def __init__(self,float beg,float end,float weight,float mean,float h=1.0,coupled=False):
        if math.isinf(mean):
            raise ValueError("fwdp.ExpS: mean not finite")
        if math.isnan(mean):
            raise ValueError("fwdpy.ExpS: mean not a number")
        self.mean=float(mean)
        super(ExpS,self).__init__(beg,end,weight,h,coupled)

class GaussianS(Sregion):
    def __init__(self,float beg,float end,float weight,float sd,float h=1.0,coupled=False):
        if math.isinf(sd):
            raise ValueError("fwdpy.GaussianS: sd not finite")
        if math.isnan(sd):
            raise ValueError("fwdpy.GaussianS: sd not a number")
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
            
