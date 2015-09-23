from libcpp.vector cimport vector
import pandas
import fwdpy
    
def process_regions(list regions):
    """
    Process a list of objects of :class:`Region`

    :param regions: a list of objects of :class:`Region`

    :return: a pandas.DataFrame consisting of the beginning ('beg'), end ('end') and weight ('weight') for each elememt in l

    A user will generally not call this function.  Rather, it is used internally by things like :func:`fwdpy.fwdpy.evolve_regions`.
    """
    starts=list()
    stops=list()
    weights=list()
    for i in range(len(regions)):
        if isinstance(regions[i],fwdpy.Region):
            starts.append(regions[i].b)
            stops.append(regions[i].e)
            weights.append(regions[i].w)
        else:
            raise ValueError("invalid callback type encountered")
    return pandas.DataFrame({'beg':starts,'end':stops,'weight':weights})

def process_sregion_callbacks( shwrappervec v, list sregions ):
    """
    Process a list of objects of :class:`Sregion`

    :param v: a :class:`shwrappervec`
    :param sregions: a list of objects of :class:`Sregion`

    :return: Nothing. This function populations v with necessary callbacks for the C++ code to run the desired model.

    A user will generally not call this function.  Rather, it is used internally by things like :func:`fwdpy.fwdpy.evolve_regions`.
    """
    cdef shmodel temp
    for i in range(len(sregions)):
        if not isinstance(sregions[i],fwdpy.Sregion):
            raise ValueError("fwdpy.process_sregion_callbacks: invalid object type")
        if isinstance(sregions[i],fwdpy.GammaS):
            make_gamma_s(&temp,sregions[i].mean,sregions[i].shape)
            make_constant_h(&temp,sregions[i].h)
            v.vec.push_back(temp)
        elif isinstance(sregions[i],fwdpy.ConstantS):
            make_constant_s(&temp,sregions[i].s)
            make_constant_h(&temp,sregions[i].h)
            v.vec.push_back(temp)
        elif isinstance(sregions[i],fwdpy.ExpS):
            make_exp_s(&temp,sregions[i].mean)
            make_constant_h(&temp,sregions[i].h)
            v.vec.push_back(temp)
        elif isinstance(sregions[i],fwdpy.UniformS):
            make_uniform_s(&temp,sregions[i].lo,sregions[i].hi)
            make_constant_h(&temp,sregions[i].h)
            v.vec.push_back(temp)
        elif isinstance(sregions[i],fwdpy.GaussianS):
            make_gaussian_s(&temp,sregions[i].sd)
            make_constant_h(&temp,sregions[i].h)
            v.vec.push_back(temp)
        else:
            raise ValueError("fwdpy.process_sregion_callbacks: unsupported Sregion type")
            
