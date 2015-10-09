from cython.operator cimport dereference as deref,preincrement as inc
from fwdpy.fwdpy cimport poptype,singlepop,metapop
import fwdpy
import pandas as pd
import numpy as np

cdef getmuts_details(cpplist[popgenmut].iterator itr,cpplist[popgenmut].iterator end, float twoN, unsigned nmuts,bint all):
    buff = np.array([np.nan]*(4*nmuts),dtype=np.float64)
    i=0
    while itr != end:
        if all is True or deref(itr).neutral is not 1:
            buff[i]=deref(itr).pos
            buff[i+1]=deref(itr).n/twoN
            buff[i+2]=deref(itr).s
            buff[i+3]=deref(itr).h
            i+=4
        inc(itr)
    a=np.ndarray(shape=(i/4,4),buffer=buff[:i])
    return pandas.DataFrame(a,columns=["pos","freq","esize","h"])

##TODO: change to const vector once Cython implements those types
cdef add_fixations( vector[popgenmut] & fixations, const vector[unsigned] & ftimes, float twoN, bint all ):
    buff = np.array([np.nan]*(5*ftimes.size()),dtype=np.float64)
    cdef vector[popgenmut].iterator fbeg = fixations.begin()
    cdef vector[popgenmut].iterator fend = fixations.end()
    i=0
    j=0
    while fbeg != fend:
        if all is True or deref(fbeg).neutral is not 1:
            buff[i]=deref(fbeg).pos
            buff[i+1]=deref(fbeg).n/twoN
            buff[i+2]=deref(fbeg).s
            buff[i+3]=deref(fbeg).h
            buff[i+4]=ftimes[j]
            i+=5
        inc(fbeg)
        j+=1
    a=np.ndarray(shape=(i/5,5),buffer=buff[:i])
    return pd.DataFrame(a,columns=["pos","count","esize","h","ftime"])

def getmuts_singlepop(singlepop pop,bint all,bint fixations):
    cdef cpplist[popgenmut].iterator itr = pop.pop.get().mutations.begin()
    cdef cpplist[popgenmut].iterator end = pop.pop.get().mutations.end()
    segregating = getmuts_details(itr,end,2*pop.popsize(),pop.pop.get().mutations.size(),all)
    if fixations is False:
        return segregating
    else:
        return pd.concat([segregating,add_fixations(pop.pop.get().fixations,pop.pop.get().fixation_times,2*pop.popsize(),all)])

def getmuts_metapop(metapop pop,bint all,bint fixations):
    cdef cpplist[popgenmut].iterator itr = pop.mpop.get().mutations.begin()
    cdef cpplist[popgenmut].iterator end = pop.mpop.get().mutations.end()
    segregating = getmuts_details(itr,end,2*sum(pop.posizes()),pop.mpop.get().mutations.size(),all)
    if fixations is False:
        return segregating
    else:
        return pd.concat([segregating,add_fixations(pop.mpop.get().fixations,pop.mpop.get().fixation_times,2*pop.popsize(),all)])
