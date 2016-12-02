cdef class MutationView(object):
    """
    An immutable view of a mutation.
    """
    def __cinit__(self,pos,n,g,ftime, s,
            h, neutral,label,key):
        self.pos=pos
        self.n=n
        self.g=g
        self.ftime=ftime
        self.s=s
        self.h=h
        self.neutral=neutral
        self.label=label
        self.mut_key=key

        #check that types are correct
        for i in [self.pos,self.s,self.h]:
            if i is not None:
                try:
                    float(i)
                except ValueError:
                    raise ValueError("incorrect value type")
        for i in [self.n,self.g,self.ftime,self.label,self.mut_key]:
            if i is not None:
                try:
                    int(i)
                except ValueError:
                    raise ValueError("incorrect value type")
    def __repr__(self):
        """
        Returns a "pretty" string representation
        """
        r = b'position: ' + format(self.pos) + b', '
        r += b'count: ' + format(self.n) + b', '
        r += b'origin time: ' + format(self.g) + b', '
        r += b'fixation time: ' + format(self.ftime) + b', '
        r += b's: ' + format(self.s) + b', '
        r += b'h: ' + format(self.h) + b', '
        r += b'neutral: ' + format(self.neutral) + b', '
        r += b'label: ' + format(self.label) + b', '
        r += b'mut_key: ' + format(self.mut_key)
        return r
    def as_dict(self):
        """
        :return: The mutation data represented as dict

        :rtype: dict
        """
        return {'position':self.pos,'count':self.n,'origin':self.g,
                'fixation':self.ftime,'s':self.s,'h':self.h,'neutral':self.neutral,
                'label':self.label,'mut_key':self.mut_key}

cdef class GameteView(object):
    """
    An immutable view of a gamete.

    The class is iterable, yielding :class:`fwdpy.views.MutationView`.  When iterating,
    neutral mutations come out first, followed by selected mutations.

    In many use cases, it may make more sense to manually iterate over neutral and 
    selected separately.
    """
    def __cinit__(self,list neutral_mutations,list selected_mutations, int count,key):
        self.neutral=neutral_mutations
        self.selected=selected_mutations
        self.n=count
        self.gam_key=key
    def __repr__(self):
        r = b'Neutral variants:\n' + b'\n'.join([str(i) for i in self.neutral]) + b'\n'
        r += b'Selected variants:\n' + b'\n'.join([str(i) for i in self.selected]) + b'\n'
        return r
    def __iter__(self):
        return iter(self.neutral+self.selected)
    def __next__(self):
        return next(self.neutral+self.selected)
    def __getitem__(self, int i):
        return (self.neutral+self.selected)[i] 
    def __len__(self):
        return len(self.neutral+self.selected)
    def num_neutral(self):
        return len(self.neutral)
    def num_selected(self):
        return len(self.selected)
    def as_list(self):
        """
        :return: All mutation data
        :rtype: list of dict
        """
        muts=[i.as_dict() for i in self]
        for m in muts:
            m['gam_key']=self.gam_key
        return muts

cdef class DiploidView(object):
    """
    Immutable view of a diploid.
    """
    def __cinit__(self,GameteView a,GameteView b,float genetic_value,float env_value,float fitness,key):
        self.first=a
        self.second=b
        self.g=genetic_value
        self.e=env_value
        self.w=fitness
        self.dip_key=key
    def as_list(self):
        """
        :return: All gamete data
        :return: list of dict
        """
        muts=self.first.as_list()+self.second.as_list()
        for i in muts:
            i['dip_key']=self.dip_key
        return muts

cdef class MultiLocusDiploidView(object):
    """
    Immutable view of a diploid from a multi-locus simulation.
    """
    def __cinit__(self,list a,list b,float genetic_value, float env_value, float fitness,key):
        self.first=a
        self.second=b
        self.g=genetic_value
        self.e=env_value
        self.w=fitness
        self.dip_key=key
    def __addkeys__(self,list l):
        locus=0
        rv=[]
        for li in l:
            for lii in li.as_list():
                lii['locus_key']=locus
                lii['dip_key']=self.dip_key
                rv.append(lii)
            locus+=1
        return rv
    def as_list(self):
        """
        :return: All mutation data for all loci
        :rtype: list of dict
        """
        loci = self.__addkeys__(self.first)
        loci2 = self.__addkeys__(self.second)
        return loci+loci2

cdef MutationView empty_MutationView():
    return MutationView(None,None,None,None,None,None,None,None,None)

include "view_mutations.pyx"
include "view_fixations.pyx"
include "view_gametes.pyx"
include "view_diploids.pyx"
include "view_diploid_traits.pyx"
