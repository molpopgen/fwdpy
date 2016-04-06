import unittest
import fwdpy as fp
import fwdpy.demography as demog
import numpy as np

class test_make_mpopvec(unittest.TestCase):
    """
    Straightfoward tests to make sure
    that copying from singlepop -> metapop 
    is working.
    """
    def testLength(self):
        s = fp.popvec(64,100)
        m = demog.make_mpopvec(s)
        self.assertEqual(len(m),len(s))
    def testSizes(self):
        s = fp.popvec(64,100)
        m = demog.make_mpopvec(s)
        for i in m:
            self.assertEqual(i.popsizes(),[100])

class test_exceptions(unittest.TestCase):
    """
    fwdpy.demography relies on fwdpp's return values
    for throwing exceptions.  See fwdpy/demograpy/demography.cc
    for details
    """
    def testCopyOutOfRange(self):
        m = fp.mpopvec(64,[100])
        ##There is only 1 deme, so index of 1 is out of range
        with self.assertRaises(IndexError):
            demog.copy_pop(m,1)
    def testMergeOutOfRange(self):
        m = fp.mpopvec(64,[100])
        #Make identical copy of
        #first deme in each pop
        demog.copy_pop(m,0)
        with self.assertRaises(IndexError):
            #Now, try to merge
            #but second index
            #out of range:
            demog.merge_pops(m,0,2)
    def testMergeDemeOntoItself(self):
        m=fp.mpopvec(64,[100,100])
        with self.assertRaises(RuntimeError):
            ##Cannot merge deme into itself:
            demog.merge_pops(m,1,1)
    def testRemoveOutOfRange(self):
        m = fp.mpopvec(64,[100,100])
        #remove the second deme
        demog.remove_pop(m,1)
        #try again = exception
        with self.assertRaises(IndexError):
            demog.remove_pop(m,1)
    def testSwapOutOfRange(self):
        m = fp.mpopvec(64,[100,100])
        with self.assertRaises(IndexError):
            demog.swap_pops(m,0,2)
    def testSplitOutOfRange(self):
        m = fp.mpopvec(64,[100])
        rng=fp.GSLrng(101)
        with self.assertRaises(IndexError):
            demog.split_pops(rng,m,2,50)
    def testAdmixOutOfRange(self):
        m = fp.mpopvec(64,[100,100])
        rng=fp.GSLrng(101)
        with self.assertRaises(IndexError):
            demog.admix_pops(rng,m,0,2,0.1,100)
    def testAdmixInvalidProportion(self):
        """
        Test requirement that admix prop is 0 <= p <= 1
        """
        m = fp.mpopvec(64,[100,100])
        rng=fp.GSLrng(101)
        with self.assertRaises(RuntimeError):
            demog.admix_pops(rng,m,0,1,-0.1,100)
