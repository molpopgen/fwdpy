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
    fwdpy.demography relies on fwdpp's to return values
    for throwing exceptions.  See fwdpy/demograpy/demography.cc
    for details
    """
    def testCopyOutOfRange(self):
        s = fp.popvec(64,100)
        m = demog.make_mpopvec(s)
        ##There is only 1 deme, so index of 1 is out of range
        with self.assertRaises(IndexError):
            demog.copy_pop(m,1)
        with self.assertRaises(IndexError):
            #This is also raise, as it'll get coverted to max(unsigned)-1:
            demog.copy_pop(m,-1)
    def testMergeOutOfRange(self):
        s = fp.popvec(64,100)
        m = demog.make_mpopvec(s)
        #Make identical copy of
        #first deme in each pop
        demog.copy_pop(m,0)
        with self.assertRaises(IndexError):
            #Now, try to merge
            #but second index
            #out of range:
            demog.merge_pops(m,0,2)
        with self.assertRaises(RuntimeError):
            ##Cannot merge deme into itself:
            demog.merge_pops(m,1,1)
    def testRemoveOutOfRange(self):
        s = fp.popvec(64,100)
        m = demog.make_mpopvec(s)
        #copy a deme
        demog.copy_pop(m,0)
        #remove the new deme
        demog.remove_pop(m,1)
        #try again = exception
        with self.assertRaises(IndexError):
            demog.remove_pop(m,1)
