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
