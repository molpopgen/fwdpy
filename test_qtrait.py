import unittest

import array
import fwdpy.qtrait

#Some basic setup
rng = fwdpy.GSLrng(101)
nlist = array.array('I',[1000]*10)

class testGaussianError(unittest.TestCase):
    def testException(self):
        with self.AssertRaises(RuntimeError):
            pops = fwdpy.qtrait.evolve_gbr(rng,1,1000,nlist[0:],0,0.001,0.,[],[fwdpy.GaussianS(0,1,1,0.25)],[],0.025)
