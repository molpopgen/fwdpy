import unittest
import fwdpy
import numpy as np

#Mock objects
nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
rregions = [fwdpy.Region(0,3,1)]
rng = fwdpy.GSLrng(100)
popsizes = np.array([1000]*5,dtype=np.uint32)

class EvolveRegions(unittest.TestCase):
    """
    Make sure that functions handle negative/non-finite rates by raising RuntimeError
    """
    def test_negativeNeutralMutationRate(self):
        with self.assertRaises(RuntimeError):
            pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],-0.001,0.0001,0.001,nregions,sregions,rregions)
    def test_nonFiniteNeutralMutationRate(self):
        with self.assertRaises(RuntimeError):
            pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],np.nan,0.0001,0.001,nregions,sregions,rregions)
    def test_negativeNonNeutralMutationRate(self):
        with self.assertRaises(RuntimeError):
            pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,-0.0001,0.001,nregions,sregions,rregions)
    def test_nonFiniteNonNeutralMutationRate(self):
        with self.assertRaises(RuntimeError):
            pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,np.nan,0.001,nregions,sregions,rregions)
    def test_negativeRecRate(self):
        with self.assertRaises(RuntimeError):
            pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.0001,-0.001,nregions,sregions,rregions)
    def test_nonFiniteRecRate(self):
        with self.assertRaises(RuntimeError):
            pops = fwdpy.evolve_regions(rng,1,1000,popsizes[0:],0.001,0.001,np.inf,nregions,sregions,rregions)

if __name__ == '__main__':
    unittest.main()
