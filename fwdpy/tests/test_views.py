import unittest
import fwdpy
import fwdpy.views
import numpy as np

nregions = [fwdpy.Region(0,1,1),fwdpy.Region(2,3,1)]
sregions = [fwdpy.ExpS(1,2,1,-0.1),fwdpy.ExpS(1,2,0.01,0.001)]
rregions = [fwdpy.Region(0,3,1)]
rng = fwdpy.GSLrng(100)
N=1000
NGENS=10000
popsizes = np.array([N],dtype=np.uint32)
popsizes=np.tile(popsizes,NGENS)
pops = fwdpy.evolve_regions(rng,1,N,popsizes[0:],0.01,0.01,0.001,nregions,sregions,rregions)

#The sum of the gamete counts must be 2*(deme size):
#mpops = fwdpy.evolve_regions_split(rng,pops,popsizes[0:],popsizes[0:],0.001,0.0001,0.001,nregions,sregions,rregions,[0]*2)

class test_singlepop_views(unittest.TestCase):
    def testNumGametes(self):
        gams = fwdpy.views.view_gametes(pops[0])
        nsingle=0
        for i in gams: nsingle += i.n
        self.assertEqual(nsingle,2000)
    def testDipsize(self):
        dips_single = fwdpy.views.view_diploids(pops[0],[0,1,2])
        self.assertEqual(len(dips_single),3)
    def testBigDips(self):
        temp = fwdpy.views.view_diploids(pops[0],list(range(pops[0].popsize())))
        self.assertEqual(len(temp),pops[0].popsize())
    def testException1(self):
        with self.assertRaises(IndexError):
            fwdpy.views.view_diploids(pops[0],[pops[0].popsize()])
    def testFixationViews(self):
        temp = fwdpy.views.view_fixations(pops[0])
#class test_metapop_views(unittest.TestCase):
#    def testNumGametes(self):
#        gams = fwdpy.views.view_gametes(mpops[0],0) 
#        nmeta=0
#        for i in gams: nmeta += i['n']
#        self.assertEqual(nmeta,2000)
#    def testDemeException1(self):
#        with self.assertRaises(IndexError):
#            [fwdpy.views.view_gametes(i,2) for i in mpops]
#    def testDemeException2(self):
#        with self.assertRaises(IndexError):
#            [fwdpy.views.view_mutations(i,2) for i in mpops]
#    def testDemeException3(self):
#        with self.assertRaises(IndexError):
#            [fwdpy.views.view_diploids(i,[0,1,2],2) for i in mpops]
#    def testIndException1(self):
#        with self.assertRaises(IndexError):
#            [fwdpy.views.view_diploids(i,[0,1,2,N],1) for i in mpops]
#    def testNone1(self):
#        with self.assertRaises(RuntimeError):
#            [fwdpy.views.view_diploids(i,[0,1,2]) for i in mpops]
#    def testNone2(self):
#        with self.assertRaises(RuntimeError):
#            [fwdpy.views.view_gametes(i) for i in mpops]
#    def testNone3(self):
#        with self.assertRaises(RuntimeError):
#            [fwdpy.views.view_mutations(i) for i in mpops]
            
if __name__ == '__main__':
    unittest.main()
