try:
    import fwdpy.qtrait
    import fwdpy.fitness
    import unittest
    import array
    #Some basic setup
    rng = fwdpy.GSLrng(101)
    nlist = array.array('I',[1000]*10)
    f = fwdpy.qtrait.SpopGBRTrait()
    pops=fwdpy.SpopVec(10,1000)
    n = fwdpy.NothingSampler(len(pops))

    class SregionErrors(unittest.TestCase):
        """
        The gene-based region functions cannot allow effect size < 0.
        """
        def testGaussian(self):
            with self.assertRaises(RuntimeError):
                fwdpy.qtrait.evolve_regions_qtrait_sampler_fitness(rng,pops,n,f,nlist[0:],0,0.001,0.,[],[fwdpy.GaussianS(0,1,1,0.25)],[],1,0.025)
        def testConstant(self):
            with self.assertRaises(RuntimeError):
                fwdpy.qtrait.evolve_regions_qtrait_sampler_fitness(rng,pops,n,f,nlist[0:],0,0.001,0.,[],[fwdpy.ConstantS(0,1,1,-0.25)],[],1,0.025)
        def testExp(self):
            with self.assertRaises(RuntimeError):
                fwdpy.qtrait.evolve_regions_qtrait_sampler_fitness(rng,pops,n,f,nlist[0:],0,0.001,0.,[],[fwdpy.ExpS(0,1,1,-0.25)],[],1,0.025)
        def testGamma(self):
            with self.assertRaises(RuntimeError):
                fwdpy.qtrait.evolve_regions_qtrait_sampler_fitness(rng,pops,n,f,nlist[0:],0,0.001,0.,[],[fwdpy.GammaS(0,1,1,-0.25,2)],[],1,0.025)
        def testUniformLo(self):
            with self.assertRaises(RuntimeError):
                fwdpy.qtrait.evolve_regions_qtrait_sampler_fitness(rng,pops,n,f,nlist[0:],0,0.001,0.,[],[fwdpy.UniformS(0,1,1,-0.25,2)],[],1,0.025)
        def testUniformHi(self):
            with self.assertRaises(RuntimeError):
                fwdpy.qtrait.evolve_regions_qtrait_sampler_fitness(rng,pops,n,f,nlist[0:],0,0.001,0.,[],[fwdpy.UniformS(0,1,1,-0.2,-0.1)],[],1,0.025)
                
except ImportError:
    pass

if __name__ == '__main__':
    unittest.main()
