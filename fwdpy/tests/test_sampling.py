import unittest
import fwdpy as fp
import numpy as np
import pandas as pd


class test_sample_details(unittest.TestCase):
    def test_CompareUsingView(self):

        nregions = [fp.Region(beg=0,end=1,weight=1)]

        sregions = [fp.ConstantS(beg=-1,end=0,weight=1,s=-0.05,h=1),
                    fp.ConstantS(beg=1,end=2,weight=1,s=-0.05,h=1)]

        recregions = [fp.Region(beg=-1,end=2,weight=1)]

        N=1000
        nlist = np.array([N]*10*N,dtype=np.uint32)

        rng = fp.GSLrng(101)

        pops = fp.evolve_regions(rng,       #The random number generator
                                 4,         #The number of pops to simulate = number of threads to use.
                                 N,         #Initial population size for each of the 4 demes
                                 nlist[0:], #List of population sizes over time.
                                 0.005,     #Neutral mutation rate (per gamete, per generation)
                                 0.01,      #Deleterious mutation rate (per gamete, per generation)
                                 0.005,     #Recombination rate (per diploid, per generation)
                                 nregions,  #Defined above
                                 sregions,  #Defined above
                                 recregions)#Defined above

        samples = [fp.get_samples(rng,i,100) for i in pops]
        
        details = [fp.get_sample_details(i[1],j) for i,j in zip(samples,pops)]

        ##Get mutation views
        mviews = fp.view_mutations(pops)

        ##For each element in details,
        ##find it in the "mutation view" for that replicate
        for i in range(len(samples)):
            ##Make sure that lengths match up
            self.assertEqual(len(samples[i][1]),len(details[i]))
            for j in range(len(samples[i][1])):
                mm=[X for X in mviews[i] if X['pos'] == samples[i][1][j][0]]
                ##Make sure each position is uniuqe
                self.assertEqual(len(mm),1)
                for mmi in mm:
                    ##Make sure that the position is equal to what we expect
                    self.assertEqual(mmi['pos'],samples[i][1][j][0])
                    ##Make sure selection coefficient matches up
                    self.assertEqual(mmi['s'],details[i].s[j])
                    EP=mmi['n']/float(2000) #"expected" frequency
                    self.assertEqual(EP,details[i].p[j])
            
if __name__ == '__main__':
    unittest.main()
