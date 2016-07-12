import unittest
import fwdpy as fp
import numpy as np
import pandas as pd

##Run some quick sims that we use for tests below:

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
##Get fixation views:
fixations = [fp.view_fixations(i) for i in pops]
##Get mutation views:
mviews = fp.view_mutations(pops)

class test_sample_details(unittest.TestCase):
    def test_CompareUsingView(self):
        """
        Evolve some population, take samples,
        get details of those samples, make sure
        those details match up what 'views' think should
        be in the population.
        """
        samples = [fp.get_samples(rng,i,100) for i in pops]
        
        details = [fp.get_sample_details(i[1],j) for i,j in zip(samples,pops)]



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

class test_Fixations(unittest.TestCase):
    def test_CountFixationsInSample1(self):
        """
        Evolve some pops, take samples w/o removing fixations,
        Make sure that all mutations in a view of fixations
        are also present in the sample.

        Ultimately, this is just a test of whether or not
        fwdpy passes the correct args to fwdpp.
        """
        #False = DO NOT remove fixed variants
        samples = [fp.ms_sample(rng,i,100,False) for i in pops] #This fxn merges all selected + neutral variants into 1 object per replicate
    

        #for each collection of fixations...
        for i in range(len(fixations)):
            #...go thru each fixation...
            FOUND=0
            for f in fixations[i]:
                #...find it in the sample object for the correct replicate...
                mm=[X for X in samples[i] if X[0] == f['pos']]
                FOUND+=len(mm)
            self.assertEqual(FOUND,len(fixations[i]))
    def test_CountFixationsInSample2(self):
        #Repeat above test, but using sampler where neutral + selected variants are in different objects
        samples = [fp.get_samples(rng,i,100,False) for i in pops]

        for i in range(len(fixations)):
            FOUNDN=0
            FOUNDS=0
            for f in fixations[i]:
                if f['neutral'] is True:
                    mm=[X for X in samples[i][0] if X[0]==f['pos']]
                    FOUNDN+=len(mm)
                else:
                    mm=[X for X in samples[i][1] if X[0]==f['pos']]
                    FOUNDS+=len(mm)
                    self.assertEqual(FOUNDN+FOUNDS,len(fixations[i]))
    def test_CountFixationsInSample3(self):
        #Finally, make sure that sample w/o including fixations works,
        #meaning that no mutations are present in the sample with
        #derived count == nsam
        samples = [fp.get_samples(rng,i,100) for i in pops]
        ##Checks that no "true fixations" are in vector
        for i in range(len(fixations)):
            FOUNDN=0
            FOUNDS=0
            for f in fixations[i]:
                if f['neutral'] is True:
                    mm=[X for X in samples[i][0] if X[0]==f['pos']]
                    FOUNDN+=len(mm)
                else:
                    mm=[X for X in samples[i][1] if X[0]==f['pos']]
                    FOUNDS+=len(mm)
                    self.assertEqual(FOUNDN+FOUNDS,0)
        ##Checks that no polymorphisms are present as fixations in sample
        for i in samples:
            for j in i:
                x=j[0][1].count(b'1')
                self.assertTrue(x<100)
                x=j[1][1].count(b'1')
                self.assertTrue(x<100)
                
    def test_CountFixationsInSample4(self):
        samples = [fp.ms_sample(rng,i,100) for i in pops]
        ##Checks that no "true fixations" are in vector
        for i in range(len(fixations)):
            FOUND=0
            for f in fixations[i]:
                mm=[X for X in samples[i] if X[0] == f['pos']]
                FOUND+=len(mm)
                self.assertEqual(FOUND,0)
        ##Checks that no polymorphisms are present as fixations in sample
        for i in samples:
            for j in i:
                self.assertTrue(j[1].count(b'1')<100)
                
if __name__ == '__main__':
    unittest.main()
