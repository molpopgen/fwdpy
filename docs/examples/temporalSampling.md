
# Temporal sampling
_fwdpy_ allows various things to be recorded over time during a simulation.  A family of objects referred to as "temporal samplers" perform these tasks.  All such objects are derived from the base class :class:`fwdpy.fwdpy.TemporalSampler`.
## Sampling nothing
Doing nothing is useful for evolving a population to equilibrium.  The relevant class is :class:`fwdpy.fwdpy.NothingSampler`.

For convenience, :func:`fwdpy.fwdpy.evolve_regions` and :func:`fwdpy.fwdpy.evolve_regions_more` and :func:`fwdpy.fwdpy.evolve_regions_fitness` all implicitly use :class:`fwdpy.fwdpy.NothingSampler`.
Let's evolve 40 populations to mutation-drift equilibrium:


```python
import fwdpy as fp
import numpy as np
import pandas as pd
nregions=[fp.Region(0,1,1)]
sregions=[fp.GammaS(0,1,0.1,0.1,0.1,1.0),
          fp.GammaS(0,1,0.9,-0.2,9.0,0.0)
         ]
recregions=nregions
N=1000
nlist=np.array([N]*(10*N),dtype=np.uint32)
mutrate_neutral=50.0/float(4*N)
recrate=mutrate_neutral
mutrate_sel=mutrate_neutral*0.2
rng=fp.GSLrng(101)
pops=fp.SpopVec(40,1000)
sampler=fp.NothingSampler(len(pops))
#This function implicitly uses a "nothing sampler"
fp.evolve_regions_sampler(rng,pops,sampler,nlist,
                          mutrate_neutral,
                          0.0,   #No selected mutations....
                          recrate,
                          nregions,sregions,recregions,
                          #Only sample every 10N generations,
                          #which is fine b/c we're not sampling anything
                          10*N)
```

## Take samples from population
Example using :class:`fwdpy.fwdpy.PopSampler`

```python
#Take sample of size n=20
sampler=fp.PopSampler(len(pops),20,rng)
fp.evolve_regions_sampler(rng,pops,sampler,
                          nlist[:N], #Evolve for N generations
                            mutrate_neutral,
                            mutrate_sel,   
                              recrate,
                              nregions,sregions,recregions,
                            #Sampler every 100 generations
                              100)
```

The output from this particular sampler type is a generator.  Let's look at the first element of the first sample:


```python
data=sampler[0]
```


```python
print data[0]
```

    (([(0.07927352748811245, '00000000000000001000'), (0.11939710704609752, '00100000000000000000'), (0.1569378236308694, '00000000100000000000'), (0.19993915781378746, '00001000000000000000'), (0.3642488766927272, '00000000000000001000'), (0.3924784071277827, '10000000000000000000'), (0.4982362166047096, '00001000000000000000'), (0.5306595861911774, '00000001000000000000'), (0.6011973915155977, '00000000000000000100'), (0.6021612668409944, '10000000000000000000'), (0.7797581860795617, '00000000001000000000'), (0.810913129709661, '00000000000100000000'), (0.8996576184872538, '00000000000000000100'), (0.9164007280487567, '00000000000000010000')], [(0.02243150118738413, '00000000001000000000'), (0.8127563807647675, '00001000000000000000'), (0.8615972911939025, '00001000000000000000')]), {'origin': [10090, 10018, 10076], 'generation': [10100, 10100, 10100], 'h': [0.0, 0.0, 0.0], 'locus': [0, 0, 0], 'p': [0.0035, 0.0475, 0.007], 's': [-0.31225451236724366, -0.12004595381424654, -0.13281077444110803], 'ftime': [4294967295, 4294967295, 4294967295], 'dcount': [1, 1, 1], 'label': [0, 0, 0]})


These "genotypes" blocks can be used to caculate summary statistics. See the example on using [pylibseq](http://molpopgen.github.io/pylibseq/) for that task.


```python
print data[1]
```

    (([(0.04461725312285125, '00000100000000000000'), (0.41906465985812247, '00000000010000000000'), (0.5870706243440509, '00000000000000100000'), (0.7365321817342192, '00000000000001000000'), (0.8885768735781312, '00010000000000000000')], []), {'origin': [4294967295], 'generation': [10200], 'h': [nan], 'locus': [0], 'p': [nan], 's': [nan], 'ftime': [4294967295], 'dcount': [4294967295], 'label': [65535]})


Each element in d[0] is a tuple:


```python
#The first element are the genotypes
data[0][0]
```




    ([(0.07927352748811245, '00000000000000001000'),
      (0.11939710704609752, '00100000000000000000'),
      (0.1569378236308694, '00000000100000000000'),
      (0.19993915781378746, '00001000000000000000'),
      (0.3642488766927272, '00000000000000001000'),
      (0.3924784071277827, '10000000000000000000'),
      (0.4982362166047096, '00001000000000000000'),
      (0.5306595861911774, '00000001000000000000'),
      (0.6011973915155977, '00000000000000000100'),
      (0.6021612668409944, '10000000000000000000'),
      (0.7797581860795617, '00000000001000000000'),
      (0.810913129709661, '00000000000100000000'),
      (0.8996576184872538, '00000000000000000100'),
      (0.9164007280487567, '00000000000000010000')],
     [(0.02243150118738413, '00000000001000000000'),
      (0.8127563807647675, '00001000000000000000'),
      (0.8615972911939025, '00001000000000000000')])




```python
#The first element in the genotypes are the neutral variants.
#The first value is the position.  The second value is a string
#of genotypes for chromosomes 1 through n.  0 = ancestral/1=derived
data[0][0][0]
```




    [(0.07927352748811245, '00000000000000001000'),
     (0.11939710704609752, '00100000000000000000'),
     (0.1569378236308694, '00000000100000000000'),
     (0.19993915781378746, '00001000000000000000'),
     (0.3642488766927272, '00000000000000001000'),
     (0.3924784071277827, '10000000000000000000'),
     (0.4982362166047096, '00001000000000000000'),
     (0.5306595861911774, '00000001000000000000'),
     (0.6011973915155977, '00000000000000000100'),
     (0.6021612668409944, '10000000000000000000'),
     (0.7797581860795617, '00000000001000000000'),
     (0.810913129709661, '00000000000100000000'),
     (0.8996576184872538, '00000000000000000100'),
     (0.9164007280487567, '00000000000000010000')]




```python
#Same format for selected variants
data[0][0][1]
```




    [(0.02243150118738413, '00000000001000000000'),
     (0.8127563807647675, '00001000000000000000'),
     (0.8615972911939025, '00001000000000000000')]




```python
#This is a dict relating to info re:
#the selected variants.
#dcount = derived freq in sample
#ftime = fixation time. 2^32-1 = has not fixed
#generation = generation when sampling occurred
#h = dominance
#origin = generation when mutation entered population
#p = population frequency
#s = effect size/selection coefficient
data[0][1]
```




    {'dcount': [1, 1, 1],
     'ftime': [4294967295, 4294967295, 4294967295],
     'generation': [10100, 10100, 10100],
     'h': [0.0, 0.0, 0.0],
     'label': [0, 0, 0],
     'locus': [0, 0, 0],
     'origin': [10090, 10018, 10076],
     'p': [0.0035, 0.0475, 0.007],
     's': [-0.31225451236724366, -0.12004595381424654, -0.13281077444110803]}




```python

```

## Tracking mutation frequencies

See the example on tracking mutation frequencies.
The relevant class is :class:`fwdpy.fwdpy.FreqSampler`.
See the example on fixation times for the use of this sampler
