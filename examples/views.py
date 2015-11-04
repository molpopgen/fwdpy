### Example of taking "views" from a simulated population
from __future__ import print_function
import fwdpy as fp
import array, pandas as pd

##Run the background selection simulation
from background_selection_setup import *

#Obtain the mutations in each population.
#Each element is a list of dictionaries with mutation information
mutations = [fp.view_mutations(i) for i in pops]

##Look at the raw info for the first mutation in each pop:
for i in mutations:
    print(i[0])

##Let's make the views nicer
