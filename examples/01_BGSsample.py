##Sample from the populations simulated in 00_BGS.py

##Run the other script
import example_00_BGS

#Use a list comprehension to get a random sample of size
#n = 20 from each replicate
samples = [fp.get_samples(rng,i,20) for i in pops]
