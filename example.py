import fwdpy
rng = fwdpy.GSLrng(100)
pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*int(1e4),50,50)
s = fwdpy.ms_sample(rng,pop,10)

#list comprehension for automatic vectorizing
D = [fwdpy.TajimasD(si) for si in s]
print (D)
print("\n\n")
print(s[0])
print("\n\n")
print(s[1])
print("\n\n")
print(s[2])

##Now, evolve them some more and end with a bottleneck + recent, partial recovery
fwdpy.evolve_pops_more_t(rng,pop,[1000]*int(1e4) + [500]*100 + [750]*10,50,50)
for i in range(len(pop)):
    print pop.generation(i)," ",pop.popsize(i)," ",pop.sane(i)
