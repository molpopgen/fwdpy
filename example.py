import fwdpy
rng = fwdpy.GSLrng(100)
pop = fwdpy.evolve_pops_t(rng,3,1000,[1000]*int(1e4),50,50)
s = fwdpy.ms_sample(rng,pop,10)

###fxn to ask if a site is a singleton
def isSingleton( site ):
    ones=0
    for i in site:
        if i == '1':
            ones += 1
    if ones == 1:
        return True
    return False

##fxn to count derived singletonss in a sample
def countDerived( sample ):
    nsing=0
    for i in range(len(sample)):
        if isSingleton(sample[i][1]):
            nsing += 1
    return nsing

#list comprehension for automatic vectorizing
D = [fwdpy.TajimasD(si) for si in s]
print (D)
nsing = [countDerived(si) for si in s]
print nsing
print("\n\n")
print(s[0])
print("\n\n")
print(s[1])
print("\n\n")
print(s[2])

##Get s,h for each mutation (all irrelevant here, but just an example...)
sh = fwdpy.get_sample_details(s,pop)

##Now, evolve them some more and end with a bottleneck + recent, partial recovery
fwdpy.evolve_pops_more_t(rng,pop,[1000]*int(1e4) + [500]*100 + [750]*10,50,50)
for i in range(len(pop)):
    print pop.generation(i)," ",pop.popsize(i)," ",pop.sane(i)
