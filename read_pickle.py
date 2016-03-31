import pickle,gzip

PIK="qtrait_pickle.dat"
f=gzip.open(PIK,"rb");

x=pickle.load(f)
for i in range(x):
    j=pickle.load(f)
    print len(j)
    print j[0]
    print j[1]
    print j[len(j)-1]
