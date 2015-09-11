import numpy as np
import math

N=1000
N2=5*N
tgrowth=500

#G is the growth rate
G = math.exp( (math.log(N2)-math.log(N))/float(tgrowth) )

nlist = np.array([N]*(10*N+tgrowth),dtype=np.uint32)

#Now, modify the list according to expoential growth rate
for i in range(tgrowth):
    nlist[10*N+i] = round( N*math.pow(G,i+1) )

print (nlist[10*N:])
