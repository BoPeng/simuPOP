#!/usr/bin/env python
#
# Purpose:
#  testing of random number generators and bernulli trials
#
# Author:
#  Bo Peng (bpeng@rice.edu)
#
from simuPOP import *


# binomial function
def binom(n, m):
  b = [0] * (n+1)
  b[0] = 1
  for i in xrange(1, n+1):
    b[i] = 1
    j = i - 1
    while j > 0:
      b[j] += b[j-1]
      j -= 1
  return b[m]

print listAllRNG()
setRNG()
setRNG("ranlux389")
setRNG("random64-libc5")

# help(RNG)

# get the random number generator.
# currently, there is one global RNG used in the system
rng = rng()
print rng.name()

# return 0,...,n
N=10000
p = .7
n = 100
val = [0]*N
for i in range(0,N):
  val[i] = rng.randBinomial(n, p)

# the theoretical ones
th = [0]*(n+1)
for k in range(0,n+1):
  th[k] = binom(n,k)*p**k*(1-p)**(n-k) *N

# draw a histogram
from scipy import *
xplt.hold('off')
xplt.histogram(val)
xplt.hold('on')
xplt.plot(range(0,n+1), th)
xplt.hold('off')

# return 0,1,2,3,4
for i in range(1,N):
  val[i] =  rng.randInt(5)
xplt.histogram(val,nbins=5, range=[0,5])  
val.count(0)
val.count(1)
  
for i in range(1,N):
  val[i] = rng.randUniform01()
xplt.histogram(val,nbins=10, range=[0,1])


# do bernulli trails
# pointer type trials can not be tested directly
p=0.01
n=1000
bt = BernulliTrials(rng(), [0.01,0.5], n)
bt.doTrial()
print bitSet(bt.trial())
print bt.viewSucc(0)
print bt.viewSucc(1)
