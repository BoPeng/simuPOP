# test sampling issues
#

from simuUtil import *

simu = simulator(
    population(subPop=[100,200], ploidy=2, loci=[5,10],
      ancestralDepth=1, maxAllele=9),
    randomMating(numOffsprings=2))
simu.evolve([
    stat( alleleFreq=[0,1], genoFreq=[0,1]),
    migrator(rate=[[0.1,0.1],[0.1,0.1]]),
    basicPenetrance(locus=0,
      penetrance={'1-1':0,'1-2':.7,'2-2':1}),
    parentsTagger(),
    ],
 preOps=[  initByFreq(alleleFreq=[.2,.8], atLoci=[0]),
    initByFreq(alleleFreq=[.2]*5, atLoci=range(1, simu.totNumLoci()))   ],
 end=4
)
pop = simu.getPopulation(0)
Dump(pop, ancestralPops=True)

# random sampling
s = RandomSample(pop, 10)
Dump(s[0])

s = RandomSample(pop,[2,8])
Dump(s[0])

# case control sampling.
s = CaseControlSample(pop, 10,10)
Dump(s[0])
listVars(s[0].dvars())

s =  CaseControlSample(pop, [5,5],[5,5])
Dump(s[0])
listVars(s[0].dvars())

# find sibpairs
s=AffectedSibpairSample(pop)
Dump(s[0])
listVars(pop.dvars())

# from each subpop
s = AffectedSibpairSample(pop,10)
Dump(s[0], ancestralPops=True)
listVars(s[0].dvars())
s = AffectedSibpairSample(pop,[5,5])
Dump(s[0])
listVars(s[0].dvars())
# 
Dump(s[0], ancestralPops=True)

# save to linkage format
# SaveLinkage(s[0], output='s0')

SaveRandFam(s[0], output='s0')
