
#
# test of sex-linked, mating, recombination etc

from simuPOP import *

pop = population(size=10, loci=[3, 4])
InitByFreq(pop, [.2,.3,.5])
Dump(pop)

simu = simulator(pop, randomMating())

TurnOnDebug(DBG_RECOMBINATOR)
TurnOnDebug(DBG_SIMULATOR)
simu.step(ops=[recombinator(rate=0.2)])
simu.step(ops=[recombinator(rate=0.2, afterLoci=[0,3])])


pop = population(size=10, loci=[3,4], sexChrom=True)
InitByFreq(pop, [.2,.3,.5])
Dump(pop)

simu = simulator(pop, randomMating())

# first, let us see if m_unusedPart is set OK
TurnOnDebug(DBG_RECOMBINATOR)
TurnOnDebug(DBG_SIMULATOR)
simu.step(ops=[recombinator(rate=0.2)])
simu.step(ops=[recombinator(rate=0.2, afterLoci=[0,3])])

#
#


pop = population(size=10, loci=[4])
InitByFreq(pop, [.2,.3,.5])
Dump(pop)

simu = simulator(pop, randomMating())

TurnOnDebug(DBG_RECOMBINATOR)
TurnOnDebug(DBG_SIMULATOR)
simu.step(ops=[recombinator(rate=0.2)])
simu.step(ops=[recombinator(rate=0.2, afterLoci=[0,3])])


pop = population(size=10, loci=[4], sexChrom=True)
InitByFreq(pop, [.2,.3,.5])
Dump(pop)

simu = simulator(pop, randomMating())

# first, let us see if m_unusedPart is set OK
TurnOnDebug(DBG_RECOMBINATOR)
TurnOnDebug(DBG_SIMULATOR)
simu.step(ops=[recombinator(rate=0.2)])
simu.step(ops=[recombinator(rate=0.2, afterLoci=[0,3])])

#
#
# Now, let us look at recombination


from simuPOP import *

pop = population(size=10, loci=[3, 4], sexChrom=False)
InitByFreq(pop, [.2, .3, .5], sex=[Male]*5+[Female]*5 )
#Dump(pop)
# set Y sex chromosome
InitByValue(pop, [1,2,3,4], atLoci=[3,4,5,6],
            indRange=[5,9], atPloidy=1)
#Dump(pop)

# now, recombine
simu = simulator(pop, randomMating())
TurnOnDebug(DBG_RECOMBINATOR)
simu.step( ops=[recombinator(rate=0.3)])
