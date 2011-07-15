#!/usr/bin/env python

#
# This example demonstrates how to call simuRareVariants.py to run a simulation,
# generate large pedigrees and save mutants of pedigrees to files.
#
#

from srv import simuRareVariants, saveMarkerInfoToFile, saveMutantsToFile
from simuPOP import sampling
import simuPOP as sim

# call function simuRareVariants to get simulated population. You can
# also call the script from command line.
#
# We only need pop file to further analysis. No .ped file is saved.
simuRareVariants(regions=['chr1:1..63000'], N=(8100, 8100, 7900, 900000),
    G=(20000, 10, 370), mu=1.8e-8, steps=[100,1,10],
    selModel='multiplicative', selDist='constant', selCoef=None,
    popFile='example.pop')

# load population
print 'Loading population example.pop'
pop = sim.loadPopulation('example.pop')

# evolve the population for a few more generations to produce pedigrees
print 'Evolving the population for three generations.'
pop.addInfoFields(['ind_id', 'father_id', 'mother_id'])
sim.tagID(pop)
# save all ancestral generations during evolution
pop.setAncestralDepth(-1)   
pop.evolve(
    matingScheme=sim.RandomMating(numOffspring=(sim.UNIFORM_DISTRIBUTION, 2, 4),
        ops=(sim.MendelianGenoTransmitter(), sim.IdTagger(), sim.PedigreeTagger())),
    gen=3
)
# what is the average number of mutants in this population?
avgMutants = (pop.popSize()*pop.totNumLoci()*2. - pop.genotype().count(0)) / pop.popSize()
print 'Average number of mutants is %.2f' % avgMutants
#
# This contains marker information for the initial population
print 'Mutant locations are saved to file sample.map'
markers = saveMarkerInfoToFile(pop, 'pedigree.map')
#
def myPenet(geno):
    # count the number of mutants
    cnt = len(geno) - geno.count(0)
    # penetrance is 0.5*number of mutants
    return min(1, 0.5 * cnt)

# apply a penetrance model ...
sim.pyPenetrance(pop, func=myPenet, loci=sim.ALL_AVAIL)
# number of affected individuals?
for gen in range(2, -1, -1):
    pop.useAncestralGen(gen)
    sim.stat(pop, numOfAffected=True, vars=('numOfAffected', 'propOfAffected'))
    print 'Generation %d has %d (%.2f%%) affected individuals.' % (gen, pop.dvars().numOfAffected, pop.dvars().propOfAffected)

# draw three generation pedigrees with at least two affected members
print 'Drawing 10 pedigrees with at least three affected members'
sample = sampling.drawThreeGenFamilySample(pop, families=10, numOffspring=(2,3), pedSize=(5, 20), numOfAffected=(3, 10))
#
# ind_id, father_id and mother_id are output, but you also need to separate pedigrees.
saveMutantsToFile(sample, 'pedigree.mut', infoFields=['ind_id', 'father_id', 'mother_id'])
print 'Mutants are saved to file pedigree.mut'




