#!/usr/bin/env python

#
# This example demonstrates how to call simuRareVariants.py to run a simulation,
# assign a quantitative trait model and draw a sample from the population.
#
#

from srv import simuRareVariants, saveMarkerInfoToFile, saveMutantsToFile
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
pop = sim.loadPopulation('example.pop')
# add an information field to every individual
pop.addInfoFields('trait')

def trait(geno):
    '''
    Define a quantitative trait model. geno is the genotype of each individual
    passed by pyQuanTrait. Note that this population is in mutational space so
    geno is the location instead of value of mutants. 0 values should be ignored.
    '''
    # get number of mutants
    count = len(geno) - geno.count(0)
    # N(count, 1) using simuPOP's random number generator.
    return sim.getRNG().randNormal(count, 1)

# apply the quantitative trait model.
sim.pyQuanTrait(pop, func=trait, loci=sim.ALL_AVAIL, infoFields='trait')
# get all trait values
traits = list(pop.indInfo('trait'))
traits.sort()
# get the first and last 500 trait values
cutoff1 = traits[499]
cutoff2 = traits[-500]
# 
print 'Cutoff values are %.4f and %.4f' % (cutoff1, cutoff2)
pop.removeIndividuals(filter=lambda trait: trait > cutoff1 and trait < cutoff2)
print pop.popSize()
markers = saveMarkerInfoToFile(pop, 'sample.map')
print 'Mutant locations are saved to file sample.map'
saveMutantsToFile(pop, 'sample.mut')
print 'Mutants are saved to file sample.mut'



