#!/usr/bin/env python

#
# This example demonstrates how to call simuRareVariants.py to run a simulation,
# assign a quantitative trait model and draw a sample from the population.
#
#

from simuRareVariants import simuRareVariants, saveMarkerInfoToFile, saveMutantsToFile
import simuPOP as sim

# call function simuRareVariants to get simulated population. You can
# also call the script from command line.
#
# A mutations.lst file will only be produced if verbose=2
simuRareVariants(regions=['chr1:1..63000'], N=(8100, 7900, 900000),
    numGen=(20000, 10, 370), mu=1.8e-8, numStep=[100,1,10],
    selModel='multiplicative', selDist='constant', selCoef=None,
    verbose=2)

lst = open('mutations.lst')
ages = {}
for line in lst.readlines():
    gen,mut,ind,mtype = map(int, line.split())
    if mtype != 2:
        ages[mut] = gen

print 'mutation generation'
for mut in ages.keys():
    print mut, gen - ages[mut]


