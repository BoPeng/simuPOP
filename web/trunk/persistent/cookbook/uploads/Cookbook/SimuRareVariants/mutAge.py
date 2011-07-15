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
# A mutations.lst file will only be produced if verbose=2
simuRareVariants(regions=['chr1:1..63000'], N=(8100, 8100, 7900, 90000),
    G=(2000, 10, 370), mu=1.8e-8, steps=[100,1,10],
    mutationModel='infinite_sites', statFile='mutAge.out',
    selModel='multiplicative', selDist='constant', selCoef=None,
    verbose=2)

lst = open('mutations.lst')
ages = {}
for line in lst.readlines():
    gen,mut,ind,mtype = map(int, line.split())
    # 0 valid, 2 relocated, 3 ignored
    if mtype != 3:
        ages[mut] = gen

print 'mutation generation'
for mut in ages.keys():
    print mut, gen - ages[mut]


