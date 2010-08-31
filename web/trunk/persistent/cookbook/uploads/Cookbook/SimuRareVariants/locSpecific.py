#!/usr/bin/env python

#
# This example demonstrates how to pass a location-specific fitness distribution to
# simuRareVariants.py

from simuRareVariants import simuRareVariants
import simuPOP as sim


def myDist(loc):
    '''Define a function that return fitness parameter of newly arising
    mutants.
    '''
    if loc % 3 != 2:
        return - sim.getRNG().randGamma(0.184, 0.160)
    else:
        # the last nucleotide of a condo is not selected
        return 0


# call function simuRareVariants to get simulated population. You can
# also call the script from command line.
#
# We only need pop file to further analysis. No .ped file is saved.
simuRareVariants(regions=['chr1:1..63000'], N=(8100, 7900, 900000),
    numGen=(20000, 10, 370), mu=1.8e-8, numStep=[100,1,10],
    selModel='multiplicative', selDist=myDist, selCoef=None,
    markerFile='locSpecific.map')
