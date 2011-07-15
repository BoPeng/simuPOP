#!/usr/bin/env python

# This example demonstrates how to pass an arbitrary distribution to
# simuRareVariants.py

from srv import simuRareVariants
import simuPOP as sim


def myDist():
    '''Define a function that return fitness parameter of newly arising
    mutants.
    '''
    rng = sim.getRNG()
    # positive selection
    if rng.randUniform() < 0.2:
        return - rng.randUniform()/1000.
    # purifying selection
    else:
        return rng.randGamma(0.184, 0.160)


# call function simuRareVariants to get simulated population. You can
# also call the script from command line.
#
# We only need pop file to further analysis. No .ped file is saved.
simuRareVariants(regions=['chr1:1..63000'], N=(8100, 8100, 7900, 90000),
    G=(2000, 10, 370), mu=1.8e-8, steps=[100,1,10],
    selModel='multiplicative', selDist=myDist, selCoef=None,
    markerFile='myDist.map')
