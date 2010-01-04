from simuPOP import *
import os, sys, exceptions

#!/usr/bin/env python
#
# Purpose:
#     This python module provides several utility functions that save a simuPOP
#     Population in SOLAR formats,
#
# License:
#     This program is freely available in the simuPOP's online cookbook
#     (http://simupop.sourceforge.net/cookbook). You can redistribute it and/or
#     modify it freely, with or without this license notice. However, this
#     license notice should be present in the online version of this file.
#
#     This program is NOT part of simuPOP and is NOT protected by simuPOP's GPL
#     license. It is provided in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY. If you notice any bug or have some new ideas, you can
#     modify this file and, as a courtesy to other simuPOP users, incoporate
#     your changes to the online version of this file. If you are uncertain
#     about your changes, please feel free to discuss your changes in the
#     simuPOP mailinglist (simupop-list@lists.sourceforge.net, subscription
#     required).
#
# Change Log:
#     2009-06-25 Bo Peng <bpeng@mdanderson.org>
#
#         Move the function SaveSolarFrqFile from simuUtil.py to the online
#         cookbook.
# 

def saveSolarFrqFile(pop, output='', loci=[], calcFreq=True):
    '''
    Output a frequency file, in a format readable by solar calcFreq
    '''
    if type(pop) == type(''):
        pop = LoadPopulation(pop)
    if output != '':
        file = output
    else:
        raise exceptions.ValueError, "Please specify output"
    # open data file and pedigree file to write.
    try:
        frqOut = open(file, "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file + ".frq to write."
    if loci == []:
        loci = range(0, pop.totNumLoci())
    if calcFreq or not pop.vars().has_key('alleleFreq'):
        stat(pop, alleleFreq=loci)
    alleleFreq = pop.dvars().alleleFreq
    for m in loci:
        try:
            print >> frqOut, pop.locusName(m),
            for a in range(len(alleleFreq[m])):
                print >> frqOut, '\t%d\t%f' % (a+1, alleleFreq[m][a]),
            print >> frqOut
        except:
            print "Can not output allele frequency for marker %s " % m
    frqOut.close()


if __name__ == '__main__':
    #for testing
    outfile = 'SolarOut.dat'
    pop = Population(size=[100]*3, loci=[2,3],
                     lociNames=['loc1-1','loc1-2','loc2-1','loc2-2','loc2-3'])
    initByFreq(pop, [.1, .3, .6])
    saveSolarFrqFile(pop, outfile)
    print open(outfile).read()


