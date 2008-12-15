#!/usr/bin/env python

'''
File: Mating_haplodiploid.py
Author: Bo Peng (bpeng@mdanderson.org)

Purpose:
  This script demonstrates how to handle recombination in haplodiploid populations.

$Date: 2008-12-14 02:19:13 -0500 (Sun, 14 Dec 2008) $
$Revision: 2114 $
$HeadURL: https://simupop.svn.sourceforge.net/svnroot/simupop/trunk/doc/cookbook/Mating_haplodiploid.py $
'''

from simuPOP import *

#
# This script demonstrate how to construct a genotype transmitter (during
# mating operator) that can be used to handle recombination in a
# haplodiploid population.
#
# Recall that, in a haplodiploid population
# - Female individuals have two sets of homologous chromosomes.
# - Male individuals have one set of homologous chromosomes.
# - Female offspring get one set of chromosome from her mother (possibly
#   with recombination), and one set of chromosome from her father (no
#   recombination).
# - Male offspring et on set of chromosome from his mother (possibly
#   with recombination).
#

class haplodiploidRecombinator(pyOperator):
    def __init__(self, intensity, rate, loci, convProb, convMode, convParam,
            *args, **kwargs):
        #
        pyOperator.__init__(self, func=self.transmitGenotype, param=None,
            stage=DuringMating, formOffGenotype=True, passOffspringOnly=False)
        self.recombinator = recombinator(intensity, rate, loci,
            convProb, convMode, convParam)
        self.initialized = False
        self.copier = cloneGenoTransmitter()
        # 

    def transmitGenotype(self, pop, off, dad, mom):
        if not self.initialized:
            self.recombinator.initialize(pop)
            self.copier.initialize(pop)
            initialized = True
        self.recombinator.transmitGenotype(mom, off, 0)
        if offspring.sex() == Female:
            self.copier.transmitGenotype(dad, 0, off, 1)
        return True


def haplodiploidRecMating(replacement=True, intensity=-1, rate=[], loci=[],
        convProb=0, convMode=CONVERT_NumMarkers, convParam=1, 
		numOffspring = 1., numOffspringFunc = None, numOffspringParam= 1, mode = MATE_NumOffspring,
		sexParam = 0.5, sexMode = MATE_RandomSex, ops = [], newSubPopSize = [],
		newSubPopSizeFunc = None, newSubPopSizeExpr = "", 
		subPop = (), weight = 0):
    return pyMating(
        chooser = randomParentsChooser(replacement),
        generator = offspringGenerator(
            [haplodiploidRecombinator(intensity, rate, loci, convProb,
                convMode, convParam)], 2, 
            numOffspring, numOffspringFunc,
            numOffspringParam, mode, sexParam, sexMode),
        newSubPopSizeExpr = newSubPopSizeExpr,
        newSubPopSizeFunc = newSubPopSizeFunc,
        subPop = subPop,
        weight = weight)


def simuHaplodiploid(N, numMito=3, gen=10):
    '''
    '''
    pop = population(N, 
        # record indexes of parents for verification purpose
        ancGen=1, infoFields=['father_idx', 'mother_idx'])

    simu = simulator(pop, haplodiploidRecMating(rate=0.01))

    simu.evolve(
        preOps=[
            # initialize alleles 0, 1, 2, 3 with different frequencies
            initByFreq([0.4] + [0.2]*3),
        ],
        ops=[
            parentsTagger(),
            dumper(structure=False),
        ],
        gen = gen
    )
    return simu.extract(0)


if __name__ == '__main__':
    simuHaplodiploid(10, 3, 2)

# A possible output:

