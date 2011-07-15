#!/usr/bin/env python

#
# This example simulates recent and remote selective sweep and draw
# trio samples.

import sys, os, math, logging

import simuOpt
simuOpt.setOptions(alleleType='binary', optimized=False, gui=False)
from simuPOP import *

import loadHapMap3, selectMarkers, simuGWAS

def downloadData(chroms, logger):
    '''
    Download and create populations from the third phase of the HapMap3 data.
    This equivalent to command

    > loadHapMap3.py --chroms=2
    '''
    if not os.path.isdir('HapMap'):
        os.mkdir('HapMap')
    for chrom in chroms:
        for popName in loadHapMap3.HapMap3_pops:
            filename = 'HapMap/HapMap3_%s_chr%d.pop' % (popName, chrom)
            if not os.path.isfile(filename):
                pop = loadHapMap3.loadHapMapPop(chrom, popName, logger)
                pop.save(filename)

def getInitPop(logger):
    '''
    Select 500 markers from the JPT+CHB population on an ENCODE region on chromosome 2
    
      > selectMarkers.py --chroms=2 --numMarkers=500 --startPos=234156563 \
        --filename=ex3_init.pop --HapMap_pops=HapMap3_JPT+CHB
    '''
    if os.path.isfile('ex3_init.pop') and os.path.isfile('ex3_init.pop.lst'):
        if logger:
            logger.info('ex3_init.pop already exists. Please remove this file if you would like to regenerate an initial population.')
        return
    if logger:
        logger.info('Select 500 markers from chromosomes 2')
    pop = selectMarkers.getHapMapMarkers(
        HapMap_dir='HapMap',
        chroms=[2],
        HapMap_pops=['HapMap3_JPT+CHB'],
        startPos=[234156563],
        numMarkers=[500],
        logger=logger)
    if logger:
        logger.info('Saving initial population to ex3_init.pop')
    pop.save('ex3_init.pop')
    if logger:
        logger.info('Saving marker information to ex3_init.pop.lst')
    selectMarkers.saveMarkerList(pop, 'ex3_init.pop.lst', logger)


DPL = 'rs2173746'

def shortsweep(logger):
    # this is equivalent to
    #
    # simuGWAS.py --initPop=ex3_init.pop --DPL=rs658054 --filename=ex3_shortsweep.pop \
    #   --curAlleleFreq=1 --trajectory='backward' --scale=1 --fitness="[1,1.05,1.1]"
    #
    # This just to make this result reproducible.
    getRNG().set(seed=1355)
    #
    filename = 'ex3_shortsweep.pop'
    if os.path.isfile(filename):
        if logger:
            logger.info('%s already exists. Please remove this file if you would like to regenerate an expanded population.' % filename)
        return
    else:
        if logger:
            logger.info('Simulating an expanded population %s from ex3_init.pop...' % filename)
    pop = loadPopulation('ex3_init.pop')
    pars = simuOpt.Params(simuGWAS.options, DPL=[DPL],
        curAlleleFreq=[0.99], trajectory='backward',
        trajFile='ex3_shortsweep.traj',
        scale=1, fitness=[1, 1.07, 1.14])
    pop = simuGWAS.simuGWAS(pars, pop, logger=logger)
    if logger:
        logger.info('Saving expanded population to %s' % filename)
    pop.save(filename)
        

def longsweep(logger):
    # this is equivalent to
    #
    # simuGWAS.py --initPop=ex3_init.pop --DPL=rs658054 --filename=ex3_longsweep.pop \
    #   --curAlleleFreq=1 --trajectory='forward' --scale=1 --fitness="[1,1.02,1.03]"
    #
    # This just to make this result reproducible.
    getRNG().set(seed=1345)
    #
    filename = 'ex3_longsweep.pop'
    if os.path.isfile(filename):
        if logger:
            logger.info('%s already exists. Please remove this file if you would like to regenerate an expanded population.' % filename)
        return
    else:
        if logger:
            logger.info('Simulating an expanded population %s from ex3_init.pop...' % filename)
    pop = loadPopulation('ex3_init.pop')
    pars = simuOpt.Params(simuGWAS.options, DPL=[DPL],
        curAlleleFreq=[0.99], trajectory='forward',
        trajFile='ex3_longsweep.traj',
        scale=1, fitness=[1, 1.02, 1.03])
    pop = simuGWAS.simuGWAS(pars, pop, logger=logger)
    if logger:
        logger.info('Saving expanded population to %s' % filename)
    pop.save(filename)
        

# IDs of the parents of selected offspring
parentalIDs = set()
# number of discarded individuals, used to calculate disease prevalence.
discardedInds = 0


alpha = -0.5
beta = -1.0
def _myPenetrance(geno):
    g = geno[0] + geno[1]
    logit = alpha + beta*g
    return math.exp(logit) / (1. + math.exp(logit))
    #return [0.5, 0.4, 0.01][g]

def _selectTrio(off, param):
    'Determine if the offspring can be kept.'
    global discardedInds, parentalIDs
    if off.affected() and not (off.father_id in parentalIDs or off.mother_id in parentalIDs):
        off.setAffected(True)
        parentalIDs |= set([off.father_id, off.mother_id])
        return True
    else:
        discardedInds += 1
        return False

def generateTrioSamples(logger=None):
    '''
    Generate trio samples to be analyzed by HelixTree.
    '''
    for popFile in ['ex3_shortsweep.pop', 'ex3_longsweep.pop']:
        if logger:
            logger.info('Loading population ' + popFile)
        pop = loadPopulation(popFile)
        locus = pop.locusByName(DPL)
        # save a map file
        map = open('ex3.map', 'w')
        print >> map, 'CHROMOSOME MARKER POSITION'
        for loc in range(pop.totNumLoci()):
            print >> map, pop.chromName(0), pop.locusName(loc), pop.locusPos(loc)/1e6
        map.close()
        dat = open('ex3.dat', 'w')
        print >> dat, 'A disease'
        for loc in range(pop.totNumLoci()):
            print >> dat, 'M', pop.locusName(loc)
        dat.close()
        pop.addInfoFields(['ind_id', 'father_id', 'mother_id'])
        # keep parental generation
        pop.setAncestralDepth(1)
        # give everyone an unique ID
        tagID(pop, reset=True)
        stat(pop, genoFreq=locus)
        print 'Genotype frequency: ', pop.dvars().genoFreq[locus]
        pop.evolve(
            preOps = PyPenetrance(func=_myPenetrance, loci=[locus]),
            matingScheme=RandomMating(
                ops=[
                    # pass genotype
                    MendelianGenoTransmitter(),
                    # assign new ID to offspring
                    IdTagger(),
                    # record the parent of each offspring
                    PedigreeTagger(),
                    # determine offspring affection status
                    PyPenetrance(func=_myPenetrance, loci=[locus]),
                    # discard the offspring if it is not affected
                    # or if parents have been chosen
                    PyOperator(func=_selectTrio, param=locus)
                ], subPopSize=1000
            ),
            gen = 1
        )
        global selectedInds, discardedInds, parentalIDs
        sample = pop.extractIndividuals(IDs=list(parentalIDs) + list(pop.indInfo('ind_id')))
        if logger:
            logger.info('Saving sample to file ' + popFile.replace('.pop', '_sample.pop'))
        sample.save(popFile.replace('.pop', '_sample.pop'))
        if logger:
            logger.info('Disease prevalence: %.4f' % (1000. / (1000 + discardedInds)))
        # save file in PED format so I need to change family ID...
        if logger:
            logger.info('Saving sample to file ' + popFile.replace('.pop', '.ped'))
        # write to merlin format
        csv = open(popFile.replace('.pop', '.ped'), 'w')
        #print >> csv, ' '.join(pop.lociNames())
        def genoString(ind):
            alleles = []
            for i in range(ind.totNumLoci()):
                alleles.extend([str(ind.allele(i, 0)+1), str(ind.allele(i, 1)+1)])
            return ' '.join(alleles)
        famid = 1
        id = 1
        for ind in sample.individuals():
            father = sample.indByID(ind.father_id)
            mother = sample.indByID(ind.mother_id)
            print >> csv, famid, id, 0, 0, '1', '2' if father.affected() else '1', genoString(father)
            print >> csv, famid, id+1, 0, 0, '2', '2' if mother.affected() else '1', genoString(mother)
            print >> csv, famid, id+2, id, id+1, '1' if ind.sex() == MALE else '2', '2' if ind.affected() else '1', genoString(ind)
            famid += 1
            id += 3
        csv.close()
        # for next population
        parentalIDs = set()
        selectedInds = 0
        discardedInds = 0

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('example3')
    downloadData([2], logger)
    getInitPop(logger)
    shortsweep(logger)
    longsweep(logger)
    generateTrioSamples(logger=logger)
    
