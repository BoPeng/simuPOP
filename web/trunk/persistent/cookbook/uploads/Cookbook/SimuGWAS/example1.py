#!/usr/bin/env python

#
# This example generates a sample using default parameters. and compares it with the HapMap
# population from which it is generated.
#

import loadHapMap3, selectMarkers, simuGWAS
import os, logging

import simuOpt
simuOpt.setOptions(gui=False, alleleType='binary')
from simuPOP import *

def downloadData(chroms, logger):
    '''
    Download and create populations from the third phase of the HapMap3 data.
    This equivalent to command

    > loadHapMap3.py --chroms=2 --dest=HapMap
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
    Select 5000 markers on chromosome 2p16.3, convering the ENr112 ENCODE
    region. This is equivalent to command:
    
     > selectMarkers.py --chroms=2 --numMarkers=5000 --startPos=51000000 --filename=ex1_init.pop
    
    The resulting population has 993 individuals in a single population. It covers
    region chr2:5102265..13037663. The average marker distance is 1587bp.
    '''
    if os.path.isfile('ex1_init.pop') and os.path.isfile('ex1_init.pop.lst'):
        if logger:
            logger.info('ex1_init.pop already exists. Please remove this file if you would like to regenerate an initial population.')
        return
    if logger:
        logger.info('Select 5000 markers from chromosome2')
    pop = selectMarkers.getHapMapMarkers(
        HapMap_dir='HapMap',
        chroms=[2], 
        HapMap_pops=selectMarkers.HapMap3_pops,
        startPos=51000000,
        numMarkers=5000,
        mergeSubPops=True,
        logger=logger)
    if logger:
        logger.info('Saving initial population to ex1_init.pop')
    pop.save('ex1_init.pop')
    if logger:
        logger.info('Saving marker information to ex1_init.pop.lst')
    selectMarkers.saveMarkerList(pop, 'ex1_init.pop.lst', logger)

def expandPop(logger):
    # Evolve the population to a population of 50,000 individuals, without a disease model.
    #
    # > simuGWAS.py --initPop=ex1_init.pop --DPL=[] --filename=ex1_expanded_1.pop --scale=1 --optimized --dumpRec=rec.log --haploCount='(0,100)'
    # > simuGWAS.py --initPop=ex1_init.pop --DPL=[] --expandSize=25000 --filename=ex1_expanded_2.pop --scale=2 --optimized --dumpRec=rec2.log --haploCount='(0,100)'
    # > simuGWAS.py --initPop=ex1_init.pop --DPL=[] --expandSize=10000 --filename=ex1_expanded_5.pop --scale=5 --optimized --dumpRec=rec5.log --haploCount='(0,100)'
    # > simuGWAS.py --initPop=ex1_init.pop --DPL=[] --expandSize=50000 --filename=ex1_expanded_55.pop --scale=5 --optimized --dumpRec=rec55.log --haploCount='(0,100)'
    #
    from simuPOP import getRNG
    getRNG().set(seed=1355)
    #
    for scale,popSize,name in [
            (1, 50000, '1'),
            (2, 25000, '2'),
            (5,10000, '5'),
            (5,50000, '55')
        ]:
        filename = 'ex1_expanded_%s.pop' % name
        if os.path.isfile(filename):
            if logger:
                logger.info('%s already exists. Please remove this file if you would like to regenerate an expanded population.' % filename)
            continue
        else:
            if logger:
                logger.info('Simulating an expanded population with scaling factor %d and ending population size %d' % (scale, popSize))
        pop = loadPopulation('ex1_init.pop')
        pars = simuOpt.Params(simuGWAS.options, migrRate=0.0001, scale=scale,
            DPL=[], expandSize=popSize,
            dumpRec='rec%d.log' % scale,
            haploCount=(0,100)
        )
        pop = simuGWAS.simuGWAS(pars, pop, logger=logger)
        if logger:
            logger.info('Saving expanded population to %s' % filename)
        pop.save(filename)
      

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('example1')
    downloadData([2], logger)
    getInitPop(logger)
    expandPop(logger)
