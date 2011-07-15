#!/usr/bin/env python

#
# This example simulates an admixed population
# 
#
#

from simuOpt import Params, setOptions
setOptions(alleleType='binary', optimized=False, gui=False)
from simuPOP import *


import sys, os, logging
import loadHapMap3, selectMarkers, simuGWAS


def downloadData(logger):
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
    Step 2: Select 2000 markers on a random regions on chromosomes 2, with minor allele frequency 0.05.
    
     > selectMarkers.py --chroms=2 --numMarkers=2000 --startPos=50000000 --filename=ex4_init.pop --minAF=0.05
       --minDist=50000 --HapMap_pops="['HapMap3_JPT+CHB','HapMap3_MKK']" --mergeSubPops=False
    '''
    if os.path.isfile('ex4_init.pop') and os.path.isfile('ex4_init.pop.lst'):
        if logger:
            logger.info('ex4_init.pop already exists. Please remove this file if you would like to regenerate an initial population.')
        return
    if logger:
        logger.info('Select 5000 markers from chromosomes 2')
    pop = selectMarkers.getHapMapMarkers(
        HapMap_dir='HapMap',
        chroms=[2],
        HapMap_pops=['HapMap3_JPT+CHB','HapMap3_MKK'],
        startPos=[50000000],
        minAF=0.05,
        minDist=50000,
        numMarkers=[2000],
        mergeSubPops=False,
        logger=logger)
    if logger:
        logger.info('Saving initial population to ex4_init.pop')
    pop.save('ex4_init.pop')
    if logger:
        logger.info('Saving marker information to ex4_init.pop.lst')
    selectMarkers.saveMarkerList(pop, 'ex4_init.pop.lst', logger)


def expandPop(logger):
    # This is equivalent to
    # > simuGWAS.py --initPop=ex3_init.pop --migrRate=0.0001 --scale=5
    #
    # This just to make this result reproducible.
    getRNG().set(seed=1355)
    #
    filename = 'ex4_expanded.pop'
    if os.path.isfile(filename):
        if logger:
            logger.info('%s already exists. Please remove this file if you would like to regenerate an expanded population.' % filename)
        return
    else:
        if logger:
            logger.info('Simulating an expanded population %s from ex4_init.pop...' % filename)
    pop = loadPopulation('ex4_init.pop')
    pars = Params(simuGWAS.options, initPop=filename, migrRate=0.0001, scale=5)
    pop = simuGWAS.simuGWAS(pars, pop, logger=logger)
    if logger:
        logger.info('Saving expanded population to ' + filename)
    pop.save(filename)        
      

def mix(logger):
    '''Load expanded population and mix using non-random mating'''
    if logger:
        logger.info('Loading population ex4_expanded.pop and mix')
    pop = loadPopulation('ex4_expanded.pop')
    pop.addInfoFields('ancestry')
    # define two virtual subpopulations by ancestry value
    pop.setVirtualSplitter(InfoSplitter(field='ancestry', cutoff = [0.5]))
    # initialize ancestry
    initInfo(pop, [0]*pop.subPopSize(0) + [1]*pop.subPopSize(1), infoFields='ancestry')
    initSex(pop)
    ops=[ MendelianGenoTransmitter(),
          InheritTagger(mode=MEAN, infoFields='ancestry')
        ]
    pop.evolve(
        preOps = Migrator(rate =[
            [0., 0], [0.05, 0]]), 
        matingScheme = HeteroMating(
            matingSchemes=[
                RandomMating(ops=ops),
                RandomMating(subPops=[(0,0)], weight=-0.80, ops=ops),
                RandomMating(subPops=[(0,1)], weight=-0.80, ops=ops)
            ],
        ),
        postOps = PyEval(r"'Generation %d\n' % gen"),
        gen=10,
    )
    # remove the second subpop
    if logger:
        logger.info('Removing MKK subpopulation and save admixed population to ex4_mixed.pop')
    pop.removeSubPops(1)
    pop.save('ex4_mixed.pop')



if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('example4')
    downloadData([2], logger)
    getInitPop(logger)
    expandPop(logger)
    mix(logger)

