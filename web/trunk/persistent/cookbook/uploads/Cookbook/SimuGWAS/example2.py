#!/usr/bin/env python

#
# This example simulates a case-control sample using a GxE model.
#

import sys, os, random, math, logging

import simuOpt
simuOpt.setOptions(gui=False, alleleType='binary', optimized=True)
from simuPOP import *

import loadHapMap3, selectMarkers, simuGWAS

def downloadData(chroms, logger):
    '''
    Download and create populations from the third phase of the HapMap3 data.
    This equivalent to command

    > loadHapMap3.py --chroms='[2, 5,10]' --dest=HapMap
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
    Select 2000 markers on a random regions on chromosomes 2, 5 and 10, using markers from the Illumina 1M chipset
    
    From command line, you could prepare a marker list file from illumina annotation file
         > cut -d, -f2 HumanHap550v3_A.csv  > HumanHap550v3_A.lst
    and then select markers
         > selectMarkers.py --markerList='HumanHap550v3_A.lst' --chroms='[2, 5,10]' \
            --numMarkers='[2000,2000,2000]' --startPos='[20000000, 20000000,40000000]' \
            --filename=ex2_init.pop
    '''
    if os.path.isfile('ex2_init.pop') and os.path.isfile('ex2_init.pop.lst'):
        if logger:
            logger.info('ex2_init.pop already exists. Please remove this file if you would like to regenerate an initial population.')
        return
    ann = open('HumanHap550v3_A.lst')
    names = []
    for line in ann:
        names.append(line.split(',')[1])
    if logger:
        logger.info('Select 6000 markers from chromosomes 2, 5 and 10')
    pop = selectMarkers.getHapMapMarkers(
        names=names,
        HapMap_dir='HapMap',
        chroms=[2, 5, 10],
        HapMap_pops=selectMarkers.HapMap3_pops,
        startPos=[20000000, 20000000, 40000000],
        numMarkers=[2000, 2000, 2000],
        mergeSubPops=True,
        logger=logger)
    if logger:
        logger.info('Saving initial population to ex2_init.pop')
    pop.save('ex2_init.pop')
    if logger:
        logger.info('Saving marker information to ex2_init.pop.lst')
    selectMarkers.saveMarkerList(pop, 'ex2_init.pop.lst', logger)


# rs4491689	2	26494285	A	G	0.2825
# rs6869003	5	27397573	C	T	0.0710
DPL=["rs4491689", "rs6869003"]

def expandPop(logger):
    # This is equivalent to
    #
    #  > simuGWAS.py --initPop=ex2_init.pop --DPL='["rs4491689","rs6869003"]' \
    #  --filename=ex2_expanded.pop --curAlleleFreq='[0.05, 0.15]' --trajectory='forward'  --mlSelModel='multiplicative' \
    #  --scale=1 --optimized --gui=False --fitness='[1, 0.996, 0.994, 1, 1.001, 1.005]'
    #
    # This just to make this result reproducible.
    getRNG().set(seed=1355)
    #
    filename = 'ex2_expanded.pop'
    if os.path.isfile(filename):
        if logger:
            logger.info('%s already exists. Please remove this file if you would like to regenerate an expanded population.' % filename)
        return
    else:
        if logger:
            logger.info('Simulating an expanded population from ex2_init.pop...')
    pop = loadPopulation('ex2_init.pop')
    pars = simuOpt.Params(simuGWAS.options, DPL=DPL,
        curAlleleFreq=[0.05, 0.15], trajectory='forward', 
        trajPlot='ex2_traj.pdf', mlSelModel='multiplicative',
        scale=1, fitness=[1, 0.996, 0.994, 1, 1.001, 1.005])
    pop = simuGWAS.simuGWAS(pars, pop, logger=logger)
    if logger:
        logger.info('Saving expanded population to ex2_expanded.pop')
    pop.save('ex2_expanded.pop')        


selectedCase = 0
selectedControl = 0
discardedInds = 0

alpha = -5
beta1 = 0.20
beta2 = 0.4
beta3 = 0.4
gamma1 = 0.2
gamma2 = 0.4

random.seed(235234)  # to keep result reproducible.
getRNG().set(seed=12345)

def _selectInds(off, param):
    'Determine if the offspring can be kept.'
    e = random.randint(0, 1)
    g1 = off.allele(param[0], 0) + off.allele(param[0], 1)
    g2 = off.allele(param[1], 0) + off.allele(param[1], 1)
    logit = alpha + beta1*g1 + beta2*g2 + beta3*g1*g2 + gamma1*e*g1 + gamma2*e*g2
    affected = random.random() < (1 / (1. + math.exp(-logit)))
    global selectedCase, selectedControl, discardedInds
    if affected:
        if selectedCase < 1000:
            off.setAffected(True)
            selectedCase += 1
            return True
    else:
        if selectedControl < 1000:
            selectedControl += 1
            off.setAffected(False)
            return True
    discardedInds += 1
    return False

def generateSample(numCase=1000, numCtrl=1000, logger=None):
    if logger:
        logger.info('Generating %d cases and %d controls...' % (numCase, numCtrl))
    pop = loadPopulation('ex2_expanded.pop')
    loci = pop.lociByNames(DPL)
    pop.evolve(
        matingScheme=RandomMating(
            ops=[
                MendelianGenoTransmitter(),
                # an individual will be discarded if _selectInds returns False
                PyOperator(func=_selectInds, param=loci)
            ], 
            subPopSize=2000
        ),
        gen = 1
    )
    pop.save('ex2_sample.pop')
    if logger:
        logger.info('Disease prevalence: %.4f' % (1000. / (2000 + discardedInds)))


if __name__ == '__main__':
    # You can change logging level to DEBUG to get more information
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('example2')
    downloadData([2, 5, 10], logger)
    getInitPop(logger)
    expandPop(logger)
    generateSample(logger=logger)
