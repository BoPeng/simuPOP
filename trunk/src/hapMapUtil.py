#!/usr/bin/env python

'''
Utility functions to manipulate HapMap data
'''

import sys, os, math

from simuOpt import setOptions
setOptions(alleleType='binary')
from simuPOP import *

# choose some markers to work with
# 
# Method 1: give a list of marker names, get a population with these
# markers
# 

def getMarkersFromName(hapmap_dir, names, chroms=[]):
    ''' get population from marker names
    hapmap_dir: where hapmap data in simuPOP format is stored. The files
        should have been prepared by scripts/loadHapMap.py.

    names: names of markers

    chroms: a list of chromosomes to look in. If empty, all 22 autosomes
        will be tried.

    return: a tuple with
        1. a population with found markers. This population has three
        subpopulations, corresponding to CEU, YRI and JPT+CHB.
        2. names of markers that can not be located in the HapMap data
    '''
    pops = []
    if len(chroms) == 0:
        chroms = range(1, 23)
    for i in chroms:
        markers = []
        print "Loading hapmap chromosome %d..." % i
        pop = LoadPopulation(os.path.join(hapmap_dir, 'hapmap_%d.bin' % i))
        for name in names:
            try:
                idx = pop.locusByName(name)
                print "Locus %s is found. Idx: %d" % (name, idx)
                markers.append(idx)
            except:
                pass
        if len(markers) > 0:
            markers.sort()
            pop.removeLoci(keep=markers)
            pops.append(pop)
        else:
            del pop
    ret = MergePopulationsByLoci(pops)
    # find out which markers are not found
    if ret.totNumLoci == len(names):
        return (ret, [])
    else:
        foundNames = ret.lociNames()
        notFound = [x for x in names if x not in foundNames]
        return (ret, notFound)


# Method 2: give a range, maximum number of markers, minimal minor allele frequency,
# and minimal distance between adjacent markers
#
# For example,
# 
# get 5000 markers
#    getMarkersFromRange(2, 0, sys.maxint, 5000, 0, 0)
# get 5000 markers with minimal allele frequency 0.1
#    getMarkersFromRange(2, 0, sys.maxint, 5000, 0.1, 0)
# get markers in range 0, 100cM, with minimal distance 0.1cM
#    getMarkersFromRange(2, 0, 100, sys.maxint. 0, 0.1)

def getMarkersFromRange(chrom, startPos, endPos, maxNum, minAF, minDist):
    '''get markers from given range
        chrom:    chromosome number (1-based index)
        startPos: starting position (in cM)
        endPos:   ending position (in cM)
        maxNum:   maximum number of markers to get
        minAF:    minimal minor allele frequency
        minDist:  minimal distance between two markers, in cM
    '''
    print "Loading hapmap population hapmap_%d.bin" % (chrom-1)
    pop = LoadPopulation(os.path.join(hapmap_dir, 'hapmap_%d.bin' % (chrom-1)))
    markers = []
    lastPos = 0
    for loc in range(pop.totNumLoci()):
        pos = pop.locusPos(loc)
        if pos < startPos or pos > endPos:
            continue
        if len(markers) > maxNum:
            break
        if pos - lastPos < minDist:
            continue
        maf = min(pop.dvars().alleleFreq[loc][0], 1 - pop.dvars().alleleFreq[loc][1])
        if maf < minAF:
            continue
        # this marker is OK
        markers.append(loc)
        lastPos = pos
    pop.removeLoci(keep=markers)
    return pop


###########################################################
#
# Evole the hapmap population
#
# NOTE1: use uniform recombination rate for now
# 
# NOTE2: resultinf population will have strong signature
#   of population expansion. The sideeffect is unknown for now.
# 
###########################################################


def evolve(pop, initMultiple=5, endingSize=1e4, expand='linear',
    mergeAt=90, endGen=100, recIntensity=0.01, mutRate=1e-7,
    step=10, keepParents=False, numOffspring=1):
    ''' evolve and expand the hapmap population
    gen: total evolution generation
    initMultiple: copy each individual initMultiple times, to avoid
        rapid loss of genotype variation when population size is small.
    endingSize: ending poplation size
    expand: expanding method, can be linear or exponential
    mergeAt: when to merge population?
    endGen: endingGeneration
    recIntensity: recombination intensity
    mutRate: mutation rate
    step: step at which to display statistics
    keepParents: whether or not keep parental generations
    numOffspring: number of offspring at the last generation
    '''
    print "Starting population size is ", pop.subPopSizes()
    if initMultiple > 1:
        print 'Propagating population to size %s' % [x*initMultiple for x in pop.subPopSizes()]
        pop.resize([x*initMultiple for x in pop.subPopSizes()], 
            propagate=True)
    #
    N0 = pop.popSize()
    N1 = endingSize
    if expand == 'linear':
        rate = (N1-N0)*1.0/(endGen+1)
    else:
        rate = math.exp(math.log(N1*1.0/N0)/(endGen+1))
    def popSizeFunc(gen, cur):
        if expand == 'linear':
            return [int(x+rate/len(cur)) for x in cur]
        else:
            return [int(x*rate) for x in cur]
    #
    simu = simulator(pop, randomMating(newSubPopSizeFunc=popSizeFunc), rep=1)
    operators = [
        mergeSubPops(subPops=[0,1,2], removeEmptySubPops=True, at=[mergeAt]),
        recombinator(intensity=recIntensity),
        kamMutator(rate=mutRate, loci=range(pop.totNumLoci())),
        # estimating Fst from all loci will be slow, so using the first 100 markers
        #stat(popSize=True, Fst=range(min(100, pop.totNumLoci())), step=step),
        stat(popSize=True, step=step),
        pyEval(r'"gen=%d, size=%s\n" % (gen, subPopSize)', step=step)
    ]
    del pop
    simu.evolve(ops=operators, end=endGen)
    if keepParents:
        print "Preparing the last generation"
        simu.addInfoFields(['father_idx', 'mother_idx'])
        simu.setAncestralDepth(1)
        simu.setMatingScheme(randomMating(numOffspring=numOffspring))
        simu.step(ops=operators + [parentsTagger()])
    return simu.getPopulation(0, True)


def sample(pop, DSL, pene, name, size):
    '''sample from the final population
    DSL: disease loci (two locus)
    pene: penetrance value, assuming a two-locus model
         penetrance
                 BB Bb bb
             AA  0  1  2
             Aa  3  4  5
             aa  6  7  8
    name: name to save sample
    size: sample size
    '''
    Stat(pop, alleleFreq=DSL)
    # assume the more common allele is wildtype
    if pop.dvars().alleleFreq[DSL[0]][0] < 0.5:
        da = 0
        wa = 1
    else:
        da = 1
        wa = 0
    if pop.dvars().alleleFreq[DSL[1]][0] < 0.5:
        db = 0
        wb = 1
    else:
        db = 1
        wb = 0
    print "Disease alleles at the DSL are %d, %d" % (da, db)
    # applying penetrance
    def peneFunc(geno):
        # 
        # 
        if da == 0:
            aa = 2 - (geno[0] + geno[1])
        else:
            aa = geno[0] + geno[1]
        if db == 0:
            bb = 2 - (geno[2] + geno[3])
        else:
            bb = geno[2] + geno[3]
        return pene[aa*3+bb]
    PyPenetrance(pop, loci=DSL, func=peneFunc)
    # draw sample
    (sample,) = CaseControlSample(pop, cases=size/2, controls=size/2)
    def comb(geno):
        return geno[0]+geno[1]
    print "Saving sample to %s " % os.path.join(name, 'sample')
    if not os.path.isdir(name):
        os.mkdir(name)
    SaveQTDT(sample, output=os.path.join(name, 'sample'), affectionCode=['0', '1'], fields=['affection'],
        combine=comb)



