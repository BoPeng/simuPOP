#!/usr/bin/env python
############################################################################
#    Copyright (C) 2004 by Bo Peng                                         
#    bpeng@mdanderson.org
#                                                                          
#    $LastChangedDate: 2007-04-13 15:55:29 -0500 (Fri, 13 Apr 2007) $          
#    $Rev: 909 $                       
#                                                                          
#    This program is free software; you can redistribute it and/or modify  
#    it under the terms of the GNU General Public License as published by  
#    the Free Software Foundation; either version 2 of the License, or     
#    (at your option) any later version.                                   
#                                                                                                                                                    
#    This program is distributed in the hope that it will be useful,             
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the                 
#    GNU General Public License for more details.                                                    
#                                                                                                                                                    
#    You should havereceived a copy of the GNU General Public License         
#    along with this program; if not, write to the                                                 
#    Free Software Foundation, Inc.,                                                                             
#    59 Temple Place - Suite 330, Boston, MA    02111-1307, USA.                         
############################################################################

'''
Utility functions to manipulate HapMap data. These functions 
are provided as samples on how to load and evolve the HapMap
dataset. They tend to change frequently so do not call
these functions directly. It is recommended that you
copy these function to your script when you need to use 
them.

'''

import sys, os, math

from simuOpt import setOptions
setOptions(alleleType='binary')
from simuPOP import *
from simuUtil import SaveQTDT

# choose some markers to work with
# 
# Method 1: give a list of marker names, get a population with these
# markers
# 

def getMarkersFromName(hapmap_dir, names, chroms=[]):
    ''' Get population from marker names. This function 
        returns a tuple with a population with found markers and names of
        markers that can not be located in the HapMap data. The returned
        population has three subpopulations, corresponding to CEU, YRI and
        JPT+CHB hapmap populations.
    
    hapmap_dir: where hapmap data in simuPOP format is stored. The files
        should have been prepared by scripts/loadHapMap.py.

    names: names of markers

    chroms: a list of chromosomes to look in. If empty, all 22 autosomes
        will be tried.

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


def getMarkersFromRange(HapMap_dir, pops, chrom, startPos, endPos, maxNum, 
    minAF, minDiffAF, minDist):
    '''Get a population with markers from given range
    
        HapMap_dir: where HapMap data in simuPOP format is stored. The files
            should have been prepared by scripts/loadHapMap.py.
    
        pops:     HapMap populations to load

        chrom:    chromosome number (1-based index)
        
        startPos: starting position (in cM)
        
        endPos:   ending position (in cM). If 0, ignore this parameter.
        
        maxNum:   maximum number of markers to get. If 0, ignore this parameter.
        
        minAF:    minimal minor allele frequency
        
        minDiffAf: minimal allele frequency between HapMap populations.
        
        minDist:  minimal distance between two adjacent markers, in cM
        
    '''
    print "Loading HapMap population hapmap_%d.bin" % chrom
    pop = LoadPopulation(os.path.join(HapMap_dir, 'hapmap_%d.bin' % chrom))
    markers = []
    lastPos = 0
    for loc in range(pop.totNumLoci()):
        pos = pop.locusPos(loc)
        if pos < startPos or (endPos > 0 and pos > endPos):
            continue
        if lastPos != 0 and pos - lastPos < minDist:
            continue
        maf = min(pop.dvars().alleleFreq[loc][0], 1 - pop.dvars().alleleFreq[loc][0])
        if maf < minAF:
            continue
        if minDiffAF > 0:
            maf1 = min(pop.dvars(pops[0]).alleleFreq[loc][0], 1 - pop.dvars(pops[0]).alleleFreq[loc][0])
            maf2 = min(pop.dvars(pops[1]).alleleFreq[loc][0], 1 - pop.dvars(pops[1]).alleleFreq[loc][0])
            if abs(maf1 - maf2) < minDiffAF:
                continue
        # this marker is OK
        markers.append(loc)
        if maxNum > 0 and len(markers) == maxNum:
            break
        lastPos = pos
    print '%d markers located' % len(markers)
    pop.removeLoci(keep=markers)
    pop.removeSubPops([x for x in range(3) if x not in pops], removeEmptySubPops=True)
    return pop


###########################################################
#
# Evole the hapmap population
#
# NOTE1: use uniform recombination rate for now.
# 
# NOTE2: no mutation at the last generation to avoid
#        mendelian inconsistency error.
# 
###########################################################
def evolveHapMap(pop, 
    endingSize, 
    gen,
    migr=noneOp(),
    expand='exponential',
    mergeAt=10000, 
    initMultiple=1, 
    recIntensity=0.01, 
    mutRate=1e-7,
    step=10, 
    keepParents=False, 
    numOffspring=1):
    ''' Evolve and expand the hapmap population
    
    gen: total evolution generation
    
    initMultiple: copy each individual initMultiple times, to avoid
        rapid loss of genotype variation when population size is small.
    
    endingSize: ending poplation size
    
    expand: expanding method, can be linear or exponential
    
    mergeAt: when to merge population?
    
    gen: generations to evolve
    
    recIntensity: recombination intensity
    
    mutRate: mutation rate
    
    step: step at which to display statistics
    
    keepParents: whether or not keep parental generations
    
    numOffspring: number of offspring at the last generation

    migr: a migrator to be used.
    
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
        rate = (N1-N0)*1.0/gen
    else:
        rate = math.exp(math.log(N1*1.0/N0)/gen)
    def popSizeFunc(gen, cur):
        if expand == 'linear':
            sz = [int(x+rate/len(cur)) for x in cur]
        else:
            sz = [int(x*rate) for x in cur]
        return sz
    #
    simu = simulator(pop, randomMating(newSubPopSizeFunc=popSizeFunc), rep=1)
    operators = [
        # mutation will be disallowed in the last generation (see later)
        kamMutator(rate=mutRate, loci=range(pop.totNumLoci())),
        mergeSubPops(subPops=range(pop.numSubPop()), removeEmptySubPops=True, at=[mergeAt]),
        recombinator(intensity=recIntensity),
        stat(popSize=True, step=step, begin=step-1),
        migr,
        pyEval(r'"gen=%d, size=%s\n" % (gen, subPopSize)', step=step, begin=step-1)
    ]
    del pop
    simu.evolve(ops=operators, end=gen-1)
    if keepParents:
        print "Preparing the last generation"
        simu.addInfoFields(['father_idx', 'mother_idx'])
        simu.setAncestralDepth(1)
        simu.setMatingScheme(randomMating(numOffspring=numOffspring))
        # evolve, but without mutation
        simu.step(ops=operators[1:] + [parentsTagger()])
    return simu.getPopulation(0, True)



#
#
# Example functions to sample from resulting population
#
#
def sample1DSL(pop, DSL, DA, pene, name, sampleSize):
    '''Sample from the final population, using a single lcous penetrance model.

    DSL: disease locus
    
    DA: disease allele
    
    pene: penetrance
    
    name: name of directory to save (it must exist)
    
    sampleSize: sample size, in this case, sampleSize/4 is the number of families
    '''
    # applying penetrance
    def peneFunc(geno):
        if DA == 1:
            return pene[geno[0] + geno[1]]
        else:
            return pene[2 - geno[0] - geno[1]]
    PyPenetrance(pop, locus=DSL, func=peneFunc)
    Stat(pop, numOfAffected=True, alleleFreq=[DSL])
    #
    dsl = open(os.path.join(name, 'sample.dsl'), 'w')
    print >> dsl, 'DSL (0 indexed), Name, Allele freq (allele 0), Number of affected (%)'
    print >> dsl, '%d, %s, %.4f, %d (%.3f)' % (DSL, pop.locusName(DSL), pop.dvars().alleleFreq[DSL][0], \
        pop.dvars().numOfAffected, pop.dvars().numOfAffected*1.0/pop.popSize())
    dsl.close()
    #
    # draw sample
    samples = AffectedSibpairSample(pop, size=sampleSize/4)
    SaveQTDT(samples[0], output=os.path.join(name, 'sample'), 
        affectionCode=['1', '2'], fields=['affection'])


def sample2DSL(pop, DSL, pene, name, size):
    '''Sample from the final population, using a two locus penetrance model
    
    DSL: disease loci (two locus)
    
    pene: penetrance value, assuming a two-locus model
    
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



