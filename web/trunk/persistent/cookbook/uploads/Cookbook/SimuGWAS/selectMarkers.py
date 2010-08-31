#!/usr/bin/env python
#
# Purpose:
#     This python module provides function getHapMapMarkers to select markers
#     from HapMap2 or 3 datasets.
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
#         Add functions mergeHapMapPops and getHapMapMarkers
# 

'''
This python module provides several utility functions that handles HapMap
populations. When used as a script, this module creates a population using
selected markers and populations.
'''

import simuOpt
simuOpt.setOptions(alleleType='binary', quiet=True, version='1.0.1')
from simuPOP import *

from types import *
import os, sys, exceptions

def mergeHapMapPops(HapMap_dir, HapMap_pops, chrom, logger=None):
    '''
    Load HapMap dataset for multiple populations and merge them.
    The important step is to find a common set of markers and make sure
    alleles are recoded if necessary. (Alleles are sometimes coded differently
    in different HapMap populations.)

    HapMap_dir
        Directory where HapMap files (in simuPOP format) are saved.

    HapMap_pops
        HapMap populations to load.

    chrom
        Which chromosome to load and combine.

    logger
        A logger to record what is going on.
    '''
    pop = None
    for HapMap_pop in HapMap_pops:
        filename = os.path.join(HapMap_dir, '%s_chr%d.pop' % \
            (HapMap_pop, chrom))
        if logger:
            logger.info('Loading HapMap population %s' % filename)
        pop1 = loadPopulation(filename)
        if pop is None:
            pop = pop1
            continue
        # need to be merged.
        markers1 = set(pop.lociNames())
        markers2 = set(pop1.lociNames())
        common_markers = markers1 & markers2
        remove1 = markers1 - common_markers
        remove2 = markers2 - common_markers
        if len(remove1) > 0:
            if logger:
                logger.info('Removing %d markers (%.2f percent) from population %s' % \
                    (len(remove1), len(remove1)*100./pop.totNumLoci(), pop.subPopNames()))
            pop.removeLoci([pop.locusByName(x) for x in remove1])
        if len(remove2) > 0:
            if logger:
                logger.info('Removing %d markers (%.2f percent) from population %s' % \
                    (len(remove2), len(remove2)*100./pop1.totNumLoci(), pop1.subPopNames()))
            pop1.removeLoci([pop1.locusByName(x) for x in remove2])
        # we need to recode alleles if there is one allele in one population and 
        # two alleles in another one at some loci.
        pop_loci = {0:[], 1:[]}
        pop_names = {0:[], 1:[]}
        pop1_loci = {0:[], 1:[]}
        pop1_names = {0:[], 1:[]}
        for i in range(pop.totNumLoci()):
            if pop.alleleNames(i) != pop1.alleleNames(i):
                alleles = list(set(pop.alleleNames(i) + pop1.alleleNames(i)))
                #if '-' in alleles:
                    # hapmap 2 has something like this for unknown allele
                #    alleles.remove('-')
                alleles.sort()
                if len(alleles) != 2:
                    raise exceptions.ValueError("Do not know how to recode alleles. Alleles "
                        "for locus %s in one population are %s, and %s in another. This is usually "
                        "caused by the use of data from different HapMap versions." % \
                        (pop.locusName(i), str(pop.alleleNames(i)), str(pop1.alleleNames(i))))
                if pop.alleleNames(i) != tuple(alleles):
                    toAllele = alleles.index(pop.alleleNames(i)[0])
                    pop_loci[toAllele].append(i)
                    pop_names[toAllele].append(alleles)
                if pop1.alleleNames(i) != tuple(alleles):
                    toAllele = alleles.index(pop1.alleleNames(i)[0])
                    pop1_loci[toAllele].append(i)
                    pop1_names[toAllele].append(alleles)
            if pop.locusName(i) != pop1.locusName(i):
                print 'Locus names are different', pop.locusName(i), pop1.locusName(i)
                sys.exit(0)
        if len(pop_loci[0]) > 0 or len(pop_loci[1]) > 0:
            # only change names
            pop.recodeAlleles(alleles=[0, 1], loci=pop_loci[0], alleleNames=pop_names[0])
            # can move allele from 0 to 1
            pop.recodeAlleles(alleles=[1, 0], loci=pop_loci[1], alleleNames=pop_names[1])
        if len(pop1_loci[0]) > 0 or len(pop1_loci[1]) > 0:
            pop1.recodeAlleles(alleles=[0, 1], loci=pop1_loci[0], alleleNames=pop1_names[0])
            pop1.recodeAlleles(alleles=[1, 0], loci=pop1_loci[1], alleleNames=pop1_names[1])
        # merge two populations
        # check coding.
        pop.addIndFrom(pop1)
    return pop



def getHapMapMarkers(HapMap_dir, names = [], chroms=[], HapMap_pops=['HapMap2_CEU'],
        startPos = [], endPos = [], numMarkers = [], minAF = 0, minDist = 0,
        mergeSubPops = True, logger=None):
    '''
    Return a population with specified HapMap markers.

    HapMap_dir
        Directory where the HapMap data has been saved, using script
        loadHapMap2.py or loadHapMap3.py from the simuPOP cookbook.

    names
        (Optional) A list of marker names. If given, only markers in this
        list will be selected.

    chroms
        A list of chromosomes to look in. If empty, all 22 autosomes
        will be tried. Chromosome index starts from 1. (1, ..., 22).

    HapMap_pops
        HapMap populations to load, can be one or both of 'CEU' and 'YRI'.

    startPos, endPos, numMarkers
        These list should be empty or match the length of ``chroms``.
        They specify the starting, ending position on each chromosome
        (in basepair), and number of markers to load.

    minAF
        Minimal minor allele frequency

    minDist
        Minimal distance between adjacent markers

    merSubPops
        Whether or not multiple populations should be merged to one.

    logger
        A logger to record what is going on.
    '''
    def paramExpandList(param, size, err=''):
        '''If parameter param is
        - a number: return a list of specified size.
        - a list of size 1: expand it to size and return.
        - a list of more than one element: raise an error if size mismatch.
        - an empty list: return param.
        '''
        if type(param) in [IntType, LongType, FloatType]:
            return [param]*size
        elif type(param) in [TupleType, ListType]:
            if len(param) == 1:
                return list(param)*size
            elif len(param) == 0:
                return param
            elif len(param) != size:
                raise exceptions.ValueError(err)
        return param
    if len(chroms) == 0:
        chs = range(1, 23)
    else:
        chs = chroms
    # read in HapMap data file
    pop = None
    geneticMap = {}
    sPos = paramExpandList(startPos, len(chs), 'Incorrect starting position')
    ePos = paramExpandList(endPos, len(chs), 'Incorrect ending position')
    nMarkers = paramExpandList(numMarkers, len(chs), 'Incorrect number of markers')
    #
    for chIdx, ch in enumerate(chs):
        markers = []
        chPop = mergeHapMapPops(HapMap_dir, HapMap_pops, ch, logger)
        # Trim markers by marker names
        if len(names) != 0:
            if logger:
                logger.info("Select markers using a list of %d markers from chromosome %s..." % (len(names), ch))
            chPop.removeLoci(keep = chPop.lociByNames(list(set(chPop.lociNames()) & set(names))))
        #
        if chPop.totNumLoci() == 0:
            continue
        if minAF > 0:
            stat(chPop, alleleFreq=range(chPop.totNumLoci()))
        # Trim by start, end position ...
        indexes = []
        lastPos = 0
        for loc in range(chPop.totNumLoci()):
            pos = chPop.locusPos(loc)
            if len(sPos) > 0 and pos < sPos[chIdx]:
                continue
            if len(ePos) > 0 and pos > ePos[chIdx]:
                continue
            if lastPos > 0 and pos - lastPos < minDist:
                continue
            if minAF > 0:
                maf = chPop.dvars().alleleFreq[loc][0]
                maf = min(maf, 1 - maf)
                if maf < minAF:
                    continue
            if len(nMarkers) > 0 and len(indexes) >= nMarkers[chIdx]:
                break
            indexes.append(loc)
            lastPos = pos
        if len(indexes) > 0:
            if logger:
                logger.info('%s markers are found on chromosome %d ' % (len(indexes), ch))
            chPop.removeLoci(keep=indexes)
            geneticMap.update(chPop.dvars().geneticMap)
            chPop.vars().clear()
            if pop is None:
                pop = chPop
            else:
                pop.addChromFrom(chPop)
        else:
            if logger:
                logger.info('No qualified marker is found on chromosome %d ' % ch)
            del chPop
    if pop.numSubPop() > 1 and mergeSubPops:
        pop.mergeSubPops(range(pop.numSubPop()))
    pop.dvars().geneticMap = geneticMap
    return pop

def saveMarkerList(pop, filename, logger=None):
    '''Save a marker list file'''
    stat(pop, alleleFreq=range(pop.totNumLoci()))
    lst = open(filename, 'w')
    print >> lst, 'name\tchrom\tpos\tallele1\tallele2\tfreq_of_allele2' 
    for ch in range(pop.numChrom()):
        chName = pop.chromName(ch)
        for loc in range(pop.chromBegin(ch), pop.chromEnd(ch)):
            alleleNames = pop.alleleNames(loc)
            if len(alleleNames) == 1:
                alleleNames = (alleleNames[0], '-')
            print >> lst, '%s\t%s\t%d\t%s\t%.4f' % (pop.locusName(loc), chName, int(pop.locusPos(loc)),
                '\t'.join(alleleNames), pop.dvars().alleleFreq[loc][1])
    lst.close()


HapMap2_pops = ['HapMap2_CEU', 'HapMap2_YRI', 'HapMap2_JPT+CHB']
HapMap3_pops = ['HapMap3_ASW', 'HapMap3_CEU', 'HapMap3_CHD', 'HapMap3_GIH', 
    'HapMap3_JPT+CHB', 'HapMap3_LWK', 'HapMap3_MEX', 'HapMap3_MKK', 'HapMap3_TSI', 'HapMap3_YRI']

options = [
    {
    'longarg': 'HapMap_dir=',
    'default': '.',
    'useDefault': True,
    'label': 'HapMap data directory',
    'description': '''Directory to store HapMap data in simuPOP format. The 
        HapMap data is expected to be downloaded and saved in simuPOP format
        using script loadHapMap_r22.py from the simuPOP online cookbook. The
        files have names such as HapMap_CEU_chr10.pop.''',
     'allowedTypes': [StringType],
     'validate': simuOpt.valueValidDir(),
    },
    {
    'longarg': 'HapMap_pops=',
    'default' : HapMap3_pops,
    'useDefault': False,
    'label' : 'Name of populations',
    'description': '''Which HapMap populations to use?''',
    'allowedTypes': [ListType, TupleType],
    'chooseFrom': HapMap2_pops + HapMap3_pops,
    'validate': simuOpt.valueListOf(HapMap2_pops + HapMap3_pops),
    },
    {
    'longarg': 'markerList=',
    'default': '',
    'useDefault': True,
    'label': 'Marker list file',
    'description': '''A file with a list of marker names. If there are more than
        one fields at a line, the rest of them are ignored.''',
    'allowedTypes': [StringType],
    'validate': simuOpt.valueOr(simuOpt.valueEqual(''), simuOpt.valueValidFile()),
    },
    {
    'longarg': 'chroms=',
    'default': [],
    'useDefault': True,
    'label': 'Chromosomes to use',
    'description': 'A list of chromosomes (1-22) to use.',
    'allowedTypes': [TupleType, ListType],
    'validate': simuOpt.valueListOf(simuOpt.valueBetween(1, 22)),
    },
    {
    'longarg': 'numMarkers=',
    'default': [],
    'useDefault': True,
    'label': 'Number of markers to use',
    'description': '''Number of markers to use for each chromosome. This
        parameter should be ignored if it is unspecified or is set to zero
        for some chromosomes.
        ''',
    'allowedTypes': [TupleType, ListType],
    'validate': simuOpt.valueOr(simuOpt.valueGT(0), simuOpt.valueListOf(simuOpt.valueGE(0)))
    },
    {
    'longarg': 'startPos=',
    'default': [],
    'useDefault': True,
    'label': 'staring position (bp)',
    'description': '''Starting position of the markers on each chromosome.
        The beginning of the chromosomes will be assumed if this parameter
        is unspecified or is set to zero.''',
    'allowedTypes': [TupleType, ListType],
    'validate': simuOpt.valueOr(simuOpt.valueGE(0), simuOpt.valueListOf(simuOpt.valueGE(0)))
    },
    {
    'longarg': 'endPos=',
    'default': [],
    'useDefault': True,
    'label': 'Ending position (bp)',
    'description': '''Ending position of the markers on each chromosome.
        The end of the chromosomes will be assumed if this parameter is
        unspecifed or is set to zero.''',
     'allowedTypes': [TupleType, ListType],
     'validate': simuOpt.valueOr(simuOpt.valueGE(0), simuOpt.valueListOf(simuOpt.valueGE(0)))
    },
    {
    'longarg': 'minAF=',
    'default': 0,
    'useDefault': True,
    'label': 'Minimal minor allele frequency',
    'description': '''Minimal allele frequency of selected markers.''',
    'allowedTypes': [IntType, LongType, FloatType],
    'validate': simuOpt.valueBetween(0, 0.5)
    },
    {
    'longarg': 'minDist=',
    'default': 0,
    'useDefault': True,
    'label': 'Minimal distance between markers',
    'description': '''Minimal distance between adjacent markers''',
    'allowedTypes': [IntType, LongType, FloatType],
    'validate': simuOpt.valueGE(0),
    },
    {
    'longarg': 'mergeSubPops',
    'default': True,
    #'useDefault': True,
    'label': 'Merge all subpopulations',
    'description': '''Merge all subpopulations''',
    'allowedTypes': [BooleanType],
    },
    {
    'longarg': 'filename=',
    'label': 'Filename to save population',
    'default': 'result.pop',
    'useDefault': False,
    'description': '''Name of the population or an absolute path to
        a file. This parameter will be ignored if an empty string or None
        is given.''',
    'allowedTypes': [StringType],
    'validate': simuOpt.valueNotEqual('')
    }
]

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger()
    pars = simuOpt.Params(options,
        'This script chooses specified markers from one or more HapMap\n'
        'populations and saves them in simuPOP format.\n',
        __doc__)
    if not pars.getParam():
        sys.exit(1)
    names = []
    if pars.markerList != '':
        mlist = open(pars.markerList)
        for line in mlist.readlines():
            if line.startswith('#') or line.strip() == '':
                continue
            names.append(line.split(',')[0].split()[0])
        if logger:
            logger.info('%d markers located from marker list file %s' %\
                (len(names), pars.markerList))
    #
    pop = getHapMapMarkers(pars.HapMap_dir, 
        names = names,
        chroms=pars.chroms, 
        HapMap_pops=pars.HapMap_pops,
        startPos = pars.startPos,
        endPos = pars.endPos,
        numMarkers = pars.numMarkers,
        minAF = pars.minAF,
        minDist = pars.minDist,
        mergeSubPops = pars.mergeSubPops,
        logger=logger)
    if logger:
        logger.info('Save population to %s and marker list to %s.lst' % \
            (pars.filename, pars.filename))
    pop.save(pars.filename)
    saveMarkerList(pop, pars.filename + '.lst', logger)
                


