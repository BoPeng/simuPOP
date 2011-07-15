#!/usr/bin/env python
#
# Purpose:
#     This python module provides function loadHapMapPop to download and import
#     the HapMap3 populations. It can also be served as a script to download
#     part or all populations of the phase 3 of the HapMap dataset.
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
#     2009-08-31 Bo Peng <bpeng@mdanderson.org>
#       - Add logging support to script loadHapMap_r22.py.
#       - Polish function loadHapMapPop and hide others.
#       - Use temporary directory to save downloaded files.
# 
'''
This script downloads and loads release 2 of hapmap phase 3 datasets in 
ftp://ftp.ncbi.nlm.nih.gov/hapmap//phasing/2009-02_phaseIII/HapMap3_r2/
also downloads the fine-scale recombination map from
http://ftp.hapmap.org/recombination/2008-03_rel22_B36/rates/
and saves the genetic distance of each marker in a dictionary
(geneticMap) in each population's local namespace.

The saved populations have the following features:

1. Different populations are saved in different files. These populations
  may not be merged directly because they have different set of markers.
  Subpopulation name is specified ('ASW', 'CEU', 'CHB'...).

2. Chromosome name is saved as "1", "2", "3", ...

3. Basepairs are used to specify physical distances of loci.

4. Alleles are saved as 0 and 1 as appear in the HapMap datafile. Allele
  names such as 'A', 'G' are saved for each marker.

5. A dictionary 'geneticMap' is used to store genetic distance of each marker.

'''

from simuOpt import *
setOptions(optimized=True, alleleType='binary', version='1.0.1')
from simuPOP import *

import os, sys, urllib, gzip, tempfile, shutil, time

HapMap3_pops = ['ASW', 'CEU', 'CHD', 'GIH', 'JPT+CHB', 'LWK', 'MEX', 'MKK', 'TSI', 'YRI']
HapMap3_pop_types = {
    'ASW': ('TRIOS', 'DUOS', 'UNRELATED'),
    'CEU': ('TRIOS', 'DUOS', 'UNRELATED'),
    'CHD': ('',),
    'GIH': ('',),
    'JPT+CHB': ('',),
    'LWK': ('',),
    'MEX': ('TRIOS', 'DUOS'),
    'MKK': ('TRIOS', 'UNRELATED'),
    'TSI': ('',),
    'YRI': ('TRIOS', 'DUOS', 'UNRELATED'),
}

HapMap3_pop_sizes = {
    'ASW': 53,
    'CEU': 113,
    'CHD': 85,
    'GIH': 88,
    'JPT+CHB': 170,
    'LWK': 90,
    'MEX': 50,
    'MKK': 143,
    'TSI': 88,
    'YRI': 113
}


release = 2
Genotype_URL = 'ftp://ftp.ncbi.nlm.nih.gov/hapmap//phasing/2009-02_phaseIII/HapMap3_r2/%s/'
Recom_URL = 'ftp://ftp.hapmap.org/hapmap/recombination/2008-03_rel22_B36/rates/'
genotype_file = 'hapmap3_r2_b36_fwd.consensus.qc.poly.chr%d_%s.%sphased.gz'
recom_file = 'genetic_map_chr%s_b36.txt'

def downloadIfNeeded(URL, path, file, logger=None):
    '''Download file from hapmap website'''
    diskfile = os.path.join(path, file)
    # this actually will not happen because files are downloaded to
    # a temporary directory.
    while True:
        if os.path.isfile(diskfile):
            return diskfile      
        try:
            urllib.urlretrieve('%s/%s' % (URL, file), diskfile)
            if logger is not None:
                logger.info('%s is downloaded.' % file)
        except:
            raise SystemError('Failed to download file %s from URL %s' \
                % (file, URL))
            time.sleep(5)
    return diskfile


def _probeInfo(datafile, logger=None):
    '''Get population size of a sample
      This function also checks if each line has desired number of alleles.
    '''
    count = 0
    ll = 0
    data = gzip.open(datafile)
    # rsID position name_A name_B name1_A name1_B     for others
    # rsID position name_A name_B name_0_A ...        for TIOS
    firstline = data.readline()
    # does not count _0_A
    numInds = len([x for x in firstline.strip().split()[2:] if x.endswith('_B')])
    name = []
    alleleNames = {}
    pos = []
    for line in data.readlines():
        data = line.split()
        name.append(data[0])
        pos.append(int(data[1]))
        alleles = list(set(data[2:]))
        alleles.sort()
        alleleNames[data[0]] = alleles
    return numInds, name, alleleNames, pos
    

def load_population(pop, diskFiles, alleleNames, logger=None):
    '''Load population from file, with type (subpopulation type)'''
    # file format:
    #
    # rsID pos ind1_A ind2_A ....
    #
    ind_base = 0
    for diskfile in diskFiles:
        data = gzip.open(diskfile)
        fields = data.readline().split()[2:]
        numInd = len([x for x in fields if x.endswith('_B')])
        indCols = len(fields) / numInd
        if indCols * numInd != len(fields):
            raise SystemError('Something wrong with individual count.')
        if logger:
            logger.info("Importing genotypes of %d individuals from %s..." % (numInd, os.path.split(diskfile)[-1]))
        for line_no,line in enumerate(data.readlines()):
            fields = line.split()
            name = fields[0]
            alleleName = alleleNames[name]
            genotype = [alleleName.index(x) for x in fields[2:]]
            for col_no,geno in enumerate(genotype):
                ind = ind_base + col_no / indCols
                ploidy = col_no % indCols
                if ploidy == 2:
                    continue
                if ind == pop.popSize():
                    print 'Warning: individual index %d greater than population size %d ' % (ind, pop.popSize())
                # always chromosome 0, because each population has only one chromosome
                assert pop.locusPos(line_no) == float(fields[1])
                assert pop.locusName(line_no) == name
                pop.individual(ind).setAllele(geno, line_no, ploidy)
        ind_base += numInd


def set_map_dist(pop, ch, dest, logger=None):
    '''Set map distance for each locus'''
    file = recom_file % ch
    downloadIfNeeded(Recom_URL, dest, file, logger)
    if logger is not None:
        logger.info('Using genetic map file %s' % file)
    dist = {}
    for line in open(os.path.join(dest, file)).readlines():
        try:
            fields = line.split()
            pos = int(fields[0])
            dist[pos] = float(fields[2])
        except:
            pass
    if logger is not None:
        logger.info('Map distance of %d markers are found' % len(dist))
    totNumLoci = pop.totNumLoci()
    # now, try to set genetic map
    map_dist = [-1]*totNumLoci;
    cnt = 0
    for loc in range(totNumLoci):
        pos = int(int(pop.locusPos(loc)))
        try:
            map_dist[loc] = dist[pos]
            cnt += 1
        except:
            pass
    if cnt != totNumLoci:
        prev = -1     # name of the previous marker with map distance
        next = -1     # name of the next marker with map distance
        loc = 0
        while (loc < totNumLoci):
            # already has value
            if map_dist[loc] != -1:
                prev = next = loc
                loc += 1
                continue
            # find the ext one
            if prev == next:
                next = -1
                for n in range(loc+1, totNumLoci):
                     if map_dist[n] != -1:
                         next = n
                         next_pos = pop.locusPos(next)
                         next_dis = map_dist[next]
                         break
                if prev == -1:
                    for n in range(next):
                        # rough estimation: distance (in cM) proportional to 0.01 recombination rate
                        map_dist[n] = map_dist[next] - (pop.locusPos(next) - pop.locusPos(n)) * 1e-8
                # if not found, this is at the end of a chromosome
                elif next == -1:
                    for n in range(loc, totNumLoci):
                        map_dist[n] = map_dist[prev] + (pop.locusPos(n) - pop.locusPos(prev)) * 1e-8
                    break
                # if found, but no previous, this is the first one
                else:
                    prev_pos = pop.locusPos(prev)
                    prev_dis = map_dist[prev]
                    for n in range(loc, next):
                        map_dist[n] = map_dist[prev] + (pop.locusPos(n) - prev_pos) / (next_pos - prev_pos) * (next_dis - prev_dis)
                prev = next
                loc = next + 1
    if totNumLoci != cnt and logger is not None:
        logger.info('Map distance of %d markers (%.2f%% of %d) are estimated' % (totNumLoci - cnt,
            (totNumLoci - cnt) * 100.0/totNumLoci, totNumLoci))
    map = {}
    for loc in range(totNumLoci):
        map[pop.locusName(loc)] = map_dist[loc]
    pop.dvars().geneticMap = map

    
def loadHapMapPop(chrom, popName, logger=None):
    '''Download and import the specified chromosome of a hapmap population.
    If a directory is specified, the loaded population will be saved in
    simuPOP format with filename HapMap_XXX_chrY.pop where XXX is population
    name and Y is chromosome number. If a file already exists, this function
    will try to load the file directly.
    
    chrom
        chromosome to download (1, 2, ..., 22.

    popName
        Name of the population, should be one of HapMap3_pops

    logger
        An optional logger object (c.f. the Python logging module) where all
        logging and debugging information is written to.

    This function returns the loaded population.
    '''

    if logger is not None:
        logger.info("Loading HapMap3 chromosome %d of population %s" % (chrom, popName))
    tmpdir = tempfile.mkdtemp()
    URL = Genotype_URL % (popName.upper())
    sampleCode = {'': 'unr.', 'DUOS':'D.', 'TRIOS':'', 'UNRELATED':'unr.'}
    totNumInds = 0
    allAlleleNames = {}
    allLociPos = []
    allLociNames = []
    diskFiles = []
    for sampleType in HapMap3_pop_types[popName]:
        datafile = genotype_file % (chrom, popName.lower(), sampleCode[sampleType])
        if sampleType == '':
            diskfile = downloadIfNeeded(URL, tmpdir, datafile, logger)
        else:
            diskfile = downloadIfNeeded(URL + '/' + sampleType, tmpdir, datafile, logger)
        (numInds, lociNames, alleleNames, lociPos) = _probeInfo(diskfile)
        diskFiles.append(diskfile)
        if logger:
            logger.info('Genotypes of %d individuals at %d loci are found' % (numInds, len(lociNames)))
        totNumInds += numInds
        if not allLociPos:
            allLociPos = lociPos
        else:
            if allLociPos != lociPos:
                raise ValueError("Loci position mismatch.")
        if not allLociNames:
            allLociNames = lociNames
        else:
            if allLociNames != lociNames:
                raise ValueError("Loci name mismatch.")
        for key in alleleNames:
            if allAlleleNames.has_key(key):
                if len(allAlleleNames[key]) == 1:
                    if len(alleleNames[key]) == 2:
                        allAlleleNames[key] = alleleNames[key]
                    elif alleleNames[key][0] != allAlleleNames[key][0]:
                        alleles = [alleleNames[key][0], allAlleleNames[key][1]]
                        alleles.sort()
                        allAlleleNames[key] = alleles
            else:
                allAlleleNames[key] = alleleNames[key]
        if sampleType == '':
            break
    #
    if logger is not None:
        logger.info("Genotypes %d individuals of %d loci (%d - %d bp) are located" % (totNumInds, 
            len(allLociNames), allLociPos[0], allLociPos[-1]))
    pop = Population(size=totNumInds, ploidy=2, loci=[len(lociPos)],
        lociPos=allLociPos, lociNames=allLociNames, chromNames=[str(chrom)],
        alleleNames=[allAlleleNames[x] for x in allLociNames], subPopNames=[popName])
    load_population(pop, diskFiles, allAlleleNames, logger)
    set_map_dist(pop, chrom, tmpdir, logger)
    pop.dvars().HapMap_rel = release
    shutil.rmtree(tmpdir)
    return pop

options = [
    {'longarg': 'dest=',
     'default': '.',
     'useDefault': True,
     'label': 'Destination directory',
     'allowedTypes': [type('')],
     'validate': valueValidDir(),
     'description': 'A directory to save HapMap data in simuPOP format.',
    },
    {'longarg': 'chroms=',
     'default': range(1, 23),
     'useDefault': True,
     'label': 'Chromosomes to download',
     'allowedTypes': [type([]), type(())],
     'chooseFrom': range(1, 23),
     'validate': valueListOf(valueBetween(1, 22)),
     'description': 'Which chromosomes to download and process',
    },
]

if __name__ == '__main__':
    pars = Params(options, 
        'This script downloads the second release of the phase 3 of the HapMap datasets\n'
        'and saves them in simuPOP format. It also downloads the fine-scale\n'
        'recombination map and saves the genetic distance of each marker in\n'
        'a dictionary (geneticMap) in the population\'s local namespace.',
        __doc__)
    if not pars.getParam():
        sys.exit(1)
    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('loadHapMap3')
    for chrom in pars.chroms:
        for sample in HapMap3_pops:
            popFile = os.path.join(pars.dest, "HapMap3_%s_chr%d.pop" % (sample, chrom))
            try:
                if os.path.isfile(popFile):
                    # test if this file is OK.
                    pop = loadPopulation(popFile)
                    if pop.popSize() == HapMap3_pop_sizes[sample]:
                        logger.info("Skipping existing population %s." % popFile)
                        continue
            except:
                # continue to load file
                pass
            pop = loadHapMapPop(chrom, sample, logger)
            logger.info("Save population to %s." % popFile)
            pop.save(popFile)
