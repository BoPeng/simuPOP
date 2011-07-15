#!/usr/bin/env python
#
# Purpose:
#     This python module provides function loadHapMapPop to download and import
#     the HapMap populations. It can also be served as a script to download
#     part or all populations from the phase 2 of the HapMap dataset.
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
This script downloads and loads release 22 of hapmap datasets in 
http://www.hapmap.org/downloads/phasing/2007-08_rel22/phased. It
also downloads the fine-scale recombination map from
http://ftp.hapmap.org/recombination/2008-03_rel22_B36/rates/
and saves the genetic distance of each marker in a dictionary
(geneticMap) in each population's local namespace.

The saved populations have the following features:

1. Different populations are saved in different files. These populations
  may not be merged directly because they have different set of markers.
  Subpopulation name is specified ('CEU', 'YRI' or 'JPT+CHT'.

2. Chromosome name is saved as "1", "2", "3", ...

3. Basepairs are used to specify physical distances of loci.

4. Alleles are saved as 0 and 1 as appear in the HapMap datafile. Allele
  names such as 'A', 'G' are saved for each marker.

5. A dictionary 'geneticMap' is used to store genetic distance of each marker.

Please refer to 
http://simupop.sourceforge.net/cookbook/pmwiki.php/Cookbook/LoadHapMap22 for 
for details about this script.
'''

from simuOpt import *
setOptions(optimized=True, alleleType='binary', version='1.0.1')
from simuPOP import *

import os, sys, urllib, gzip, tempfile, shutil

# URL and revision number. You can choose to download other hapmap files
release = 22
# This link does not have JPT+CHT FILES, they are older (2007-11-20)
#Genotype_URL = 'http://www.hapmap.org/downloads/phasing/2007-08_rel22/phased/'
# Data from this link may not be correct but they are newer (2008-6-25)
Genotype_URL = 'ftp://ftp.hapmap.org/hapmap/phasing/2007-08_rel22/phased/'
Recom_URL = 'ftp://ftp.hapmap.org/hapmap/recombination/2008-03_rel22_B36/rates/'
legend_file = 'genotypes_chr%d_%s_r22_nr.b36_fwd_legend.txt.gz'
genotype_file = {
    'CEU': 'genotypes_chr%d_CEU_r22_nr.b36_fwd.phase.gz',
    'YRI': 'genotypes_chr%d_YRI_r22_nr.b36_fwd.phase.gz',
    'JPT+CHB': 'genotypes_chr%d_JPT+CHB_r22_nr.b36_fwd.phased.gz'  # Come on, what is going on? (Note the d in phased)
}
recom_file = 'genetic_map_chr%s_b36.txt'

def downloadIfNeeded(URL, path, file, logger=None):
    '''Download file from hapmap website'''
    diskfile = os.path.join(path, file)
    # this actually will not happen because files are downloaded to
    # a temporary directory.
    if os.path.isfile(diskfile):
        return        
    if logger is not None:
        logger.info('Downloading %s ...' % file)
    try:
        urllib.urlretrieve('%s/%s' % (URL, file), diskfile)
    except:
        raise SystemError('Failed to download file %s from URL %s' \
            % (file, URL))


def _getLegend(ch, sample, dest, logger=None):
    '''Loci information is retrieved from legend files'''
    lociPos = []
    lociName = []
    alleleNames = []
    file = legend_file % (ch, sample)
    downloadIfNeeded(Genotype_URL, dest, file, logger)
    try:
        legend = gzip.open(os.path.join(dest, file))
        legend.readline()  # skip first line
        for line in legend.readlines():
            fields = line.split()
            lociName.append(fields[0])
            lociPos.append(float(fields[1]))
            alleleNames.append((fields[2], fields[3]))
    except Exception, e:
        if logger is not None:
            logger.error('Failed to read file %s' % file)
            logger.error('You may want to remove your local copy and let this script re-download it.')
        raise e
    return (lociPos, lociName, alleleNames)


def _getPopSize(numLoci, ch, sample, dest, logger=None):
    '''Get population size of a sample
      This function also checks if each line has desired number of alleles.
    '''
    count = 0
    ll = 0
    genotype = genotype_file[sample] % ch
    downloadIfNeeded(Genotype_URL, dest, genotype, logger)
    for line in gzip.open(os.path.join(dest, genotype)).readlines():
        if (ll == 0 and len(line.split()) != numLoci) or \
            (ll != 0 and len(line) != ll):
            if logger is not None:
                logger.error("Number of loci does not match in %s " % ceu)
                logger.error("Number of loci: %d, number of fields: %d" % (numLoci, len(line.split())))
            sys.exit(1)
        ll = len(line)
        count += 1
    return count/2
    

def load_population(pop, ch, sample, dest, logger=None):
    '''Load population from file, with type (subpopulation type)'''
# For the CEU and YRI the haplotypes are arranged as follows:
#  
# row 1 - trio 1 parent 1 transmitted haplotype
# row 2 - trio 1 parent 1 untransmitted haplotype
# row 3 - trio 1 parent 2 transmitted haplotype
# row 4 - trio 1 parent 2 untransmitted haplotype
# row 5 - trio 2 parent 1 transmitted haplotype
# row 6 - trio 2 parent 1 untransmitted haplotype
# row 7 - trio 2 parent 2 transmitted haplotype
# row 8 - trio 2 parent 2 untransmitted haplotype
# .
# .
# For the JPT+CHB the haplotypes are arranged as
#  
# row 1 - individual 1 haplotype 1
# row 2 - individual 1 haplotype 2
# row 3 - individual 2 haplotype 1
# row 4 - individual 2 haplotype 2
# row 5 - individual 3 haplotype 1
# row 6 - individual 3 haplotype 2
# .
# .
# We are loading row by row, so actually only load the parents of CEU and YRI
# populations.
# 
    file = genotype_file[sample] % ch
    downloadIfNeeded(Genotype_URL, dest, file, logger)
    for line_no,line in enumerate(gzip.open(os.path.join(dest, file)).readlines()):
        genotype = [int(x) for x in line.split()]
        ind = int(line_no / 2)
        ploidy = line_no % 2
        # always chromosome 0, because each population has only one chromosome
        pop.individual(ind).setGenotype(genotype, ploidy)


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
        Name of the population, should be one of 'CEU', 'YRI' or 'JPT+CHB'

    logger
        An optional logger object (c.f. the Python logging module) where all
        logging and debugging information is written to.

    This function returns the loaded population.
    '''

    if logger is not None:
        logger.info("Loading HapMap chromosome %d of population %s" % (chrom, popName))
    tmpdir = tempfile.mkdtemp()
    (lociPos, lociName, alleleNames) = _getLegend(chrom, popName, tmpdir, logger)
    if logger is not None:
        logger.info("%d loci (%d - %d bp) are located" % (len(lociPos), lociPos[0], lociPos[-1]))
    popSize = _getPopSize(len(lociPos), chrom, popName, tmpdir, logger)
    if logger is not None:
        logger.info('Sample size is %d' % popSize)
    pop = Population(size=popSize, ploidy=2, loci=[len(lociPos)],
        lociPos=lociPos, lociNames=lociName, chromNames=[str(chrom)],
        alleleNames=alleleNames, subPopNames=[popName])
    load_population(pop, chrom, popName, tmpdir, logger)
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
        'This script downloads the 22 release of the HapMap datasets\n'
        'and saves them in simuPOP format. It also downloads the fine-scale\n'
        'recombination map and saves the genetic distance of each marker in\n'
        'a dictionary (geneticMap) in the population\'s local namespace.',
        __doc__)
    if not pars.getParam():
        sys.exit(1)
    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('loadHapMap2')
    for chrom in pars.chroms:
        for sample in ['CEU', 'YRI', 'JPT+CHB']:
            popFile = os.path.join(pars.dest, "HapMap2_%s_chr%d.pop" % (sample, chrom))
            try:
                if os.path.isfile(popFile):
                    # test if this file is OK.
                    loadPopulation(popFile)
                    logger.info("Population %s already exists. Please remove it first if you would like to regenerate this file." % popFile)
                    continue
            except:
                # continue to load file
                pass
            pop = loadHapMapPop(chrom, sample, logger)
            logger.info("Save population to %s." % popFile)
            pop.save(popFile)
