#!/usr/bin/env python

'''
This script downloads and loads hapmap data in 
http://www.hapmap.org/downloads/phasing/2006-07_phaseII/consensus/
I chose consensus data because these SNPs exist in all three populations,
and I chose phased data because this information is important for the 
followup applications.

Estimated fine-scale recombination rates between adjacent markers are
also downloaded from http://ftp.hapmap.org/recombination/2006-10_rel21_phaseI+II/rates/
The result is saved in a population dictionary as genDist, which stores
the map distance of each locus.

Run this script in the directory. You will obtain files named hapmap_XX.bin,
where XX is chromosome number. These populations can be loaded and used in
followup applications. You should of course tailor the population (choose SNP
markers etc) using applicable simuPOP population related functions. If you
would like to combine two populations, use pop.mergePopulationByLoci(pop1)

NOTE:

1. Actual allele information is lost. In the legent+freq files, alleles coded
   in 0 and 1 are specified. This information is not available in the simuPOP
   population.

2. There is no missing data in the fwd_phased_consensus files.

3. Sex chromosome is currently ignored.

4. For CEU and YRI populations, only parents (60 each) are loaded. This is
   because we want as many individuals as possible, but no duplicate
   genotype of offspring.

5. Allele frequency is clculated (and stored) in the population files.
   This makes it easy to pick locus according to allele frequency. Note
   that these frequencies are counted directly from the sample. No estimation
   is used.

6. Loci positions are translated from sequence location to cM by deviding indices
   by 1000000, roughly in the unit of cM.

'''

from simuOpt import setOptions
setOptions(optimized=True, alleleType='binary')
from simuPOP import *

import os, sys, urllib, gzip, exceptions

# URL and revision number. You can choose to download other hapmap files
Genotype_URL = 'http://www.hapmap.org/downloads/phasing/2006-07_phaseII/consensus/'
Recom_URL = 'http://ftp.hapmap.org/recombination/2006-10_rel21_phaseI+II/rates/'
rev = 21
legend_freq_file = 'genotypes_chr%d_r%d_nr_fwd_consensus_legend+freq.gz'
genotype_file = 'genotypes_chr%d_%s_r%d_nr_fwd_phased_consensus.gz'
recom_file = 'genetic_map_chr%s.txt.gz'

def downloadIfNeeded(URL, path, file):
    '''Download file from hapmap website'''
    diskfile = os.path.join(path, file)
    if os.path.isfile(diskfile):
        return        
    print 'Downloading %s to %s ...' % (file, path)
    urllib.urlretrieve('%s/%s' % (URL, file), diskfile)
    if not os.path.isfile(diskfile):
        raise exceptions.SystemError('Failed to download file %s from URL %s' \
            % (file, URL))


def getLoci(ch, dest):
    '''Loci information is retrieved from _consensus_legend+freq files'''
    lociPos = []
    lociName = []
    file = legend_freq_file % (ch, rev)
    downloadIfNeeded(Genotype_URL, dest, file)
    try:
        for line in gzip.open(os.path.join(dest, file)).readlines():
            fields = line.split()
            lociName.append(fields[0])
            # translate pos from index to cM. This is tentative
            lociPos.append(float(fields[1])/1000000.)
    except:
        print 'Failed to read file ', file
        print 'You may want to remove your local copy and let this script re-download it.'
        print
        raise
    return (lociPos, lociName)


def getPopSize(numLoci, ch, dest):
    '''Get population size of each subpopulations
      This function also checks if each line has desired number of alleles.
    '''
    count = [0]*3
    ll = 0
    ceu = genotype_file % (ch, 'CEU', rev)
    downloadIfNeeded(Genotype_URL, dest, ceu)
    for line in gzip.open(os.path.join(dest, ceu)).readlines():
        if (ll == 0 and len(line.split()) != numLoci) or \
            (ll != 0 and len(line) != ll):
            print "Number of loci does not match in %s " % ceu
            print "Number of loci: %d, number of fields: %d" % (numLoci, len(line.split()))
            sys.exit(1)
        ll = len(line)
        count[0] += 1
    ll = 0
    yri = genotype_file % (ch, 'YRI', rev)
    downloadIfNeeded(Genotype_URL, dest, yri)
    for line in gzip.open(os.path.join(dest, yri)).readlines():
        if (ll == 0 and len(line.split()) != numLoci) or \
            (ll != 0 and len(line) != ll):
            print "Number of loci does not match in %s " % yri
            print "Number of loci: %d, number of fields: %d" % (numLoci, len(line.split()))
            sys.exit(1)
        ll = len(line)
        count[1] += 1
    ll = 0
    jpt_chb = genotype_file % (ch, 'JPT+CHB', rev)
    downloadIfNeeded(Genotype_URL, dest, jpt_chb)
    for line in gzip.open(os.path.join(dest, jpt_chb)).readlines():
        if (ll == 0 and len(line.split()) != numLoci) or \
            (ll != 0 and len(line) != ll):
            print "Number of loci does not match in %s " % jpt_chb
            print "Number of loci: %d, number of fields: %d" % (numLoci, len(line.split()))
            sys.exit(1)        
        ll = len(line)
        count[2] += 1
    return [x/2 for x in count]       
    

def load_population(pop, ch, type, dest):
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
    subPop = {'CEU':0, 'YRI':1, 'JPT+CHB':2}[type]
    file = genotype_file % (ch, type, rev)
    downloadIfNeeded(Genotype_URL, dest, file)
    print 'from %s...' % file
    for line_no,line in enumerate(gzip.open(os.path.join(dest, file)).readlines()):
        genotype = [int(x) for x in line.split()]
        ind = line_no / 2
        ploidy = line_no % 2
        ind = pop.individual(ind, subPop)
        for i,g in enumerate(genotype):
            # always chromosome 0, because each population has only one chromosome
            ind.setAllele(g, i, ploidy)


def set_map_dist(pop, ch, dest):
    '''Set map distance for each locus'''
    file = recom_file % ch
    downloadIfNeeded(Recom_URL, dest, file)
    print 'Reading genetic map file %s...' % file,
    dist = {}
    for line in gzip.open(os.path.join(dest, file)).readlines():
        try:
            fields = line.split()
            pos = int(fields[0])
            dist[pos] = float(fields[2])
        except:
            pass
    # now, try to set genetic map
    map = {}
    for loc in range(pop.totNumLoci()):
        pos = int(1000000*pop.locusPos(loc))
        name = pop.locusName(loc)
        try:
            map[name] = dist[pos]
        except:
            pass
    print ' map distance of %d markers are found' % len(map)
    prev = -1     # name of the previous marker with map distance
    next = -1     # name of the next marker with map distance
    cnt = 0
    for loc in range(pop.totNumLoci()):
        pos = pop.locusPos(loc)
        name = pop.locusName(loc)
        if map.has_key(name):
            prev = next = loc
        else:
            cnt += 1
            # if it is before the first one, map distance is 0. 
            if prev == -1:
                map[name] = 0.
                continue
            # if it is after the last one, map distance is the last one
            elif next == -1:
                map[name] = map[pop.locusName(prev)]
                continue
            # neither prev of next is ''
            if prev == next: # find the next!!
                next = -1
                for n in range(loc+1, pop.totNumLoci()):
                    if map.has_key(pop.locusName(n)):
                        next = n
                        prev_pos = pop.locusPos(prev)
                        next_pos = pop.locusPos(next)
                        prev_dis = map[pop.locusName(prev)]
                        next_dis = map[pop.locusName(next)]
                        break
            # if not found
            if next == -1:
                map[name] = map[pop.locusName(prev)]
                continue
            else:
                map[name] = map[pop.locusName(prev)] + (pos - prev_pos) / (next_pos - prev_pos) * (next_dis - prev_dis)
    print 'Map distance of %d markers (%.2f%% of %d) are estimated' % (cnt, cnt*100.0/pop.totNumLoci(), pop.totNumLoci())
    pop.dvars().genDist = map

    
def loadHapMap(chroms, dest='.'):
    '''Download, import and save hapmap data of given chromosomes'''
    ps = [0,0,0]
    for ch in chroms:
        popFile = "hapmap_%d.bin" % ch
        print "\n\nLoading HapMap chromosome %d" % ch
        (lociPos, lociName) = getLoci(ch, dest)
        print "%d loci (%.2f - %.2f cM) are located" % (len(lociPos), lociPos[0], lociPos[-1])
        popSize = getPopSize(len(lociPos), ch, dest)
        print 'Sample sizes are %d, %d, %d' % tuple(popSize)
        if ps[0] == 0:
            ps = popSize
        else:
            if ps[0] != popSize[0] or ps[1] != popSize[1] or ps[2] != popSize[2]:
                print "Population size does not match across chromosomes"
                sys.exit(1)
        print 'Creating population %s' % popFile
        pop = population(size=popSize, ploidy=2, loci=[len(lociPos)],
            lociPos=lociPos, lociNames=lociName, chromNames=[str(ch)])
        print 'Loading CEU population',
        load_population(pop, ch, type='CEU', dest=dest)
        print 'Loading YRI population',
        load_population(pop, ch, type='YRI', dest=dest)
        print 'Loading JPT+CHB population',
        load_population(pop, ch, type='JPT+CHB', dest=dest)
        print 'Setting genetic distance...'
        set_map_dist(pop, ch, dest)
        print 'Calculating allele frequency ...'
        Stat(pop, alleleFreq=range(pop.totNumLoci()))
        print "Saving population to %s..." % popFile
        SavePopulation(pop, os.path.join(dest, popFile))
        

if __name__ == '__main__':
    loadHapMap(range(1, 23))
