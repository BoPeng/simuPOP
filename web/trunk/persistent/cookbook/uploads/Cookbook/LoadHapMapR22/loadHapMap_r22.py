#!/usr/bin/env python

'''
This script downloads and loads release 22 of hapmap datasets in 
http://www.hapmap.org/downloads/phasing/2007-08_rel22/phased.

Please refer to 
http://simupop.sourceforge.net/cookbook/pmwiki.php/Cookbook/LoadHapMap22 for 
for details about this script.
'''

from simuOpt import *
setOptions(optimized=True, alleleType='binary')
from simuPOP import *
from simuUtil import simuProgress

import os, sys, urllib, gzip, exceptions

# URL and revision number. You can choose to download other hapmap files
release = 22
# This link does not have JPT+CHT FILES, they are older (2007-11-20)
#Genotype_URL = 'http://www.hapmap.org/downloads/phasing/2007-08_rel22/phased/'
# Data from this link may not be correct but they are newer (2008-6-25)
Genotype_URL = 'http://ftp.hapmap.org/phasing/2007-08_rel22/phased/'
Recom_URL = 'http://ftp.hapmap.org/recombination/2008-03_rel22_B36/rates/'
legend_file = 'genotypes_chr%d_%s_r22_nr.b36_fwd_legend.txt.gz'
genotype_file = {
    'CEU': 'genotypes_chr%d_CEU_r22_nr.b36_fwd.phase.gz',
    'YRI': 'genotypes_chr%d_YRI_r22_nr.b36_fwd.phase.gz',
    'JPT+CHB': 'genotypes_chr%d_JPT+CHB_r22_nr.b36_fwd.phased.gz'  # Come on, what is going on?
}
recom_file = 'genetic_map_chr%s_b36.txt'

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


def getLoci(ch, sample, dest):
    '''Loci information is retrieved from legend files'''
    lociPos = []
    lociName = []
    file = legend_file % (ch, sample)
    downloadIfNeeded(Genotype_URL, dest, file)
    try:
        legend = gzip.open(os.path.join(dest, file))
        legend.readline()  # skip first line
        for line in legend.readlines():
            fields = line.split()
            lociName.append(fields[0])
            lociPos.append(float(fields[1]))
    except Exception, e:
        print 'Failed to read file ', file
        print 'You may want to remove your local copy and let this script re-download it.'
        print
        raise e
    return (lociPos, lociName)


def getPopSize(numLoci, ch, sample, dest):
    '''Get population size of a sample
      This function also checks if each line has desired number of alleles.
    '''
    count = 0
    ll = 0
    genotype = genotype_file[sample] % ch
    downloadIfNeeded(Genotype_URL, dest, genotype)
    for line in gzip.open(os.path.join(dest, genotype)).readlines():
        if (ll == 0 and len(line.split()) != numLoci) or \
            (ll != 0 and len(line) != ll):
            print "Number of loci does not match in %s " % ceu
            print "Number of loci: %d, number of fields: %d" % (numLoci, len(line.split()))
            sys.exit(1)
        ll = len(line)
        count += 1
    return count/2
    

def load_population(pop, ch, sample, dest):
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
    downloadIfNeeded(Genotype_URL, dest, file)
    progress = simuProgress('Load genotype from %s ' % file, pop.popSize())
    for line_no,line in enumerate(gzip.open(os.path.join(dest, file)).readlines()):
        genotype = [int(x) for x in line.split()]
        ind = line_no / 2
        ploidy = line_no % 2
        # always chromosome 0, because each population has only one chromosome
        pop.individual(ind).setGenotype(genotype, ploidy)
        progress.update(ind + 1)


def set_map_dist(pop, ch, dest):
    '''Set map distance for each locus'''
    file = recom_file % ch
    downloadIfNeeded(Recom_URL, dest, file)
    print 'Reading genetic map file %s...' % file,
    dist = {}
    for line in open(os.path.join(dest, file)).readlines():
        try:
            fields = line.split()
            pos = int(fields[0])
            dist[pos] = float(fields[2])
        except:
            pass
    print ' map distance of %d markers are found' % len(dist)
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
        progress = simuProgress('Estimating map distance of unspecified loci', totNumLoci)
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
            progress.update(loc)
        progress.done()
    print 'Map distance of %d markers (%.2f%% of %d) are estimated' % (totNumLoci - cnt,
        (totNumLoci - cnt) * 100.0/totNumLoci, totNumLoci)
    map = {}
    for loc in range(totNumLoci):
        map[pop.locusName(loc)] = map_dist[loc]
    pop.dvars().genDist = map

    
def loadHapMap(ch, sample, dest='.'):
    '''Download, import and save hapmap data of given chromosomes'''
    ps = None
    popFile = os.path.join(dest, "HapMap_%s_chr%d.pop" % (sample, ch))
    if os.path.isfile(popFile):
        print "Population %s already exists. Please remove it first if you would like to regenerate this file." % popFile
        return
    print "\n\nLoading HapMap chromosome %d for population %s" % (ch, sample)
    (lociPos, lociName) = getLoci(ch, sample, dest)
    print "%d loci (%.2f - %.2f cM) are located" % (len(lociPos), lociPos[0], lociPos[-1])
    popSize = getPopSize(len(lociPos), ch, sample, dest)
    print 'Sample size is %d' % popSize
    if ps is None:
        ps = popSize
    elif ps != popSize:
        print "Population size does not match across chromosomes"
        sys.exit(1)
    print 'Creating population %s' % popFile
    pop = population(size=popSize, ploidy=2, loci=[len(lociPos)],
        lociPos=lociPos, lociNames=lociName, chromNames=[str(ch)],
        subPopNames=[sample])
    load_population(pop, ch, sample, dest=dest)
    set_map_dist(pop, ch, dest)
    print "Saving population to %s..." % popFile
    pop.dvars().HapMap_rel = release
    pop.save(popFile)
        

options = [
    {'longarg': 'dest=',
     'default': '.',
     'useDefault': True,
     'label': 'Destination directory',
     'allowedTypes': [type('')],
     'validate': valueValidDir()
    },
    {'longarg': 'chroms=',
     'default': range(1, 23),
     'useDefault': True,
     'label': 'Chromosomes to download',
     'allowedTypes': [type([]), type(())],
     'chooseFrom': range(1, 23),
     'validate': valueListOf(valueBetween(1, 22))
    },
]

if __name__ == '__main__':
    pars = simuOpt(options, __doc__)
    if not pars.getParam():
        sys.exit(1)
    for chrom in pars.chroms:
        for sample in ['CEU', 'YRI', 'JPT+CHB']:
            loadHapMap(chrom, sample, pars.dest)
