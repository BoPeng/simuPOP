#!/usr/bin/env python

'''
This script loads hapmap data in 
http://www.hapmap.org/downloads/phasing/2006-07_phaseII/consensus/
I chose consensus data because these SNPs exist in all three populations,
and I chose phased data because this information is important for the 
followup applications.

How to prepare the data:

1. Download all data in http://www.hapmap.org/downloads/phasing/2006-07_phaseII/consensus/
   You can save this page as index.html, and run
      > wget -r -l 1 -i index.html --force-html
   if you can use wget

2. Decompress all .gz file. If you use tcsh, run commands like
      > for f (*.gz)
      >     gzip -d $f
      > end

3. Run this script in the directory. You will obtain files named hapmap_XX.bin,
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

import os, sys

# revision number
rev = 21
legend_freq_file = 'genotypes_chr%d_r%d_nr_fwd_consensus_legend+freq'
genotype_file =  'genotypes_chr%d_%s_r%d_nr_fwd_phased_consensus'

def getLoci(ch):
    '''Loci information is retrieved from _consensus_legend+freq files'''
    lociPos = []
    lociName = []
    file = legend_freq_file % (ch, rev)
    if not os.path.isfile(file):
        print 'File %s does not exist' % file
        sys.exit(1)
    for line in open(file).readlines():
        fields = line.split()
        lociName.append(fields[0])
        # translate pos from index to cM. This is tentative
        lociPos.append(float(fields[1])/1000000.)
    return (lociPos, lociName)


def getPopSize(numLoci, ch):
    '''Get population size of each subpopulations
      This function also checks if each line has desired number of alleles.
    '''
    count = [0]*3
    ll = 0
    ceu = genotype_file % (ch, 'CEU', rev)
    for line in open(ceu).readlines():
        if (ll == 0 and len(line.split()) != numLoci) or \
            (ll != 0 and len(line) != ll):
            print "Number of loci does not match in %s " % ceu
            print "Number of loci: %d, number of fields: %d" % (numLoci, len(line.split()))
            sys.exit(1)
        ll = len(line)
        count[0] += 1
    ll = 0
    yri = genotype_file % (ch, 'YRI', rev)
    for line in open(yri).readlines():
        if (ll == 0 and len(line.split()) != numLoci) or \
            (ll != 0 and len(line) != ll):
            print "Number of loci does not match in %s " % yri
            print "Number of loci: %d, number of fields: %d" % (numLoci, len(line.split()))
            sys.exit(1)
        ll = len(line)
        count[1] += 1
    ll = 0
    jpt_chb = genotype_file % (ch, 'JPT+CHB', rev)
    for line in open(jpt_chb).readlines():
        if (ll == 0 and len(line.split()) != numLoci) or \
            (ll != 0 and len(line) != ll):
            print "Number of loci does not match in %s " % jpt_chb
            print "Number of loci: %d, number of fields: %d" % (numLoci, len(line.split()))
            sys.exit(1)        
        ll = len(line)
        count[2] += 1
    return [x/2 for x in count]       
    

def load_population(pop, ch, type):
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
    print 'from %s...' % file
    for line_no,line in enumerate(open(file).readlines()):
        genotype = [int(x) for x in line.split()]
        ind = line_no / 2
        ploidy = line_no % 2
        ind = pop.individual(ind, subPop)
        for i,g in enumerate(genotype):
            # always chromosome 0, because each population has only one chromosome
            ind.setAllele(g, i, ploidy)


if __name__ == '__main__':
    ps = [0,0,0]
    for ch in range(1, 23):
        popFile = "hapmap_%d.bin" % ch
        print "Processing chromosome %d..." % ch
        (lociPos, lociName) = getLoci(ch)
        print "    Number of loci", len(lociPos)
        print "    Getting population size and verifying data file..."
        popSize = getPopSize(len(lociPos), ch)
        print "    Population size %s (CEU and YRI only counts parents)" % popSize
        if ps[0] == 0:
            ps = popSize
        else:
            if ps[0] != popSize[0] or ps[1] != popSize[1] or ps[2] != popSize[2]:
                print "Population size does not match across chromosomes"
                sys.exit(1)
        print "    Creating population"
        pop = population(subPop=popSize, ploidy=2, loci=[len(lociPos)],
            lociPos=lociPos, lociNames=lociName)
        print "    Loading CEU population...",
        load_population(pop, ch, type='CEU')
        print "    Loading YRI population...",
        load_population(pop, ch, type='YRI')
        print "    Loading JPT+CHB population...",
        load_population(pop, ch, type='JPT+CHB')
        print "    Calculating allele frequency..."
        Stat(pop, alleleFreq=range(pop.totNumLoci()))
        print "    Saving population to %s..." % popFile
        SavePopulation(pop, popFile)
        
    
