#!/usr/bin/env python
from simuPOP import *
import os, sys, exceptions

def savePOPGENEDatFile(pop, output='', title='', pad=1):
    '''
    Output a data file readable by POPGENE
    pad - shift genotype from 0 to pad, 1 to 1+pad...
          since POPGENE genotype can not be 0 
    '''
    if type(pop) == type(''):
        pop = LoadPopulation(pop)
    # open data file to write
    try:
        datOut = open(output, "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file + " to write."
    # header section of POPGENE data file, default to output
    if title != '':
        print >> datOut, '/* ' + title + ' */'
    else:
        print >> datOut, '/* ' + output + ' */'
    print >> datOut, 'Number of populations = %d' % pop.numSubPop()
    print >> datOut, 'Number of loci = %d' % pop.totNumLoci()
    print >> datOut, 'Locus name :'
    # output locus names
    print >> datOut, ' '.join(pop.lociNames()) + '\n'
    # output genotype, only handle haploid and diploid
    if pop.ploidy() == 1:
        for sp in range(0, pop.numSubPop()):
            gt = ''
            # obtain genotype
            # POPGENE genotype can not be 0, so transfer 0 to pad, 1 to 1+pad...
            for g in pop.genotype(sp):
                if len(pop.alleleNames()) == 0:
                    gt += '%d' % (g + pad)
                else:
                    gt += pop.alleleName(g)
            # add space after each genotype
            gt = ' '.join(gt)
            # add enter after each individual
            gt = '\n'.join(gt[i:i+2*pop.totNumLoci()] for i in range(len(gt)) 
                if i % (2*pop.totNumLoci()) == 0)
            print >> datOut, gt
            print >> datOut, '\n'
    elif pop.ploidy() == 2:
        for sp in range(0, pop.numSubPop()):
            for ind in pop.individuals(sp):
                gt = ''
                for loc in range(pop.totNumLoci()):
                    if len(pop.alleleNames()) == 0:
                        gt += '%d%d' % (ind.allele(loc, 0) + pad,
                            ind.allele(loc, 1) + pad)
                    else:
                        gt += pop.alleleName(ind.allele(loc, 0)) + \
                            pop.alleleName(ind.allele(loc, 1))
                    gt += ' '
                print >> datOut, gt
            print >> datOut, '\n'
    else:
        raise exceptions.ValueError, "POPGENE only supports haploid and diploid."
    datOut.close()

if __name__ == '__main__':
    # for testing
    # haploid
    outfile = 'PopGeneOut_hap.dat'
    pop = Population(size=[2, 3], ploidy=1, loci=[2, 3])
    initByFreq(pop, [.3, .4, .3])
    savePOPGENEDatFile(pop, outfile, 'testing haploid')
    print open(outfile).read()    
    # diploid
    outfile = 'PopGeneOut_dip.dat'
    pop = Population(size=[4, 5], ploidy=2, loci=[2, 4], alleleNames=['A','B'])
    initByFreq(pop, [.3, .7])
    savePOPGENEDatFile(pop, outfile, 'testing diploid')
    print open(outfile).read()
