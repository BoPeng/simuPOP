#!/usr/bin/env python
#
# Purpose:
#     This python module provides several utility functions that save a simuPOP
#     Population in merlin formats,
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
#         Move functions SaveMerlinDatFile, SaveMerlinMapFile and
#         SaveMerlinPedFile from simuUtil.py to the online cookbook.
# 

from simuPOP import *

def saveMerlinDatFile(pop, output='', outputExpr='', loci=[], fields=[], outputAffection=False):
    '''Output a .dat file readable by merlin'''
    if type(pop) == type(''):
        pop = LoadPopulation(pop)
    if output != '':
        file = output
    elif outputExpr != '':
        file = eval(outputExpr, globals(), pop.vars())
    else:
        raise exceptions.ValueError, "Please specify output or outputExpr"
    # open data file and pedigree file to write.
    try:
        datOut = open(file + ".dat", "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file + ".dat to write."
    if loci == []:
        loci = range(0, pop.totNumLoci())
    if outputAffection:
        print >> datOut, 'A\taffection'
    for f in fields:
        print >> datOut, 'T\t%s' % f
    for marker in loci:
        print >> datOut, 'M\t%s' % pop.locusName(marker)
    datOut.close()



def saveMerlinMapFile(pop, output='', outputExpr='', loci=[]):
    '''Output a .map file readable by merlin'''
    if type(pop) == type(''):
        pop = LoadPopulation(pop)
    if output != '':
        file = output
    elif outputExpr != '':
        file = eval(outputExpr, globals(), pop.vars())
    else:
        raise exceptions.ValueError, "Please specify output or outputExpr"
    # open data file and pedigree file to write.
    try:
        mapOut = open(file + ".map", "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file + ".map to write."
    if loci == []:
        loci = range(0, pop.totNumLoci())
    print >> mapOut, 'CHROMOSOME MARKER POSITION'
    for marker in loci:
        print >> mapOut, '%d\t%s\t%f' % (pop.chromLocusPair(marker)[0] + 1,
            pop.locusName(marker), pop.locusPos(marker))
    mapOut.close()



def saveMerlinPedFile(pop, output='', outputExpr='', loci=[], fields=[], header=False,
    outputAffection=False, affectionCode=['U', 'A'], combine=None, shift=1, **kwargs):
    '''Output a .ped file readable by merlin'''
    if type(pop) == type(''):
        pop = LoadPopulation(pop)
    if output != '':
        file = output
    elif outputExpr != '':
        file = eval(outputExpr, globals(), pop.vars())
    else:
        raise exceptions.ValueError, "Please specify output or outputExpr"
    # open data file and pedigree file to write.
    try:
        pedOut = open(file + ".ped", "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file + ".ped to write."
    if loci == []:
        loci = range(0, pop.totNumLoci())
    if header:
        print >> pedOut, "famID ID fa mo sex",
        if outputAffection:
            print >> pedOut, "affection",
        for f in fields:
            print >> pedOut, f,
        for marker in loci:
            print >> pedOut, pop.locusName(marker),
        print >> pedOut
    sexCode = {Male:1, Female:2}
    affectedCode = {False: affectionCode[0], True: affectionCode[1]}
    #
    pldy = pop.ploidy()
    def writeInd(ind, famID, id, fa, mo):
        print >> pedOut, '%d %d %d %d %d' % (famID, id, fa, mo, sexCode[ind.sex()]),
        if outputAffection:
            print >> pedOut, affectedCode[ind.affected()],
        for f in fields:
            print >> pedOut, '%.3f' % ind.info(f),
        if combine is None:
            # for efficiency, assuming diploid population
            for marker in loci:
                print >> pedOut, "%d %d" % (ind.allele(marker, 0) + shift, ind.allele(marker, 1) + shift),
        else:
            for marker in loci:
                print >> pedOut, "%d" % combine([ind.allele(marker, 0), ind.allele(marker, 1)]),
        print >> pedOut
    #
    # number of Pedigrees
    # get unique pedgree id numbers
    from sets import Set
    if pop.hasInfoField('pedindex') and pop.ancestralDepth() >= 1:
        peds = Set(pop.indInfo('pedindex'))
        # do not count peds
        peds.discard(-1)
        if len(peds) == 0:
            print 'Warning: no valid Pedigree, please check the pedindex field of your population'
        #
        newPedIdx = 1
        #
        for ped in peds:
            id = 1
            # -1 means no parents
            pastmap = {-1:0}
            # go from generation 2, 1, 0 (for example)
            for anc in range(pop.ancestralDepth(), -1, -1):
                newmap = {-1:0}
                pop.useAncestralPop(anc)
                # find all individual in this pedigree
                for i in range(pop.popSize()):
                    ind = pop.individual(i)
                    if ind.info('pedindex') == ped:
                        dad = int(ind.info('father_idx'))
                        mom = int(ind.info('mother_idx'))
                        if dad == mom and dad != -1:
                            print "Something wrong with pedigree %d, father and mother idx are the same: %s" % \
                                (ped, dad)
                        writeInd(ind, newPedIdx, id, pastmap.setdefault(dad, 0), pastmap.setdefault(mom, 0))
                        newmap[i] = id
                        id += 1
                pastmap = newmap
            newPedIdx += 1
    else:
        # rare case: no pedigree structure, only output the last generation without parents
        for idx, ind in enumerate(pop.individuals()):
            writeInd(ind, idx+1, 1, 0, 0)
    pedOut.close()


# save in merlin qtdt format
def saveQTDT(pop, output='', outputExpr='', loci=[], header=False,
    affectionCode=['U', 'A'], fields=[], combine=None, shift=1, **kwargs):
    """
    save population in Merlin/QTDT format. The population must have pedindex,
    father_idx and mother_idx information fields.

    pop
        Population to be saved. If pop is a filename, it will be loaded.

    output
        base filename.

    outputExpr
        expression for base filename, will be evaluated in pop's
        local namespace.

    affectionCode
        code for unaffected and affected. '1', '2' are default,
        but 'U', and 'A' or others can be specified.

    loci
        loci to output

    header
        whether or not put head line in the ped file.

    fields
        information fields to output

    combine
        an optional function to combine two alleles of a diploid
        Individual.

    shift
        if combine is not given, output two alleles directly, adding
        this value (default to 1).
    """
    if type(pop) == type(''):
        pop = LoadPopulation(pop)
    #
    if loci == []:
        loci = range(0, pop.totNumLoci())
    # write dat file
    if 'affection' in fields:
        outputAffection = True
        fields.remove('affection')
    else:
        outputAffection = False
    saveMerlinDatFile(pop, output, outputExpr, loci, fields, outputAffection)
    # write map file
    saveMerlinMapFile(pop, output, outputExpr, loci)
    # write ped file
    saveMerlinPedFile(pop, output, outputExpr, loci, fields, header,
        outputAffection, affectionCode, combine, shift)


