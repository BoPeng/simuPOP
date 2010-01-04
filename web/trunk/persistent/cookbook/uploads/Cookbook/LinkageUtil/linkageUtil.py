#!/usr/bin/env python
#
# Purpose:
#     This python module provides several utility functions that save a simuPOP
#     Population in linkage formats,
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
#         Move functions SaveLinkageFile, and LoadLinakgeFile from simuUtil.py
#         to the online cookbook.
# 

from simuPOP import *

def saveLinkage(pop, output='', outputExpr='', loci=[], shift=1, combine=None,
        fields = [], recombination=0.00001, penetrance=[0,0.25,0.5],
        affectionCode=['1', '2'],  pre=True, daf=0.001):
    """
    save population in Linkage format. Currently only
    support affected sibpairs sampled with affectedSibpairSample
    operator.

    pop
        Population to be saved. Must have ancestralDepth 1.
        paired individuals are sibs. Parental population are corresponding
        parents. If pop is a filename, it will be loaded.

    output
        Output.dat and output.ped will be the data and pedigree file.
        You may need to rename them to be analyzed by LINKAGE. This allows
        saving multiple files.

    outputExpr
        expression version of output.

    affectionCode
        default to '1': unaffected, '2': affected

    pre
        True. pedigree format to be fed to makeped. Non-pre format it is likely to
        be wrong now for non-sibpair families.

    Note
        the first child is always the proband.

    """
    if type(pop) == type(''):
        pop = LoadPopulation(pop)
    if output != '':
        file = output
    elif outputExpr != '':
        file = eval(outputExpr, globals(), pop.vars() )
    else:
        raise exceptions.ValueError, "Please specify output or outputExpr"
    # open data file and pedigree file to write.
    try:
        datOut = open(file + ".dat", "w")
        if pre:
            pedOut = open(file + ".pre", "w")
        else:
            pedOut = open(file + ".ped", "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file + ".dat/.ped to write."
    #
    if loci == []:
        loci = range(pop.totNumLoci())
    #
    # file is opened.
    # write data file
    # nlocus
    # another one is affection status
    # risklocus (not sure. risk is not to be calculated)
    # sexlink autosomal: 0
    # nprogram whatever
    # mutsys: all loci are mutational? 0 right now
    # mutmale
    # mutfemale
    # disequil: assume in LD? Yes.
    datOut.write( '''%d 0 0 5 << nlocus, risklocus, sexlink, nprogram
0 0 0 0 << mutsys, mutmale, mutfemale, disequil
'''    % (len(loci)+1) )
    # order of loci, allegro does not welcome comments after this line.
    # we need one more than the number of loci (including disease marker)
    datOut.write( ' '.join( [str(m+1) for m in range(len(loci) + 1)]) + "\n")
    # describe affected status
    datOut.write( "1 2 << affection status code, number of alleles\n")
    datOut.write( "%f %f << gene frequency\n" % ( 1-daf, daf) )
    datOut.write( "1 << number of factors\n")
    datOut.write( "%f %f %f << penetrance\n" % tuple(penetrance) )
    # describe each locus
    stat(pop, alleleFreq=loci)
    af = pop.dvars().alleleFreq
    for marker in loci:
        # now, 3 for numbered alleles
        numAllele = len(af[marker])
        print >> datOut, '3 %d << %s' % (numAllele, pop.locusName(marker))
        datOut.write( ''.join(['%.6f ' % af[marker][ale] for ale in range(numAllele)]) + ' << gene frequencies\n' )
    # sex-difference
    # interference
    datOut.write('0 0 << sex difference, interference\n')
    # recombination
    if type(recombination) in [type([]), type(())]:
        datOut.write( ' '.join(['%f '% x for x in recombination]) + ' << recombination rates \n ')
    else:
        datOut.write( ''.join(['%f '%recombination]*len(loci)) + ' << recombination rates \n ')
    # I do not know what they are
    datOut.write( "1 0.1 0.1\n")
    # done!
    datOut.close()
    # write pedigree file (affected sibpairs)
    # sex: in linkage, male is 1, female is 2
    sexCode = {Male:1, Female:2}
    affectedCode = {False: affectionCode[0], True: affectionCode[1]}
    # alleles string, since simuPOP allele starts from 0, add 1 to avoid
    # being treated as missing data.
    pldy = pop.ploidy()
    def writeInd(ind, famID, id, fa, mo):
        if pre:
            print >> pedOut, '%d %d %d %d %s %s' % (famID, id, fa, mo, sexCode[ind.sex()], affectedCode[ind.affected()]),
        else:
            if fa == 0:
                print >> pedOut, '%d %d %d 3 0 0 %d %s 0 %s' % (famID, id, fa, mo, sexCode[ind.sex()], affectedCode[ind.affected()]),
            else:
                print >> pedOut, '%d %d %d 0 4 4 %d %s 1 %s' % (famID, id, fa, mo, sexCode[ind.sex()], affectedCode[ind.affected()]),
        for marker in loci:
            if combine is None:
                for p in range(pldy):
                    print >> pedOut, " %d" % (ind.allele(marker, p) + shift),
            else:
                print >> pedOut, " %d" % combine([ind.allele(marker, p) for p in range(pldy)]),
        print >> pedOut
    #
    # get unique pedgree id numbers
    from sets import Set
    peds = Set(pop.indInfo('pedindex'))
    # do not count peds
    peds.discard(-1)
    if len(peds) == 0:
        print 'Warning: no valid Pedigree, please check the pedindex field of your population'
    #
    newPedIdx = 1
    for ped in peds:
        id = 1
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
    pedOut.close()


