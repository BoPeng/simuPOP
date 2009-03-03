#!/usr/bin/env python

#
# $File: simuUtil.py $
# $LastChangedDate$
# $Rev$
#
# This file is part of simuPOP, a forward-time population genetics
# simulation environment. Please visit http://simupop.sourceforge.net
# for details.
#
# Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#


"""
simuPOP utilities.

This module provides some commonly used operators
and format conversion utilities.

"""

import exceptions, operator, types, os, sys, re

from simuPOP import *
import pprint
from simuOpt import simuOptions

def ViewVars(var, gui=None):
    '''
    list a variable in tree format, either in text format or in a
        wxPython window.

    var
        any variable to be viewed. Can be a dw object returned
        by dvars() function

    level
        level of display.

    name
        only view certain variable

    subPop
        whether or not display info in subPop

    '''
    if gui is None:
        gui = simuOptions['GUI']
    #
    if gui == False or gui == 'Tkinter':
        pprint.pprint(var)
        return

    try:
        import wx, wx.py.filling as fill
    except ImportError:
        pprint.pprint(var)
        return

    app = wx.PySimpleApp()
    wx.InitAllImageHandlers()
    if var==None:
        fillFrame = fill.FillingFrame()
    else:
        if type(var) == type( dw({}) ):
            fillFrame = fill.FillingFrame(rootObject=var.__dict__,
                rootLabel='var')
        else:
            fillFrame = fill.FillingFrame(rootObject=var,
                rootLabel='var')
    fillFrame.Show(True)
    app.MainLoop()


# migration rate matrix generators
def MigrIslandRates(r, n):
    '''migration rate matrix

::

         x m/(n-1) m/(n-1) ....
         m/(n-1) x ............
         .....
         .... m/(n-1) m/(n-1) x

    where x = 1-m
    '''
    # n==1?
    if n == 1:
        return [[1]]
    #
    m = []
    for i in range(0,n):
        m.append([r/(n-1.)]*n)
        m[-1][i] = 1-r
    return m


def MigrSteppingStoneRates(r, n, circular=False):
    '''migration rate matrix, circular stepping stone model (X=1-m)

::

           X   m/2               m/2
           m/2 X   m/2           0
           0   m/2 x   m/2 ......0
           ...
           m/2 0 ....       m/2  X

or non-circular

::

           X   m/2               m/2
           m/2 X   m/2           0
           0   m/2 X   m/2 ......0
           ...
           ...              m   X
    '''
    if n < 2:
        raise exceptions.ValueError("Can not define step stone model for n < 2")
    elif n == 2:
        return [[1-r,r],[r,1-r]]
    # the normal case (n>2)
    m = []
    for i in range(0, n):
        m.append([0]*n)
        m[i][i] = 1-r
        m[i][(i+1)%n] = r/2.
        m[i][(i+n-1)%n] = r/2.
    if not circular:
        m[0][1] = r
        m[0][-1] = 0
        m[n-1][0] = 0
        m[n-1][n-2] = r
    return m




#########################################################################
###
### The following are file import / export (mostly) functions
###
### These functions will observe the same interface for convenience
### some options are not needed, but should be provided (and safely
### ignored.)
###
### 1. pop: population to save, can be a string, in which case
###    the population will be loaded from a file.
### 2. output and outputExpr: output filename or pattern.
### 3. loci: loci to save, default to all loci
###    If you want to save all loci on a chromosome, use
###       loci = range(pop.chromBegin(ch), pop.chromEnd(ch))
### 4. shift: value add to allele number
### 5. combine: if combine alleles, function to use
### 6. fields: information fields to save
###
### X. additional parameters for each file format
###
###
#########################################################################

# save file in FSTAT format
def SaveFstat(pop, output='', outputExpr='', maxAllele=0, loci=[], shift=1,
    combine=None):
    if output != '':
        file = output
    elif outputExpr != '':
        file = eval(outputExpr, globals(), pop.vars() )
    else:
        raise exceptions.ValueError, "Please specify output or outputExpr"
    # open file
    try:
        f = open(file, "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file + " to write."
    #
    # file is opened.
    np = pop.numSubPop()
    if np > 200:
        print "Warning: Current version (2.93) of FSTAT can not handle more than 200 samples"
    if loci == []:
        loci = range(pop.totNumLoci())
    nl = len(loci)
    if nl > 100:
        print "Warning: Current version (2.93) of FSTAT can not handle more than 100 loci"
    if maxAllele != 0:
        nu = maxAllele
    else:
        nu = pop.maxAllele()
    if nu > 999:
        print "Warning: Current version (2.93) of FSTAT can not handle more than 999 alleles at each locus"
        print "If you used simuPOP_la library, you can specify maxAllele in population constructure"
    if nu < 10:
        nd = 1
    elif nu < 100:
        nd = 2
    elif nu < 1000:
        nd = 3
    else: # FSTAT can not handle this now. how many digits?
        nd = len(str(nu))
    # write the first line
    f.write( '%d %d %d %d\n' % (np, nl, nu, nd) )
    # following lines with loci name.
    for loc in loci:
        f.write( pop.locusName(loc) +"\n");
    gs = pop.totNumLoci()
    for sp in range(0, pop.numSubPop()):
        # genotype of subpopulation sp, individuals are
        # rearranged in perfect order
        gt = pop.arrGenotype(sp, True)
        for ind in range(0, pop.subPopSize(sp)):
            f.write("%d " % (sp+1))
            p1 = 2*gs*ind          # begining of first hemo copy
            p2 = 2*gs*ind + gs     # second
            for al in loci: # allele
                # change from 0 based allele to 1 based allele
                if combine is None:
                    ale1 = gt[p1+al] + shift
                    ale2 = gt[p2+al] + shift
                    f.write('%%0%dd%%0%dd ' % (nd, nd) % (ale1, ale2))
                else:
                    f.write('%%0%dd' % nd % combine([gt[p1+al], gt[p2+al]]))
            f.write( "\n")
    f.close()


def saveFstat(output='', outputExpr='', **kwargs):
    'operator version of the function SaveFstat'
    # deal with additional arguments
    parm = ''
    for (k,v) in kwargs.items():
        parm += str(k) + '=' + str(v) + ', '
    # pyEval( exposePop=1, param?, stmts="""
    # saveInFSTATFormat( pop, rep=rep?, output=output?, outputExpr=outputExpr?)
    # """)
    opt = '''pyEval(exposePop=1, %s
        stmts=r\'\'\'SaveFstat(pop, rep=rep, output=r"""%s""",
        outputExpr=r"""%s""" )\'\'\')''' % ( parm, output, outputExpr)
    # print opt
    return eval(opt)

# used to parse name
import re


def LoadFstat(file, loci=[]):
    '''load population from fstat file 'file'
    since fstat does not have chromosome structure
    an additional parameter can be given
    '''
    try:
        f = open(file, "r")
    except exceptions.IOError:
        raise exceptions.IOError("Can not open file " + file + " to read.")
    #
    # file is opened. get basic parameters
    try:
        # get numSubPop(), totNumLoci(), maxAllele(), digit
        [np, nl, nu, nd] = map(int, f.readline().split())
    except exceptions.ValueError:
        raise exceptions.ValueError("The first line does not have 4 numbers. Are you sure this is a FSTAT file?")

    # now, ignore nl lines, if loci is empty try to see if we have info here
    # following lines with loci name.
    numLoci = loci
    lociNames = []
    if loci != []: # ignore allele name lines
        if nl != reduce(operator.add, loci):
            raise exceptions.ValueError("Given number of loci does not add up to number of loci in the file")
        for al in range(0, nl):
            lociNames.append(f.readline().strip() )
    else:
        scan = re.compile(r'\D*(\d+)\D*(\d+)')
        for al in range(0, nl):
            lociNames.append( f.readline().strip())
            # try to parse the name ...
            try:
                #print "mating ", lociNames[-1]
                ch,loc = map(int, scan.match(lociNames[-1]).groups())
                # get numbers?
                #print ch, loc
                if len(numLoci)+1 == ch:
                    numLoci.append( loc )
                else:
                    numLoci[ ch-1 ] = loc
            except exceptions.Exception:
                pass
        # if we can not get numbers correct, put all loci in one chromosome
        if reduce(operator.add, numLoci, 0) != nl:
            numLoci = [nl]
    #
    # now, numLoci should be valid, we need number of population
    # and subpopulations
    maxAllele = 0
    gt = []
    for line in f.readlines():
        gt.append( line.split() )
    f.close()
    # subpop size?
    subPopIndex = map(lambda x:int(x[0]), gt)
    # count subpop.
    subPopSize = [0]*subPopIndex[-1]
    for i in range(0, subPopIndex[-1]):
        subPopSize[i] = subPopIndex.count(i+1)
    if len(subPopSize) != np:
        raise exceptions.ValueError("Number of subpop does not match")
    if reduce(operator.add, subPopSize) != len(gt):
        raise exceptions.ValueError("Population size does not match")
    # we have all the information, create a population
    pop = population( subPop=subPopSize, loci = numLoci, ploidy=2,
        lociNames=lociNames)
    #
    gs = pop.totNumLoci()
    popGT = pop.arrGenotype(True)
    for ind in range(0, len(gt)):
        p1 = 2*gs*ind                # begining of first hemo copy
        p2 = 2*gs*ind + gs     # second
        for al in range(0, gs): # allele
            ale = int(gt[ind][al+1])
            popGT[2*gs*ind + al] = ale/(10**nd)
            popGT[2*gs*ind + gs + al] = ale%(10*nd)
            if popGT[2*gs*ind + al] > maxAllele:
                maxAllele = popGT[2*gs*ind + al]
            if popGT[2*gs*ind + gs + al] > maxAllele:
                maxAllele = popGT[2*gs*ind + gs + al]
    pop.setMaxAllele(maxAllele)
    return pop


def LoadGCData(file, loci=[]):
    '''HIDDEN
    '''
    # open file
    try:
        f = open(file, "r")
    except exceptions.IOError:
        raise exceptions.IOError("Can not open file " + file + " to read.")
    gt = []
    for line in f.readlines():
        gt.append( line.split() )
    f.close()
    # now we have a 2-d matrix of strings
    # population size?
    popSize = len(gt)
    # number of alleles
    numAllele = (len(gt[0]))/2-1
    #
    # loci number
    if reduce(operator.add, loci,0.) == numAllele:
        lociNum = loci
    else:
        lociNum = [numAllele]
    # create population
    pop = population(size=popSize, ploidy=2, loci=lociNum, maxAllele=2)
    #
    gs = pop.totNumLoci()
    popGT = pop.arrGenotype(True)
    for ind in range(0, len(gt)):
        pop.individual(ind).setAffected( int(gt[ind][1]))
        p1 = 2*gs*ind                # begining of first hemo copy
        p2 = 2*gs*ind + gs     # second
        for al in range(0, gs): # allele
            popGT[2*gs*ind + al] = int(gt[ind][al*2+2])
            popGT[2*gs*ind + gs + al] = int(gt[ind][al*2+3])
    return pop

#
def SaveLinkage(pop, output='', outputExpr='', loci=[], shift=1, combine=None,
        fields = [], recombination=0.00001, penetrance=[0,0.25,0.5],
        affectionCode=['1', '2'],  pre=True, daf=0.001):
    """
    save population in Linkage format. Currently only
    support affected sibpairs sampled with affectedSibpairSample
    operator.

    pop
        population to be saved. Must have ancestralDepth 1.
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
    Stat(pop, alleleFreq=loci)
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
        print 'Warning: no valid pedigree, please check the pedindex field of your population'
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


# operator version of saveLinkage
def saveLinkage(output='', outputExpr='', **kwargs):
    "An operator to save population in linkage format"
    # deal with additional arguments
    parm = ''
    for (k,v) in kwargs.items():
        parm += str(k) + '=' + str(v) + ', '
    # pyEval( exposePop=1, param?, stmts="""
    # saveInFSTATFormat( pop, rep=rep?, output=output?, outputExpr=outputExpr?)
    # """)
    opt = '''pyEval(exposePop=1, %s
        stmts=r\'\'\'SaveLinkage(pop, rep=rep, output=r"""%s""",
        outputExpr=r"""%s""" )\'\'\')''' % ( parm, output, outputExpr)
    # print opt
    return eval(opt)


def SaveSolarFrqFile(pop, output='', outputExpr='', loci=[], calcFreq=True):
    '''Output a frequency file, in a format readable by solar
    calcFreq
        whether or not calculate allele frequency
    '''
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
        frqOut = open(file + ".frq", "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file + ".frq to write."
    if loci == []:
        loci = range(0, pop.totNumLoci())
    if calcFreq or not pop.vars().has_key('alleleFreq'):
        Stat(pop, alleleFreq=loci)
    alleleFreq = pop.dvars().alleleFreq
    for m in loci:
        try:
            print >> frqOut, pop.locusName(m),
            for a in range(len(alleleFreq[m])):
                print >> frqOut, '\t%d\t%f' % (a+1, alleleFreq[m][a]),
            print >> frqOut
        except:
            print "Can not output allele frequency for marker %s " % m
    frqOut.close()


def SaveMerlinDatFile(pop, output='', outputExpr='', loci=[], fields=[], outputAffection=False):
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


def SaveMerlinMapFile(pop, output='', outputExpr='', loci=[]):
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


def SaveMerlinPedFile(pop, output='', outputExpr='', loci=[], fields=[], header=False,
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
    # number of pedigrees
    # get unique pedgree id numbers
    from sets import Set
    if pop.hasInfoField('pedindex') and pop.ancestralDepth() >= 1:
        peds = Set(pop.indInfo('pedindex'))
        # do not count peds
        peds.discard(-1)
        if len(peds) == 0:
            print 'Warning: no valid pedigree, please check the pedindex field of your population'
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
def SaveQTDT(pop, output='', outputExpr='', loci=[], header=False,
    affectionCode=['U', 'A'], fields=[], combine=None, shift=1, **kwargs):
    """
    save population in Merlin/QTDT format. The population must have pedindex,
    father_idx and mother_idx information fields.

    pop
        population to be saved. If pop is a filename, it will be loaded.

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
        individual.

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
    SaveMerlinDatFile(pop, output, outputExpr, loci, fields, outputAffection)
    # write map file
    SaveMerlinMapFile(pop, output, outputExpr, loci)
    # write ped file
    SaveMerlinPedFile(pop, output, outputExpr, loci, fields, header,
        outputAffection, affectionCode, combine, shift)



def SaveCSV(pop, output='', outputExpr='', fields=['sex', 'affection'],
        loci=[], combine=None, shift=1,
        sexCode={Male: '1', Female: '2'}, affectionCode={True: '1', False: '2'}, **kwargs):
    """save file in CSV format

    fileds
        information fields, 'sex' and 'affection' are special fields that
        is treated differently.

    genotype
        list of loci to output, default to all.

    combine
        how to combine the markers. Default to None.
        A function can be specified, that takes the form::

             def func(markers):
                 return markers[0]+markers[1]

    shift
        since alleles in simuPOP is 0-based, shift=1 is usually needed to
        output alleles starting from allele 1. This parameter is ignored if
        combine is used.

    """
    if output != '':
        file = output
    elif outputExpr != '':
        file = eval(outputExpr, globals(), pop.vars() )
    else:
        raise exceptions.ValueError, "Please specify output or outputExpr"
    if loci == []:
        loci = range(0, pop.totNumLoci())
    try:
        out = open( file, "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file +" to write."
    # keep the content of pieces in strings first
    content = [''] * pop.numChrom()
    # write out header
    print >> out, 'id, ', ', '.join(fields), ', ',
    if combine is None:
        print >> out, ', '.join(['%s_1, %s_2' % (pop.locusName(loc), pop.locusName(loc)) for loc in loci])
    else:
        print >> out, ', '.join(['%s' % pop.locusName(loc) for loc in loci])
    # write out
    id = 1
    pldy = pop.ploidy()
    for ind in pop.individuals():
        print >> out, id,
        for f in fields:
            if f == 'sex':
                print >> out, ',', sexCode[ind.sex()],
            elif f == 'affection':
                print >> out, ',', affectionCode[ind.affected()],
            else:
                print >> out, ',', ind.info(f),
        for marker in loci:
            if combine is None:
                for p in range(pldy):
                    out.write(", %d" % (ind.allele(marker, p) + shift))
            else:
                out.write(", %d" % combine([ind.allele(marker, p) for p in range(pldy)]))
        print >> out
        id += 1
    out.close()


class _baseProgressBar:
    def __init__(self, message, totalCount):
        '''
        message
            Title of the progress bar

        totalCount
            Total expected steps.

        done
            Message displayed when the job is finished.
        '''
        self.message = message
        self.totalCount = totalCount
        self.percent = 0
        self.completed = False

    def update(self, count):
        '''
        Update the progreebar.
        '''
        count = min(count, self.totalCount)
        self.progress = int(round(100*count/self.totalCount))
        if self.progress <= self.percent:
            return False
        else:
            return True

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        if self.completed:
            return False
        else:
            self.completed = True
            return True

class _textProgressBar(_baseProgressBar):
    def __init__(self, message, totalCount, progressChar='.', block=2, done=' Done.\n'):
        '''
        message
            Title of the progress bar

        totalCount
            Total expected steps.

        progressChar
            Character to be displayed for each progress.

        block
            display progress at which interval (in terms of percentage)?

        done
            Message displayed when the job is finished.
        '''
        _baseProgressBar.__init__(self, message, totalCount)
        self.percent = 0
        self.progressChar = progressChar
        self.block = block
        self.doneMsg = done
        sys.stdout.write(message)
        sys.stdout.flush()

    def update(self, count):
        '''
        Update the text progreebar.
        '''
        if not _baseProgressBar.update(self, count):
            return
        for p in range(self.percent + 1, self.progress + 1):
            if p == 100:
                self.done()
            elif p % 10 == 0:
                sys.stdout.write(str(p/10))
            elif p % self.block == 0:
                sys.stdout.write(self.progressChar)
        sys.stdout.flush()
        self.percent = self.progress
        if self.percent == 100:
            self.done()

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        if not _baseProgressBar.done(self):
            return
        sys.stdout.write(self.doneMsg)
        sys.stdout.flush()


class _tkProgressBar(_baseProgressBar):
    def __init__(self, message, totalCount):
        '''
        totalCount
            Total expected steps.

        progressChar
            Character to be displayed for each progress.

        block
            display progress at which interval (in terms of percentage)?

        done
            Message displayed when the job is finished.
        '''
        _baseProgressBar.__init__(self, message, totalCount)
        import Tkinter as tk
        self.width = 300
        self.height = 30
        self.max = 100
        self.fillColor = 'blue'
        self.labelColor = 'black'
        self.label = 'Progress'
        #
        self.app = tk.Tk()
        self.app.title(self.label)
        self.frame = tk.Frame(self.app, bd=0)
        self.canvas = tk.Canvas(self.frame, bd=0, width=self.width+40,
            height = self.height + 70, highlightthickness=0)
        self.label = self.canvas.create_text(20, 20, 
            text='', anchor="w", fill=self.labelColor, font=('Verdana', 10))
        self.scale = self.canvas.create_rectangle(
            20, 50, self.width + 20, 50 + self.height, fill=self.fillColor)
        self.rect = self.canvas.create_rectangle(
            20, 50, self.width + 20, 50 + self.height)
        self.canvas.pack(side='top', fill='x', expand='yes', padx=0)
        self.update(0)
        self.frame.pack(padx=0, pady=0)

    def update(self, count):
        '''
        Update the progreebar.
        '''
        if not _baseProgressBar.update(self, count):
            return
        #
        self.canvas.coords(self.scale, 20, 50,
            20 + self.progress * 1.0 / self.max * self.width, 50 + self.height)
        # Now update the colors
        self.canvas.itemconfig(self.scale, fill=self.fillColor)
        self.canvas.itemconfig(self.label, fill=self.labelColor)
        # And update the label
        if self.progress > 0:
            self.canvas.itemconfig(self.label, text=self.message + "\n%d%% completed." % self.progress)
        else:
            self.canvas.itemconfig(self.label, text=self.message)
        self.canvas.update_idletasks()
        self.app.update()
        #
        self.percent = self.progress
        if self.percent == 100:
            self.done()

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        if not _baseProgressBar.done(self):
            return
        self.app.destroy()
        del self.app


class _wxProgressBar(_baseProgressBar):
    def __init__(self, message, totalCount):
        '''
        totalCount
            Total expected steps.

        progressChar
            Character to be displayed for each progress.

        block
            display progress at which interval (in terms of percentage)?

        done
            Message displayed when the job is finished.
        '''
        _baseProgressBar.__init__(self, message, totalCount)
        import wx
        self.app = wx.PySimpleApp(0)
        self.dialog = wx.ProgressDialog(
            'Progress', self.message + '\n', self.totalCount,
            style = \
                # wx.PD_CAN_ABORT | \
                # wx.PD_CAN_SKIP | \
                wx.PD_ELAPSED_TIME | \
                # wx.PD_ESTIMATED_TIME | \
                wx.PD_AUTO_HIDE | \
                wx.PD_REMAINING_TIME
            )
        self.dialog.Update(0)

    def update(self, count):
        '''
        Update the progreebar.
        '''
        if not _baseProgressBar.update(self, count):
            return
        self.dialog.Update(count, self.message + "\n%d%% completed." % self.progress)
        self.percent = self.progress
        if self.percent == 100:
            self.done()

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        if not _baseProgressBar.done(self):
            return
        self.dialog.Destroy()
        del self.app



class simuProgress:
    '''
    This class defines a very simple text based progress bar. It will display a
    character (default to "``.``") for each change of progress (default to 2%),
    and a number (1, 2, ..., 9) for each 10% of progress, and print a message
    (default to "``Done.\\n``") when the job is finished.

    This class is used as follows::

        progress = simuProgress("Start simulation", 500)
        for i in range(500):
            progress.update(i+1)
        # if you would like to make sure the done message is displayed.
        progress.done()
    '''
    def __init__(self, message, totalCount, progressChar='.', block=2, done=' Done.\n', gui=None):
        '''
        message
            Title of the progress bar.

        totalCount
            Total expected steps.

        progressChar
            Character to be displayed for each progress.

        block
            display progress at which interval (in terms of percentage)?

        done
            Message displayed when the job is finished.
        '''
        if gui is None:
            self.gui = simuOptions['GUI']
        else:
            self.gui = gui
        if self.gui in ['wxPython', True]:
            try:
                import wx
            except ImportError:
                self.gui = 'Tkinter'
        if self.gui == 'Tkinter':
            try:
                import Tkinter
            except ImportError:
                self.gui = False
        if self.gui == 'wxPython':
            self.progressBar = _wxProgressBar(message, totalCount)
        elif self.gui == 'Tkinter':
            self.progressBar = _tkProgressBar(message, totalCount)
        else:
            self.progressBar = _textProgressBar(message, totalCount, progressChar, block, done)

    def update(self, count):
        '''
        Update the progreebar.
        '''
        self.progressBar.update(count)

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        self.progressBar.done()


class pySubset(pyOperator):
    '''
    This operator rearranges and removes individuals according to their values at
    an information field. Individuals with positive values at this information
    field are moved to the subpopulation specified by the integer value of this
    value. Individuals with negative values are removed. There is no function
    form of this operator because this operator is essentially a wrapper around
    function ``population::setSubPopByIndInfo(field)``.
    '''
    def __init__(self, field, *args, **kwargs):
        '''
        Create a pySubset operator that rearranges and removes individuals
        according to their values at an information field *field*.
        '''
        self.field = field
        pyOperator.__init__(self, func=self.apply,
            param=self.field, *args, **kwargs)

    def apply(self, pop, field):
        pop.setSubPopByIndInfo(field)
        return True

    def __repr_():
        return "<simuPOP::pySubset>"

    def clone():
        return pySubset(self.field)

  
if __name__ == "__main__":
    pass
