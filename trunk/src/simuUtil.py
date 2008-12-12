#!/usr/bin/env python

############################################################################
#    Copyright (C) 2004 by Bo Peng
#    bpeng@mdanderson.org
#
#    $LastChangedDate$
#    $Rev$
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
#    GNU General Public License for more details.
#
#    You should havereceived a copy of the GNU General Public License
#    along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place - Suite 330, Boston, MA    02111-1307, USA.
############################################################################


"""
simuPOP utilities.

This module provides some commonly used operators
and format conversion utilities.

"""

import exceptions, operator, types, os, sys, getopt, re, math, tempfile, shutil
import copy, random

from simuPOP import *


# mating schemes

def cloneMating(numOffspring = 1., numOffspringFunc = None,
        maxNumOffspring= 1., mode = MATE_NumOffspring,
		sexParam = 0.5, sexMode = MATE_RandomSex, newSubPopSize = [],
		newSubPopSizeFunc = None, newSubPopSizeExpr = "", 
		subPop = None, weight = 0):
    '''
    Note that
   \li selection is not considered (fitness is ignored)
   \li sequentialParentMating is used. If offspring (virtual) subpopulation size
   is smaller than parental subpopulation size, not all parents will be cloned.
   If offspring (virtual) subpopulation size is larger, some parents will be
   cloned more than once.
   \li numOffspring interface is respected.
   \li during mating operators are applied.
    '''
    return pyMating(
        chooser = sequentialParentChooser(),
        generator = offspringGenerator(numOffspring,
            numOffspringFunc, maxNumOffspring, mode,
            sexParam, sexMode,
            cloneGenoTransmitter()),
        newSubPopSizeExpr = newSubPopSizeExpr,
        newSubPopSizeFunc = newSubPopSizeFunc,
        subPop = subPop,
        weight = weight)


def binomialSelection(numOffspring = 1., numOffspringFunc = None,
        maxNumOffspring= 1., mode = MATE_NumOffspring,
		sexParam = 0.5, sexMode = MATE_RandomSex, newSubPopSize = [],
		newSubPopSizeFunc = None, newSubPopSizeExpr = "", 
		subPop = None, weight = 0):
    '''a mating scheme that uses binomial selection, regardless of sex
   No sex information is involved (binomial random selection). Offspring is chosen from parental generation
   by random or according to the fitness values.
   In this mating scheme,
   \li \c numOffspring protocol is honored;
   \li population size changes are allowed;
   \li selection is possible;
   \li haploid population is allowed.
   <applicability>all ploidy</applicability>
    '''
    return pyMating(
        chooser = randomParentChooser(),
        generator = offspringGenerator(numOffspring,
            numOffspringFunc, maxNumOffspring, mode,
            sexParam, sexMode,
            cloneGenoTransmitter()),
        newSubPopSizeExpr = newSubPopSizeExpr,
        newSubPopSizeFunc = newSubPopSizeFunc,
        subPop = subPop,
        weight = weight)


def randomMating(numOffspring = 1., numOffspringFunc = None,
        maxNumOffspring= 1., mode = MATE_NumOffspring,
		sexParam = 0.5, sexMode = MATE_RandomSex, newSubPopSize = [],
		newSubPopSizeFunc = None, newSubPopSizeExpr = "", 
		subPop = None, weight = 0):
    '''
    A mating scheme of basic sexually random mating

   In this scheme, sex information is considered for each individual,
   and ploidy is always 2. Within each subpopulation, males and females
   are randomly chosen. Then randomly get one copy of chromosomes from
   father and mother. If only one sex exists in a subpopulation, a
   parameter (\c contWhenUniSex) can be set to determine the behavior.
   Default to continuing without warning.
    '''
    return pyMating(
        chooser = randomParentsChooser(true, false, Male, 1, Male, 0, ''),
        generator = offspringGenerator(numOffspring,
            numOffspringFunc, maxNumOffspring, mode,
            sexParam, sexMode,
            mendelianGenoTransmitter()),
        newSubPopSizeExpr = newSubPopSizeExpr,
        newSubPopSizeFunc = newSubPopSizeFunc,
        subPop = subPop,
        weight = weight)


def monogamousMating(replenish=False, numOffspring = 1., numOffspringFunc = None,
        maxNumOffspring= 1., mode = MATE_NumOffspring,
		sexParam = 0.5, sexMode = MATE_RandomSex, newSubPopSize = [],
		newSubPopSizeFunc = None, newSubPopSizeExpr = "", 
		subPop = None, weight = 0):
    '''
   This mating scheme is identical to random mating except that parents
   are chosen without replacement. Under this mating scheme, offspring share
   the same mother must share the same father. In case that all parental
   pairs are exhausted, parameter \c replenish=True allows for the replenishment
   of one or both sex groups.
    '''
    return pyMating(
        chooser = randomParentsChooser(false, replenish, Male, 1, Male, 0, ''),
        generator = offspringGenerator(numOffspring,
            numOffspringFunc, maxNumOffspring, mode,
            sexParam, sexMode,
            mendelianGenoTransmitter()),
        newSubPopSizeExpr = newSubPopSizeExpr,
        newSubPopSizeFunc = newSubPopSizeFunc,
        subPop = subPop,
        weight = weight)


def polygamousMating(polySex=Male, polyNum=1, replacement =False,
        replenish=False, numOffspring = 1., numOffspringFunc = None,
        maxNumOffspring= 1., mode = MATE_NumOffspring,
		sexParam = 0.5, sexMode = MATE_RandomSex, newSubPopSize = [],
		newSubPopSizeFunc = None, newSubPopSizeExpr = "", 
		subPop = None, weight = 0):
    '''
   This mating scheme is composed of a random parents chooser that allows for
   polygamous mating, and a mendelian offspring generator. In this mating scheme,
   a male (or female) parent will have more than one sex partner (\c numPartner).
   Parents returned from this parents chooser will yield the same male (or female)
   parents, each with varying partners.
    '''
    return pyMating(
        chooser = randomParentsChooser(replacement, replenish,
            polySex, polyNum, Male, 0, ''),
        generator = offspringGenerator(numOffspring,
            numOffspringFunc, maxNumOffspring, mode,
            sexParam, sexMode,
            mendelianGenoTransmitter()),
        newSubPopSizeExpr = newSubPopSizeExpr,
        newSubPopSizeFunc = newSubPopSizeFunc,
        subPop = subPop,
        weight = weight)


def alphaMating(alphaSex=Male, alphaNum=0, alphaFiels='',
        numOffspring = 1., numOffspringFunc = None,
        maxNumOffspring= 1., mode = MATE_NumOffspring,
		sexParam = 0.5, sexMode = MATE_RandomSex, newSubPopSize = [],
		newSubPopSizeFunc = None, newSubPopSizeExpr = "", 
		subPop = None, weight = 0):
    '''
     Only a number of alpha individuals can mate with individuals of opposite sex.

   This mating scheme is composed of an random parents chooser with alpha individuals,
   and a Mendelian offspring generator. That is to say, a certain number of alpha
   individual (male or female) are determined by \c alphaNum or an information field. Then,
   only these alpha individuals are able to mate with random individuals of
   opposite sex.

   	   \param alphaSex the sex of the alpha individual, i.e. alpha male
	           or alpha female who be the only mating individuals in their
	           sex group.
	   \param alphaNum Number of alpha individuals. If \c infoField is
	           not given, \c alphaNum random individuals with \c alphaSex
	           will be chosen. If selection is enabled, individuals with higher+  
               fitness values have higher probability to be selected. There is
	           by default no alpha individual (\c alphaNum = 0).
	   \param alphaField if an information field is given, individuals
	           with non-zero values at this information field are alpha individuals.
	           Note that these individuals must have \c alphaSex.

	   Please refer to class \c mating for descriptions of other parameters.
	   Note: If selection is enabled, it works regularly on on-alpha sex, but
	           works twice on alpha sex. That is to say, \c alphaNum alpha indiviudals
	           are chosen selectively, and selected again during mating.

    '''
    return pyMating(
        chooser = randomParentsChooser(True, False, Male, 1, alphaSex, alphaNum, alphaField),
        generator = offspringGenerator(numOffspring,
            numOffspringFunc, maxNumOffspring, mode,
            sexParam, sexMode,
            mendelianGenoTransmitter()),
        newSubPopSizeExpr = newSubPopSizeExpr,
        newSubPopSizeFunc = newSubPopSizeFunc,
        subPop = subPop,
        weight = weight)


def haplodiploidMating(alphaSex = Female, alphaNum = 1, alphaField = '',
		numOffspring = 1., numOffspringFunc = None, maxNumOffspring= 1., mode = MATE_NumOffspring,
		sexParam = 0.5, sexMode = MATE_RandomSex, newSubPopSize = [],
		newSubPopSizeFunc = None, newSubPopSizeExpr = "", 
		subPop = None, weight = 0):
    '''
    This mating scheme is composed of an alphaParentsChooser and a
    haplodiploidOffspringGenerator. The alphaParentChooser chooses a single
    Female randomly or from a given information field. This female will
    mate with random males from the colony. The offspring will have one of the
    two copies of chromosomes from the female parent, and the first copy
    of chromosomes from the male parent. Note that if a recombinator
    is used, it should disable recombination of male parent.
    '''
    return pyMating(
        chooser = randomParentsChooser(True, False, Male, 1,
            alphaSex, alphaNum, alphaField),
        generator = offspringGenerator(numOffspring,
            numOffspringFunc, maxNumOffspring, mode,
            sexParam, sexMode,
            haplodiploidGenoTransmitter()),
        newSubPopSizeExpr = newSubPopSizeExpr,
        newSubPopSizeFunc = newSubPopSizeFunc,
        subPop = subPop,
        weight = weight)


def selfMating(numOffspring = 1., numOffspringFunc = None,
        maxNumOffspring= 1., mode = MATE_NumOffspring,
		sexParam = 0.5, sexMode = MATE_RandomSex, newSubPopSize = [],
		newSubPopSizeFunc = None, newSubPopSizeExpr = "", 
		subPop = None, weight = 0):
    '''
    In this mating scheme, a parent is choosen randomly, acts
    both as father and mother in the usual random mating. The parent
    is chosen randomly, regardless of sex. If selection is turned on,
    the probability that an individual is chosen is proportional to
    his/her fitness.
    '''
    return pyMating(
        chooser = randomParentChooser(replacement=True, replenish=False),
        generator = offspringGenerator(numOffspring,
            numOffspringFunc, maxNumOffspring, mode,
            sexParam, sexMode,
            selfingGenoTransmitter()),
        newSubPopSizeExpr = newSubPopSizeExpr,
        newSubPopSizeFunc = newSubPopSizeFunc,
        subPop = subPop,
        weight = weight)


def consanguineousMatingMating(relativeFields = [], func = None, param = None,
        replacement = False, replenish = True,
        maxNumOffspring = 1., mode = MATE_NumOffspring,
		sexParam = 0.5, sexMode = MATE_RandomSex, newSubPopSize = [],
		newSubPopSizeFunc = None, newSubPopSizeExpr = "", 
		subPop = None, weight = 0):
    '''
       In this mating scheme, a parent is choosen randomly and mate with a
   relative that has been located and written to a number of information
   fields.

   	   This mating scheme randomly choose a parent and then choose his/her spouse from indexes
	   stored in \c infoFields.

	   \param relativeFields The information fields that stores indexes to other individuals
	    in a population. If more than one valid (positive value) indexes exist, a random
	    index will be chosen. (c.f. \c infoParentsChooser ) If there is no individual
	    having any valid index, the second parent will be chosen randomly from the
	    whole population.

	   \param func A python function that can be used to prepare the indexes of these
	    information fields. For example, functions population::locateRelatives and/or
	    population::setIndexesOfRelatives can be used to locate certain types of relatives
	    of each individual.

	   \param param An optional parameter that can be passed to \c func.

	   Please refer to \c infoParentsChooser and \c mendelianOffspringGenerator for
	   other parameters.
    '''
    # FIXME: lack a mechanism to call preparePopulation(pop)
    return pyMating(
        chooser = infoParentsChooser(relativeFields, replacement, replenish),
        generator = offspringGenerator(numOffspring,
            numOffspringFunc, maxNumOffspring, mode,
            sexParam, sexMode,
            mendelianGenoTransmitter()),
        newSubPopSizeExpr = newSubPopSizeExpr,
        newSubPopSizeFunc = newSubPopSizeFunc,
        subPop = subPop,
        weight = weight)



def pedigreeMating(ped, generator=None, newSubPopSize = [],
		newSubPopSizeFunc = None, newSubPopSizeExpr = "", 
		subPop = None, weight = 0):
    '''
//    In this scheme, a pedigree is given and the mating scheme will
//    choose parents and produce offspring strictly following the pedigree.
//    Parameters setting number of offspring per mating event, and
//    size of the offspring generations are ignored.
//
//    To implement this mating scheme in pyMating,
//    1.) a newSubPopSizeFunc should be given to return the exact subpopulation
//      size, returned from pedigree.subPopSizes(gen).
//    2.) use pedigreeChooser to choose parents
//    3.) use a suitable offspring generator to generate offspring.
//
//    This pedigreeMating helps you do 1 and 2, and use a mendelianOffspringGenerator
//    as the default offspring generator. You can use another offspring generator
//    by setting the generator parameter. Note that the offspring generator can
//    generate one and only one offspring each time.
    '''
    return pyMating(
        chooser = pedigreeParentsChooser(ped),
        generator = generator,
        newSubPopSizeExpr = newSubPopSizeExpr,
        newSubPopSizeFunc = newSubPopSizeFunc,
        subPop = subPop,
        weight = weight)



def _listVars(var, level=-1, name='', subPop=True, indent=0, curLevel=0):
    '''called by listVars. Will list variables recursively'''
    if type(var) == type( dw({}) ):
        var = var.__dict__
    # all level or level < specified maximum level
    if level < 0 or (level > 0 and curLevel < level):
        # list is list or typle type
        if type(var) == types.ListType or type(var) == types.TupleType:
            index = 0
            for x in var:
                # literals
                if type(x) != types.ListType and type(x) != types.DictType:
                    # this will save a huge amount of output for sparse matrix
                    # generated by Stat(LD=[]) etc.
                    if x != None:
                        if type(var) == types.ListType:
                            print ' '*indent, '['+str(index)+']\t', x
                        else:
                            print ' '*indent, '('+str(index)+')\t', x
                # nested stuff
                elif type(x) == types.ListType or type(x) == types.DictType:
                    if type(var) == types.ListType:
                        print ' '*indent, '['+str(index)+']\n',
                    else:
                        print ' '*indent, '('+str(index)+')\n',
                    _listVars(x, level, name, False, indent+2, curLevel + 1)
                index += 1
        elif type(var) == types.DictType:
            # none array first
            for x in var.items():
                if not type(x[1]) in [types.ListType, types.DictType, types.TupleType]:
                    if name == '' or x[0] == name:
                        print ' '*indent, x[0], ':\t', x[1]
            # array but not subPop
            for x in var.items():
                if x[0] != 'subPop' and type(x[1]) in [types.ListType, types.DictType, types.TupleType]:
                    if name == '' or x[0] == name:
                        print ' '*indent, x[0], ':\n',
                        _listVars(x[1], level, name, False, indent+2, curLevel + 1)
            # subPop
            if subPop == True and var.has_key('subPop'):
                print ' '*indent, 'subPop\n',
                _listVars(var['subPop'], level, name, False, indent+2, curLevel + 1)
        else:
            print ' '*indent, var
    else: # out of the range of level
        if type(var) == types.ListType or type(var) == types.TupleType:
            print ' '*indent, 'list of length', len(var)
        elif type(var) == types.DictType:
            print ' '*indent, 'dict with keys [',
            for num in range(0,len(var.keys())):
                if type(var.keys()[num]) == types.StringType:
                    print "'"+ var.keys()[num] + "',",
                else:
                    print var.keys()[num], ",",
                if num != len(var.keys())-1 and num%4 == 3:
                    print '\n' + ' '*(indent+5),
            print ']'
        else:
            print ' '*indent, var


def ListVars(var, level=-1, name='', subPop=True, useWxPython=True):
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

    useWxPython
        if True, use terminal output even if wxPython is available.
    '''
    if not useWxPython:
        _listVars(var, level, name, subPop, 0, 0)
        return

    # a wxPython version of listVars
    try:
        import wx, wx.py.filling as fill
    except:
        _listVars(var, level, name, subPop, 0, 0)
        return

    app = wx.App()
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
    app.SetTopWindow(fillFrame)
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
        m.append( [r/(n-1.)]*n)
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
        loci=[], combine=None, shift=1, **kwargs):
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
    # for each family
    def sexCode(ind):
        if ind.sex() == Male:
            return 1
        else:
            return 2
    # disease status: in linkage affected is 2, unaffected is 1
    def affectedCode(ind):
        if ind.affected():
            return 1
        else:
            return 2
    # write out header
    print >> out, 'id, ', ', '.join(fields), ', ',
    if combine is None:
        print >> out, ', '.join(['marker%s_1, marker%s_2' % (marker, marker) for marker in loci])
    else:
        print >> out, ', '.join(['marker%s' % marker for marker in loci])
    # write out
    id = 1
    pldy = pop.ploidy()
    for ind in pop.individuals():
        print >> out, id,
        for f in fields:
            if f == 'sex':
                print >> out, ', ', sexCode(ind),
            elif f == 'affection':
                print >> out, ', ', affectedCode(ind),
            else:
                print >> out, ', ', ind.info(f),
        for marker in loci:
            if combine is None:
                for p in range(pldy):
                    print >> out, ", %d" % (ind.allele(marker, p) + shift),
            else:
                print >> out, ", %d" % combine([ind.allele(marker, p) for p in range(pldy)]),
        print >> out
        id += 1
    out.close()


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
    def __init__(self, message, totalCount, progressChar='.', block=2, done=' Done.\n'):
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
        self.totalCount = totalCount
        self.percent = 0
        self.progressChar = progressChar
        self.block = block
        self.doneMsg = done
        self.completed = False
        sys.stdout.write(message)
        sys.stdout.flush()

    def update(self, count):
        '''
        Update the progreebar.
        '''
        #
        count = min(count, self.totalCount)
        #
        completed = int(round(100*count/self.totalCount))
        if completed <= self.percent:
            return
        for p in range(self.percent + 1, completed + 1):
            if p == 100:
                self.done()
            elif p % 10 == 0:
                sys.stdout.write(str(p/10))
            elif p % self.block == 0:
                sys.stdout.write(self.progressChar)
        self.percent = completed
        sys.stdout.flush()

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        if self.completed:
            return
        sys.stdout.write(self.doneMsg)
        sys.stdout.flush()
        self.completed = True


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


class _sample(pyOperator):
    '''
    Ascertainment/sampling refers to ways of selecting individuals from a
    population. This base class defines the common interface of all
    ascertainment operators, including how samples are saved and returned.
    Individual ascertainment operators (derived class) only need to
    write *prepareSample* and *drawSample* functions.
    '''
    def __init__(self, times = 1, name = '', nameExpr = '',
	       saveAs = '', saveAsExpr = '', *args, **kwargs):
        '''
        Create an operator that draws a certain type of samples from a
        population *times* times. The samples are saved in the population's
        local namespace if *name* or *nameExpr* is given, and are saved as
        diskfiles if *saveAs* or *saveAsExpr* is given. *nameExpr* or
        *saveAsExpr* are evaluated at the population's local namespace.
        '''
        self.times = times
        self.name = name
        self.nameExpr = nameExpr
        self.saveAs = saveAs
        self.saveAsExpr = saveAsExpr
        self.samples = []
        self.pedigree = None
        self.repr = '<simuPOP::sample>'
        pyOperator.__init__(self, func=self.drawSamples, *args, **kwargs)

    def prepareSample(self, pop):
        '''
        This function is usually used to prepare a pedigree object so that
        samples can be drawn.
        '''
        raise SystemError('Please re-implement this prepareSample function in the derived class.')
        return True

    def drawSample(self, pop):
        '''
        Draw and return a sample, using population *pop*, and *self.pedigree*
        prepared in prepareSample.
        '''
        raise SystemError('Please re-implement this drawSample function in the derived class.')
        return True

    def drawSamples(self, pop):
        if not self.prepareSample(pop) or self.times <= 0:
            return True

        self.samples = []
        for t in range(self.times):
            sample = self.drawSample(pop)
            self.samples.append(sample)
            # svae sample to local namespace
            if self.nameExpr != '':
                name = eval(self.nameExpr, globals(), pop.vars())
            elif self.name != '':
                name = self.name
            else:
                name = None
            if name is not None:
                if not pop.vars().has_key(name):
                    pop.dvars().name = []
                elif type(pop.vars()[name]) != type([]):
                    raise ValueError("Variable %s already exsits in population's local namespace." % name)
                pop.vars()[name].append(sample)
            # save to a file
            if self.saveAsExpr != '':
                saveAs = eval(self.saveExpr, globals(), pop.vars())
            elif self.saveAs != '':
                saveAs = self.saveAs
            else:
                saveAs = None
            if saveAs is not None:
                sample.save(saveAs)
        return True

    def __repr__(self):
        return self.repr

    def clone(self):
        return copy.copy(self)


class randomSample(_sample):
    '''
    This operator draws random individuals from a population repeatedly and
    forms a number of random samples. These samples can be put in the
    population's local namespace, or save to disk files. The function form
    of this operator returns a list of samples directly.
    '''
    def __init__(self, size, *args, **kwargs):
        '''
        Draw *size* random samples from a population *times* times. *size* can
        be a number or a list of numbers. In the former case, individuals are
        drawn from the whole population and the samples has only one
        subpopulation. In the latter case, a given number of individuals are
        drawn from each subpopulation and the result sample has the same number
        of subpopulation as the population from which samples are drawn. The
        samples are saved in the population's local namespace if *name* or
        *nameExpr* is given, and are saved as diskfiles if *saveAs* or
        *saveAsExpr* is given.
        '''
        _sample.__init__(self, *args, **kwargs)
        self.size = size
        self.repr = '<simuPOP::randomSample>'

    def prepareSample(self, pop):
        self.pedigree = pedigree(pop)
        self.pedigree.addInfoField('sample', -1)
        return True

    def drawSample(self, pop):
        if type(self.size) not in [type(()), type([])]:
            size = self.size
            if size > pop.popSize():
                print 'Warning: sample size %d is greater than population size %d.' % (size, pop.popSize())
                size = pop.popSize()
            # randomly choose self.size individuals
            values = [0] * size + [-1] * (pop.popSize() - size)
            random.shuffle(values)
            self.pedigree.setIndInfo(values, 'sample')
        else:
            for sp in range(pop.numSubPop()):
                size = self.size[sp]
                if size > pop.subPopSize(sp):
                    print 'Warning: sample size (%d) at subpopulation %d is greater than subpopulation size %d ' \
                        % (size, sp, pop.subPopSize(sp))
                values = [sp] * size + [-1] * (pop.subPopSize(sp) - size)
                random.shuffle(values)
                self.pedigree.setIndInfo(values, 'sample', sp)
        return pop.extract(field='sample', ped=self.pedigree)


def RandomSample(pop, *args, **kwargs):
     s = randomSample(*args, **kwargs)
     s.apply(pop)
     return s.samples
 
RandomSample.__doc__ = "Function version of operator randomSample whose __init__function is \n" + randomSample.__init__.__doc__


class caseControlSample(_sample):
    '''
    This operator chooses random cases and controls from a population
    repeatedly. These samples can be put in the population's local namespace,
    or save to disk files. The function form of this operator returns a list
    of samples directly.
    '''
    def __init__(self, cases, controls, *args, **kwargs):
        '''
        Draw *cases* affected and *controls* unaffected individuals from a
        population repeatedly. *cases* can be a number or a list of numbers.
        In the former case, affected individuals are drawn from the whole
        population. In the latter case, a given number of individuals are
        drawn from each subpopulation. The same hold for *controls*. The
        resulting samples have two subpopulations that hold cases and controls
        respectively. The samples are saved in the population's local namespace
        if *name* or *nameExpr* is given, and are saved as diskfiles if
        *saveAs* or *saveAsExpr* is given.
        '''
        _sample.__init__(self, *args, **kwargs)
        self.cases = cases
        self.controls = controls
        self.repr = '<simuPOP::caseControlSample>'

    def prepareSample(self, pop):
        self.pedigree = pedigree(pop)
        self.pedigree.addInfoField('sample', -1)
        self.pedigree.setVirtualSplitter(affectionSplitter())
        return True

    def drawSample(self, pop):
        if type(self.cases) not in [type(()), type([])] and type(self.controls) not in [type(()), type([])]:
            Stat(pop, numOfAffected=True)
            allCases = pop.dvars().numOfAffected
            allControls = pop.dvars().numOfUnaffected
            #
            cases = self.cases
            if cases > allCases:
                print 'Warning: number of cases %d is greater than number of affected individuals %d.' \
                    % (cases, allCases)
                cases = allCases
            #
            controls = self.controls
            if controls > allControls:
                print 'Warning: number of controls %d is greater than number of affected individuals %d.' \
                    % (controls, allControls)
                controls = allControls
            # caseControlly choose self.size individuals
            value_cases = [0] * cases + [-1] * (allCases - cases)
            value_controls = [1] * controls + [-1] * (allControls - controls)
            random.shuffle(value_cases)
            random.shuffle(value_controls)
            # assign information fields
            idx_cases = 0
            idx_controls = 0
            for sp in range(pop.numSubPop()):
                self.pedigree.setIndInfo(value_cases[idx_cases:], 'sample', (sp, 1))
                self.pedigree.setIndInfo(value_controls[idx_controls:], 'sample', (sp, 0))
                idx_cases += self.pedigree.subPopSize((sp, 1))
                idx_controls += self.pedigree.subPopSize((sp, 0))
        else:
            if len(self.cases) != pop.numSubPop():
                raise ValueError('If an list of cases is given, it should be specified for all subpopulations')
            if len(self.controls) != pop.numSubPop():
                raise ValueError('If an list of controls is given, it should be specified for all subpopulations')
            for sp in range(pop.numSubPop()):
                allCases = self.pedigree.subPopSize((sp, 1))
                allControls = self.pedigree.subPopSize((sp, 0))
                #
                cases = self.cases[sp]
                if cases > allCases:
                    print 'Warning: number of cases %d is greater than number of affected individuals %d in subpopulation %d.' \
                        % (cases, allCases, sp)
                    cases = allCases
                #
                controls = self.controls[sp]
                if controls > allControls:
                    print 'Warning: number of controls %d is greater than number of affected individuals %d in subpopulation %d.' \
                        % (controls, allControls, sp)
                    controls = allControls
                # 
                value_cases = [0] * cases + [-1] * (allCases - cases)
                value_controls = [1] * controls + [-1] * (allControls - controls)
                random.shuffle(value_cases)
                random.shuffle(value_controls)
                # assign information fields
                self.pedigree.setIndInfo(value_cases, 'sample', (sp, 1))
                self.pedigree.setIndInfo(value_controls, 'sample', (sp, 0))
        return pop.extract(field='sample', ped=self.pedigree)


def CaseControlSample(pop, *args, **kwargs):
     s = caseControlSample(*args, **kwargs)
     s.apply(pop)
     return s.samples
 
CaseControlSample.__doc__ = "Function version of operator caseControlSample whose __init__function is \n" + caseControlSample.__init__.__doc__


class affectedSibpairSample(_sample):
    '''
    This operator chooses affected sibpairs and their parents from a population
    repeatedly. These samples can be put in the population's local namespace,
    or save to disk files. The function form of this operator returns a list
    of samples directly.\n
    
    The population to be sampled needs to have at least one ancestral
    generation. In addition, parents of each offspring is needed so information
    fields, most likely *father_idx* and *mother_idx* should be used to track
    parents in the parental generation. An during mating operator
    *parentsTagger* is designed for such a purpose. In addition, because it is
    very unlikely for two random offspring to share parents, affected sibpairs
    can only be ascertained from populations that are generated using a mating
    scheme that produes more than one offspring at each mating event.
    '''
    def __init__(self, size, infoFields=['father_idx', 'mother_idx'], *args, **kwargs):
        '''
        Draw *size* families, including two affected siblings and their parents
        from a population repeatedly. The population to be sampled must have
        at least one ancestral generation. It should also have two information
        fields specified by parameter *infoFields* (Default to ``['father_idx',
        'mother_idx']``. Parameter *size* can be a number or a list of numbers.
        In the former case, affected sibpairs are drawn from the whole
        population. In the latter case, a given number of affected sibpairs are
        drawn from each subpopulation. In both cases, affected sibpairs in the
        resulting sample form their own subpopulations (of size two). The
        samples are saved in the population's local namespace if *name* or
        *nameExpr* is given, and are saved as diskfiles if *saveAs* or
        *saveAsExpr* is given.
        '''
        _sample.__init__(self, *args, **kwargs)
        self.size = size
        self.fields = infoFields
        if len(self.fields) != 2:
            raise ValueError('Two information fields that indicate indexes of parents in the parental generation is needed')
        self.repr = '<simuPOP::affectedSibpairSample>'

    def prepareSample(self, pop):
        if pop.ancestralGens() < 1:
            raise ValueError('No ancestral generation if found.')
        for field in self.fields:
            if field not in pop.infoFields():
                raise ValueError('Information field %s not found in population' % field)
        #
        self.pedigree = pedigree(pop, infoFields=self.fields, ancGen=1)
        self.pedigree.addInfoFields(['sample', 'pedindex', 'offspring0', 'offspring1', 'spouse'], -1)
        # locate all affected siblings
        self.pedigree.locateRelatives(REL_Offspring, ['offspring0', 'offspring1'])
        self.pedigree.locateRelatives(REL_Spouse, ['spouse'])
        # look for affected siblings from the parental generation
        self.pedigree.useAncestralGen(1)
        parent0 = self.pedigree.infoIdx(self.fields[0])
        parent1 = self.pedigree.infoIdx(self.fields[1])
        pedindex = self.pedigree.infoIdx('pedindex')
        offspring0 = self.pedigree.infoIdx('offspring0')
        offspring1 = self.pedigree.infoIdx('offspring1')
        spouse = self.pedigree.infoIdx('spouse')
        #
        pedCount = 0
        self.validPeds = [[] for x in range(self.pedigree.numSubPop())]
        for selfIdx, ind in enumerate(self.pedigree.individuals()):
            # if this individual is used
            # or if no valid spouse
            # or if no valid first offspring
            # or if no valid second offspring
            if ind.intInfo(pedindex) != -1 \
                or ind.intInfo(spouse) == -1 \
                or ind.intInfo(offspring0) == -1 \
                or ind.intInfo(offspring1) == -1:
                continue
            # if spouse has been used
            spouseIdx = ind.intInfo(spouse)
            spouseInd = self.pedigree.individual(spouseIdx)
            if spouseInd.intInfo(pedindex) != -1:
                continue
            # if the first offspring has been used, or if parents do not match, or if
            # not affected.
            offspring0Ind = self.pedigree.ancestor(ind.intInfo(offspring0), 0)
            if not offspring0Ind.affected() \
                or offspring0Ind.intInfo(pedindex) != -1 \
                or offspring0Ind.intInfo(parent0) not in [selfIdx, spouseIdx] \
                or offspring0Ind.intInfo(parent1) not in [selfIdx, spouseIdx]:
                continue
            # if the second offspring has been used, or if parents do not match, or
            # if not affected
            offspring1Ind = self.pedigree.ancestor(ind.intInfo(offspring1), 0)
            if not offspring1Ind.affected() \
                or offspring1Ind.intInfo(pedindex) != -1 \
                or offspring1Ind.intInfo(parent1) not in [selfIdx, spouseIdx] \
                or offspring1Ind.intInfo(parent1) not in [selfIdx, spouseIdx]:
                continue
            # good pedigree
            ind.setInfo(pedCount, pedindex)
            spouseInd.setInfo(pedCount, pedindex)
            offspring0Ind.setInfo(pedCount, pedindex)
            offspring1Ind.setInfo(pedCount, pedindex)
            # count the number of pedigrees
            self.validPeds[self.pedigree.subPopIndPair(selfIdx)[0]].append(pedCount)
            pedCount += 1
        return True

    def drawSample(self, pop):
        #
        pedindex = self.pedigree.infoIdx('pedindex')
        sample = self.pedigree.infoIdx('sample')
        #
        # clear information sample in case this operator is applied twice
        pop.setIndInfo([0], sample)
        #
        pedCount = sum([len(x) for x in self.validPeds])
        chosenPeds = [False] * pedCount
        #
        if type(self.size) not in [type(()), type([])]:
            size = self.size
            if size > pedCount:
                print 'Warning: number of requested sibpairs %d is greater than what exists (%d).' \
                    % (size, pedCount)
                size = pedCount
            #
            values = range(pedCount)
            random.shuffle(values)
            for v in values[:size]:
                chosenPeds[v] = True
        else:
            if len(self.size) != pop.numSubPop():
                raise ValueError('If an list of sizes is given, it should be specified for all subpopulations')
            for sp in range(pop.numSubPop()):
                allPeds = len(self.validPeds[sp])
                #
                size = self.size[sp]
                if size > allPeds:
                    print 'Warning: number of requested sibpairs %d is greater than what exists (%d) in subpopulation %d.' \
                        % (size, allPeds, sp)
                    size = allPeds
                #
                random.shuffle(self.validPeds[sp])
                for v in self.validPeds[sp][:size]:
                    chosenPeds[v] = True
        # assign genotype
        for gen in range(1, -1, -1):
            self.pedigree.useAncestralGen(gen)
            for ind in self.pedigree.individuals():
                ped = ind.intInfo(pedindex)
                if ped != -1 and chosenPeds[ped]:
                    ind.setInfo(ped, sample)
        sample = pop.extract(field='sample', ancGen=1, ped=self.pedigree)
        sample.removeEmptySubPops()
        return sample


def AffectedSibpairSample(pop, size, *args, **kwargs):
     s = affectedSibpairSample(size, *args, **kwargs)
     s.apply(pop)
     return s.samples
 
AffectedSibpairSample.__doc__ = "Function version of operator affectedSibpairSample whose __init__function is \n" + affectedSibpairSample.__init__.__doc__



# 
# 
# def new_caseControlSample(self, cases=[], controls=[], *args, **kwargs):
#     if type(cases) in [types.IntType, types.LongType]:
#         ca = [cases]
#         spSample = False
#     else:
#         ca = cases
#         spSample = True
#     if type(controls) in [types.IntType, types.LongType]:
#         ct = [controls]
#         spSample = False
#     else:
#         ct = controls
#         spSample = True
#     cppModule.caseControlSample_swiginit(self,
#         cppModule.new_caseControlSample(cases=ca, controls=ct,
#             spSample=spSample, *args, **kwargs))
# 
# new_caseControlSample.__doc__ = caseControlSample.__init__.__doc__
# del caseControlSample.__init__
# caseControlSample.__init__ = new_caseControlSample
# 
# 
# def new_affectedSibpairSample(self,size=[], *args, **kwargs):
#     if type(size) in [types.IntType, types.LongType]:
#         sz=[size]
#     else:
#         sz = size
#     cppModule.affectedSibpairSample_swiginit(self,
#         cppModule.new_affectedSibpairSample(size=sz, *args, **kwargs))
# 
# new_affectedSibpairSample.__doc__ = affectedSibpairSample.__init__.__doc__
# del affectedSibpairSample.__init__
# affectedSibpairSample.__init__ = new_affectedSibpairSample
# 
# 
# def new_largePedigreeSample(self, size=[], *args, **kwargs):
#     if type(size) in [types.IntType, types.LongType]:
#         sz= [size]
#     else:
#         sz = size
#     cppModule.largePedigreeSample_swiginit(self,
#         cppModule.new_largePedigreeSample(size=sz, *args, **kwargs))
# 
# new_largePedigreeSample.__doc__ = largePedigreeSample.__init__.__doc__
# del largePedigreeSample.__init__
# largePedigreeSample.__init__ = new_largePedigreeSample
# 
# 
# def new_nuclearFamilySample(self, size=[], *args, **kwargs):
#     if type(size) in [types.IntType, types.LongType]:
#         sz= [size]
#     else:
#         sz = size
#     cppModule.nuclearFamilySample_swiginit(self,
#         cppModule.new_nuclearFamilySample(size=sz, *args, **kwargs))
# 
# new_nuclearFamilySample.__doc__ = nuclearFamilySample.__init__.__doc__
# del nuclearFamilySample.__init__
# nuclearFamilySample.__init__ = new_nuclearFamilySample
# 


 
# def CaseControlSample(pop, *args, **kwargs):
#     s = caseControlSample(*args, **kwargs)
#     s.apply(pop)
#     return s.samples(pop)
# 
# if caseControlSample.__init__.__doc__ is not None:
#     CaseControlSample.__doc__ = "Function version of operator caseControlSample whose __init__function is \n" + caseControlSample.__init__.__doc__
# 
# def PySample(pop, *args, **kwargs):
#     s = pySample(*args, **kwargs)
#     s.apply(pop)
#     return s.samples(pop)
# 
# if pySample.__init__.__doc__ is not None:
#     PySample.__doc__ = "Function version of operator pySample whose __init__function is \n" + pySample.__init__.__doc__
# 
# def AffectedSibpairSample(pop, *args, **kwargs):
#     s = affectedSibpairSample(*args, **kwargs)
#     s.apply(pop)
#     return s.samples(pop)
# 
# if affectedSibpairSample.__init__.__doc__ is not None:
#     AffectedSibpairSample.__doc__ = "Function version of operator affectedSibpairSample whose __init__function is \n" + affectedSibpairSample.__init__.__doc__
# 
# def LargePedigreeSample(pop, *args, **kwargs):
#     s = largePedigreeSample(*args, **kwargs)
#     s.apply(pop)
#     return s.samples(pop)
# 
# if largePedigreeSample.__init__.__doc__ is not None:
#     LargePedigreeSample.__doc__ = "Function version of operator largePedigreeSample whose __init__function is \n" + largePedigreeSample.__init__.__doc__
# 
# def NuclearFamilySample(pop, *args, **kwargs):
#     s = nuclearFamilySample(*args, **kwargs)
#     s.apply(pop)
#     return s.samples(pop)
# 
# if nuclearFamilySample.__init__.__doc__ is not None:
#     NuclearFamilySample.__doc__ = "Function version of operator nuclearFamilySample whose __init__function is \n" + nuclearFamilySample.__init__.__doc__
# 
# def PySubset(pop, *args, **kwargs):
#     s = pySubset(*args, **kwargs)
#     s.apply(pop)
# 
# if pySubset.__init__.__doc__ is not None:
#     PySubset.__doc__ = "Function version of operator pySubset whose __init__function is \n" + pySubset.__init__.__doc__
# 
  
if __name__ == "__main__":
    pass
