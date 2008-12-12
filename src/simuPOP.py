#!/usr/bin/env python

############################################################################
#    Copyright (C) 2004 by Bo Peng
#    bpeng@rice.edu
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
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
############################################################################


# get options
from simuOpt import simuOptions
import os, sys

if simuOptions['Optimized']:
    if simuOptions['AlleleType'] == 'short':
        from simuPOP_op import *
    elif simuOptions['AlleleType'] == 'long':
        from simuPOP_laop import *
    elif simuOptions['AlleleType'] == 'binary':
        from simuPOP_baop import *
    else:
        from simuPOP_op import *
else:
    if simuOptions['AlleleType'] == 'short':
        from simuPOP_std import *
    elif simuOptions['AlleleType'] == 'long':
        from simuPOP_la import *
    elif simuOptions['AlleleType'] == 'binary':
        from simuPOP_ba import *
    else:
        from simuPOP_std import *


if not simuOptions['Quiet']:
    print "simuPOP : Copyright (c) 2004-2008 Bo Peng"
    # compile date, compiler etc are macros that are replaced during compile time.
    if simuVer() == '9.9.9':
        # this is the subversion version of simuPOP
        print ("Developmental Version (%s) for Python %s" % (ModuleDate(), ModulePyVersion() ))
    else:
        # this is the released version
        print ("Version %s (Revision %d, %s) for Python %s" % (simuVer(), simuRev(), ModuleDate(),
            ModulePyVersion() ))
    print ModuleCompiler()
    print "Random Number Generator is set to %s with random seed 0x%08x" % (rng().name(), rng().seed())
    # MaxAllele + 1 since 0 is one of the allelic states
    if Optimized():
        print "This is the optimized %s allele version with %d maximum allelic states." % (AlleleType(), MaxAllele()+1)
    else:
        print "This is the standard %s allele version with %d maximum allelic states." % (AlleleType(), MaxAllele()+1)
    print "For more information, please visit http://simupop.sourceforge.net,"
    print "or email simupop-list@lists.sourceforge.net (subscription required)."

    if simuOptions['Debug'] != []:
        for g in simuOptions['Debug']:
            if g not in ['', None]:
                print "Turn on debug '%s'" % g
                TurnOnDebug(g)

# Other definitions that does not really belong to simuUtil.py
#

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



