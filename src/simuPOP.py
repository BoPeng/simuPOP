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

"""
simuPOP core module.
"""

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
    # Turn on general debug information when not in 'quiet' mode
    # This will print out error messages when needed.
    TurnOnDebug('DBG_GENERAL')

if simuOptions['Debug'] != []:
    for g in simuOptions['Debug']:
        if g not in ['', None]:
            print "Turn on debug '%s'" % g
            TurnOnDebug(g)

# Other definitions that does not really belong to simuUtil.py
#


class dw(object):
    def __init__(self, var):
        try:
            self.__dict__ = var
        except exceptions.TypeError:
            raise exceptions.TypeError("The returned value is not a dictionary.\nNote: simu.vars() is a list so simu.dvars() is not allowed. \n    Use simu.dvars(rep) for population namespace.")
    def clear(self):
        self.__dict__.clear()
    def __repr__(self):
        return str(self.__dict__)

def dvars(self, *args, **kwargs):
    return dw(self.vars(*args, **kwargs))

population.dvars = dvars
simulator.dvars = dvars

# expose the clone() method to Python copy module.
def deepcopy(self, memo):
    return self.clone()

population.__deepcopy__ = deepcopy
simulator.__deepcopy__ = deepcopy
baseOperator.__deepcopy__ = deepcopy

#
# functions to corresponding operators
def Dump(pop, *args, **kwargs):
    dumper(*args, **kwargs).apply(pop)

if dumper.__init__.__doc__ is not None:
    Dump.__doc__ = "Function version of operator dump whose __init__ function is \n" + dumper.__init__.__doc__

def InitSex(pop, *args, **kwargs):
    initSex(*args, **kwargs).apply(pop)

if initSex.__init__.__doc__ is not None:
    InitSex.__doc__ = "Function version of operator initSex whose __init__ function is \n" + initSex.__init__.__doc__

def InitByFreq(pop, *args, **kwargs):
    initByFreq(*args, **kwargs).apply(pop)

if initByFreq.__init__.__doc__ is not None:
    InitByFreq.__doc__ = "Function version of operator initByFreq whose __init__ function is \n" + initByFreq.__init__.__doc__

def InitByValue(pop, *args, **kwargs):
    initByValue(*args, **kwargs).apply(pop)

if initByValue.__init__.__doc__ is not None:
    InitByValue.__doc__ = "Function version of operator initByValue whose __init__ function is \n" + initByValue.__init__.__doc__

def PyEval(pop, *args, **kwargs):
    pyEval(*args, **kwargs).apply(pop)

if pyEval.__init__.__doc__ is not None:
    PyEval.__doc__ = "Function version of operator pyEval whose __init__ function is \n" + pyEval.__init__.__doc__

def PyExec(pop, *args, **kwargs):
    pyExec(*args, **kwargs).apply(pop)

if pyExec.__init__.__doc__ is not None:
    PyExec.__doc__ = "Function version of operator pyExec whose __init__ function is \n" + pyExec.__init__.__doc__

def Stat(pop, *args, **kwargs):
    stat(*args, **kwargs).apply(pop)

if stat.__init__.__doc__ is not None:
    Stat.__doc__ = "Function version of operator stat whose __init__ function is \n" + stat.__init__.__doc__

def KamMutate(pop, *args, **kwargs):
    kamMutator(*args, **kwargs).apply(pop)

if kamMutator.__init__.__doc__ is not None:
    KamMutate.__doc__ = "Function version of operator kamMutator whose __init__ function is \n" + kamMutator.__init__.__doc__

def SmmMutate(pop, *args, **kwargs):
    smmMutator(*args, **kwargs).apply(pop)

if smmMutator.__init__.__doc__ is not None:
    SmmMutate.__doc__ = "Function version of operator smmMutator whose __init__ function is \n" + smmMutator.__init__.__doc__

def GsmMutate(pop, *args, **kwargs):
    gsmMutator(*args, **kwargs).apply(pop)

if gsmMutator.__init__.__doc__ is not None:
    GsmMutate.__doc__ = "Function version of operator gsmMutator whose __init__ function is \n" + gsmMutator.__init__.__doc__

def PyMutate(pop, *args, **kwargs):
    pyMutator(*args, **kwargs).apply(pop)

if pyMutator.__init__.__doc__ is not None:
    PyMutate.__doc__ = "Function version of operator pyMutator whose __init__ function is \n" + pyMutator.__init__.__doc__

def PointMutate(pop, *args, **kwargs):
    pointMutator(*args, **kwargs).apply(pop)

if pointMutator.__init__.__doc__ is not None:
    PointMutate.__doc__ = "Function version of operator pointMutator whose __init__ function is \n" + pointMutator.__init__.__doc__

def Migrate(pop, *args, **kwargs):
    migrator(*args, **kwargs).apply(pop)

if migrator.__init__.__doc__ is not None:
    Migrate.__doc__ = "Function version of operator migrator whose __init__ function is \n" + migrator.__init__.__doc__

def PyMigrate(pop, *args, **kwargs):
    pyMigrator(*args, **kwargs).apply(pop)

if pyMigrator.__init__.__doc__ is not None:
    PyMigrate.__doc__ = "Function version of operator pyMigrate whose __init__ function is \n" + pyMigrator.__init__.__doc__

def SplitSubPop(pop, *args, **kwargs):
    splitSubPop(*args, **kwargs).apply(pop)

if splitSubPop.__init__.__doc__ is not None:
    SplitSubPop.__doc__ = "Function version of operator splitSubPops whose __init__ function is \n" + splitSubPop.__init__.__doc__

def MergeSubPops(pop, *args, **kwargs):
    mergeSubPops(*args, **kwargs).apply(pop)

if mergeSubPops.__init__.__doc__ is not None:
    MergeSubPops.__doc__ = "Function version of operator mergeSubPops whose __init__ function is \n" + mergeSubPops.__init__.__doc__

def ResizeSubPops(pop, *args, **kwargs):
    resizeSubPops(*args, **kwargs).apply(pop)

if resizeSubPops.__init__.__doc__ is not None:
    ResizeSubPops.__doc__ = "Function version of operator resizeSubPops whose __init__ function is \n" + resizeSubPops.__init__.__doc__

def RemoveSubPops(pop, *args, **kwargs):
    pop.removeSubPops(*args, **kwargs)

if population.removeSubPops.__doc__ is not None:
    RemoveSubPops.__doc__ = "Function versionof member function population::removeSubPop with help info:\n" + population.removeSubPops.__doc__

def MapSelect(pop, *args, **kwargs):
    mapSelector(stage=PostMating, *args, **kwargs).apply(pop)

if mapSelector.__init__.__doc__ is not None:
    MapSelect.__doc__ = "Function version of operator mapSelect whose __init__ function is \n" + mapSelector.__init__.__doc__

def MaSelect(pop, *args, **kwargs):
    maSelector(stage=PostMating, *args, **kwargs).apply(pop)

if maSelector.__init__.__doc__ is not None:
    MaSelect.__doc__ = "Function version of operator maSelect whose __init__ function is \n" + maSelector.__init__.__doc__

def MlSelect(pop, *args, **kwargs):
    mlSelector(stage=PostMating, *args, **kwargs).apply(pop)

if mlSelector.__init__.__doc__ is not None:
    MlSelect.__doc__ = "Function version of operator mlSelect whose __init__ function is \n" + mlSelector.__init__.__doc__

def PySelect(pop, *args, **kwargs):
    pySelector(stage=PostMating, *args, **kwargs).apply(pop)

if pySelector.__init__.__doc__ is not None:
    PySelect.__doc__ = "Function version of operator pySelect whose __init__ function is \n" + pySelector.__init__.__doc__

def MapPenetrance(pop, *args, **kwargs):
    mapPenetrance(stage=PostMating, *args, **kwargs).apply(pop)

if mapPenetrance.__init__.__doc__ is not None:
    MapPenetrance.__doc__ = "Function version of operator mapPenetrance whose __init__ function is \n" + mapPenetrance.__init__.__doc__

def MaPenetrance(pop, *args, **kwargs):
    maPenetrance(stage=PostMating, *args, **kwargs).apply(pop)

if maPenetrance.__init__.__doc__ is not None:
    MaPenetrance.__doc__ = "Function version of operator maPenetrance whose __init__ function is \n" + maPenetrance.__init__.__doc__

def MlPenetrance(pop, *args, **kwargs):
    mlPenetrance(stage=PostMating, *args, **kwargs).apply(pop)

if mlPenetrance.__init__.__doc__ is not None:
    MlPenetrance.__doc__ = "Function version of operator mlPenetrance whose __init__ function is \n" + mlPenetrance.__init__.__doc__

def PyPenetrance(pop, *args, **kwargs):
    pyPenetrance(stage=PostMating, *args, **kwargs).apply(pop)

if pyPenetrance.__init__.__doc__ is not None:
    PyPenetrance.__doc__ = "Function version of operator pyPenetrance whose __init__ function is \n" + pyPenetrance.__init__.__doc__

def MapQuanTrait(pop, *args, **kwargs):
    mapQuanTrait(*args, **kwargs).apply(pop)

if mapQuanTrait.__init__.__doc__ is not None:
    MapQuanTrait.__doc__ = "Function version of operator mapQuanTrait whose __init__ function is \n" + mapQuanTrait.__init__.__doc__

def MaQuanTrait(pop, *args, **kwargs):
    maQuanTrait(*args, **kwargs).apply(pop)

if maQuanTrait.__init__.__doc__ is not None:
    MaQuanTrait.__doc__ = "Function version of operator maQuanTrait whose __init__ function is \n" + maQuanTrait.__init__.__doc__

def MlQuanTrait(pop, *args, **kwargs):
    mlQuanTrait(*args, **kwargs).apply(pop)

if mlQuanTrait.__init__.__doc__ is not None:
    MlQuanTrait.__doc__ = "Function version of operator mlQuanTrait whose __init__ function is \n" + mlQuanTrait.__init__.__doc__

def PyQuanTrait(pop, *args, **kwargs):
    pyQuanTrait(*args, **kwargs).apply(pop)

if pyQuanTrait.__init__.__doc__ is not None:
    PyQuanTrait.__doc__ = "Function version of operator pyQuanTrait whose __init__ function is \n" + pyQuanTrait.__init__.__doc__

def InfoEval(pop, *args, **kwargs):
    infoEval(*args, **kwargs).apply(pop)

if infoEval.__init__.__doc__ is not None:
    InfoEval.__doc__ = "Function version of operator infoEval whose __init__function is \n" + infoEval.__init__.__doc__

def InfoExec(pop, *args, **kwargs):
    infoExec(*args, **kwargs).apply(pop)

if infoExec.__init__.__doc__ is not None:
    InfoExec.__doc__ = "Function version of operator infoExec whose __init__function is \n" + infoExec.__init__.__doc__


# mating schemes

def cloneOffspringGenerator(ops=[], *args, **kwargs):
    'An offspring generator that uses cloneGenoTransmitter()'
    return  offspringGenerator([cloneGenoTransmitter()] + ops, 0, *args, **kwargs)

def mendelianOffspringGenerator(ops=[], *args, **kwargs):
    'An offspring generator that uses mendelianGenoTransmitter()'
    return  offspringGenerator([mendelianGenoTransmitter()] + ops, 2, *args, **kwargs)

def haplodiploidOffspringGenerator(ops=[], *args, **kwargs):
    'An offspring generator that uses haplodiploidGenoTransmitter()'
    return  offspringGenerator([haplodiploidGenoTransmitter()] + ops, 2, *args, **kwargs)

def selfingOffspringGenerator(ops=[], *args, **kwargs):
    'An offspring generator that uses selfingGenoTransmitter()'
    return  offspringGenerator([selfingGenoTransmitter()] + ops, 1, *args, **kwargs)

def cloneMating(numOffspring = 1., sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0):
    '''
    Note that
    selection is not considered (fitness is ignored)
    sequentialParentMating is used. If offspring (virtual) subpopulation size
   is smaller than parental subpopulation size, not all parents will be cloned.
   If offspring (virtual) subpopulation size is larger, some parents will be
   cloned more than once.
    numOffspring interface is respected.
    during mating operators are applied.
    '''
    return pyMating(
        chooser = sequentialParentChooser(),
        generator = cloneOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def randomSelection(numOffspring = 1., sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0):
    '''a mating scheme that uses binomial selection, regardless of sex
   No sex information is involved (binomial random selection). Offspring is chosen from parental generation
   by random or according to the fitness values.
   In this mating scheme,
    numOffspring protocol is honored;
    population size changes are allowed;
    selection is possible;
    haploid population is allowed.
    '''
    return pyMating(
        chooser = randomParentChooser(),
        generator = cloneOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def randomMating(numOffspring = 1., sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0):
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
        chooser = randomParentsChooser(replacement=True),
        generator = mendelianOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def monogamousMating(numOffspring = 1., sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0):
    '''
   This mating scheme is identical to random mating except that parents
   are chosen without replacement. Under this mating scheme, offspring share
   the same mother must share the same father. In case that all parental
   pairs are exhausted, parameter \c replenish=True allows for the replenishment
   of one or both sex groups.
    '''
    return pyMating(
        chooser = randomParentsChooser(replacement=False),
        generator = mendelianOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def polygamousMating(polySex=Male, polyNum=1,
        numOffspring = 1., sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0):
    '''
   This mating scheme is composed of a random parents chooser that allows for
   polygamous mating, and a mendelian offspring generator. In this mating scheme,
   a male (or female) parent will have more than one sex partner (\c numPartner).
   Parents returned from this parents chooser will yield the same male (or female)
   parents, each with varying partners.
    '''
    return pyMating(
        chooser = polyParentsChooser(polySex, polyNum),
        generator = mendelianOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def alphaMating(alphaSex=Male, alphaNum=0, alphaField='',
        numOffspring = 1., 
		 sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0):
    '''
     Only a number of alpha individuals can mate with individuals of opposite sex.

   This mating scheme is composed of an random parents chooser with alpha individuals,
   and a Mendelian offspring generator. That is to say, a certain number of alpha
   individual (male or female) are determined by \c alphaNum or an information field. Then,
   only these alpha individuals are able to mate with random individuals of
   opposite sex.

   alphaSex the sex of the alpha individual, i.e. alpha male
	           or alpha female who be the only mating individuals in their
	           sex group.
   alphaNum Number of alpha individuals. If \c infoField is
	           not given, \c alphaNum random individuals with \c alphaSex
	           will be chosen. If selection is enabled, individuals with higher+  
               fitness values have higher probability to be selected. There is
	           by default no alpha individual (\c alphaNum = 0).
   alphaField if an information field is given, individuals
	           with non-zero values at this information field are alpha individuals.
	           Note that these individuals must have \c alphaSex.

   Please refer to class \c mating for descriptions of other parameters.
   Note: If selection is enabled, it works regularly on on-alpha sex, but
	           works twice on alpha sex. That is to say, \c alphaNum alpha indiviudals
	           are chosen selectively, and selected again during mating.

    '''
    return pyMating(
        chooser = alphaParentsChooser(alphaSex, alphaNum, alphaField),
        generator = mendelianOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def haplodiploidMating(replacement=True,
		numOffspring = 1., 
		 sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0):
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
        chooser = randomParentsChooser(replacement),
        generator = haplodiploidOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def selfMating(replacement=True, numOffspring = 1.,	 sexMode = RandomSex,
        ops = [], subPopSize = [],
		subPop = (), weight = 0):
    '''
    In this mating scheme, a parent is choosen randomly, acts
    both as father and mother in the usual random mating. The parent
    is chosen randomly, regardless of sex. If selection is turned on,
    the probability that an individual is chosen is proportional to
    his/her fitness.
    '''
    return pyMating(
        chooser = randomParentChooser(replacement),
        generator = selfingOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def consanguineousMating(infoFields = [], func = None, param = None,
        replacement = False, numOffspring = 1.,	 sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0):
    '''
   This mating scheme randomly choose a parent and then choose his/her spouse from indexes
   stored in \c infoFields.

    relativeFields The information fields that stores indexes to other individuals
	    in a population. If more than one valid (positive value) indexes exist, a random
	    index will be chosen. (c.f. \c infoParentsChooser ) If there is no individual
	    having any valid index, the second parent will be chosen randomly from the
	    whole population.

    func A python function that can be used to prepare the indexes of these
	    information fields. For example, functions population::locateRelatives and/or
	    population::setIndexesOfRelatives can be used to locate certain types of relatives
	    of each individual.

	 param An optional parameter that can be passed to \c func.

	   Please refer to \c infoParentsChooser and \c mendelianOffspringGenerator for
	   other parameters.
    '''
    return pyMating(
        chooser = infoParentsChooser(infoFields, func, param, replacement),
        generator = mendelianOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize, subPop = subPop,
        weight = weight)


def controlledRandomMating(loci=[], alleles=[], freqFunc=None,
        numOffspring = 1., sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0):
    '''
    This is the controlled random mating scheme described in
    *Peng 2007 (PLoS Genetics)*. Basically, a *freqFunc*
    is passed to this mating scheme and set the allele frequencies of given
    alleles at given loci at the offspring generation.
    \n
    The offspring generation is conceptually populated in two steps.
    At the first step, only families with disease alleles are accepted
    until the expected number of disease alleles are met. At the second
    step, only families with wide type alleles are accepted to populate
    the rest of the offspring generation.
    \n
    \param loci loci at which allele frequencies are monitored (controlled)
	\param alleles alleles at given loci. It should have the same length as \c loci
	\param freqFunc a Python function that accepts a generation number and returns
	  expected allele frequencies at given loci
	\param acceptScheme internal use only
    '''
    return pyMating(chooser = randomParentsChooser(True),
        generator = controlledOffspringGenerator(loci, alleles, freqFunc,
            [mendelianGenoTransmitter()] + ops, 2, numOffspring,
            sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)



# Utility functions

def MergePopulations(pops, subPopSizes=[], keepAncestralPops=-1):
    'merge several populations with the same genotypic structure and create a new population'
    if len(pops) == 0:
        raise exceptions.ValueError("MergePopuations: empty population list is given")
    elif len(pops) == 1:
        return pops[0].clone()
    # to avoid repeated merging, it is better to merge files two by two
    merged = []
    for i in range(len(pops)):
        if i*2 < len(pops):
            merged.append(pops[i*2].clone())
        if i*2 + 1 < len(pops):
            merged[i].mergePopulation(pops[i*2+1], keepAncestralPops=keepAncestralPops)
    # merge merged with minimal population copying
    while True:
        count = 0
        for i in range(len(merged)):
            if merged[i] is not None:
                count += 1
                for j in range(i+1, len(merged)):
                    if merged[j] is not None:
                        merged[i].mergePopulation(merged[j], keepAncestralPops=keepAncestralPops)
                        merged[j] = None
                        break
        if count == 1:
            break
    res = merged[0]
    if len(subPopSizes) != 0:
        if sum(subPopSizes) != res.popSize():
            raise exceptions.ValueError("MergePopulations: can not change total population size")
        res.setSubPopStru(subPopSizes)
    return res


def MergePopulationsByLoci(pops, newNumLoci=[], newLociPos=[], byChromosome=False):
    'merge several populations of the same size by loci and create a new population'
    if len(pops) == 0:
        raise exceptions.ValueError("MergePopuations: empty population list is given")
    elif len(pops) == 1:
        return pops[0].clone()
    # to avoid repeated merging, it is better to merge files two by two
    merged = []
    for i in range(len(pops)):
        if i*2 < len(pops):
            merged.append(pops[i*2].clone())
        if i*2 + 1 < len(pops):
            merged[i].mergePopulationByLoci(pops[i*2+1], [], [], byChromosome)
    # merge merged with minimal population copying
    while True:
        count = 0
        for i in range(len(merged)):
            if merged[i] is not None:
                count += 1
                for j in range(i+1, len(merged)):
                    if merged[j] is not None:
                        merged[i].mergePopulationByLoci(merged[j], [], [], byChromosome)
                        merged[j] = None
                        break
        if count == 1:
            break
    res = merged[0]
    if len(newNumLoci) != 0:
        if sum(newNumLoci) != res.totNumLoci():
            raise exceptions.ValueError("MergePopulationsByLoci: can not change total number of loci")
        if sum(newLociPos) != res.totNumLoci():
            raise exceptions.ValueError("MergePopulationsByLoci: can not change total number of loci")
        res.rearrangeLoci(newNumLoci, newLociPos)
    return res



