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
from exceptions import ImportError

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

if simuOptions['Revision'] is not None and simuRev() < simuOptions['Revision']:
    raise ImportError('simuPOP version %s (revision %d) is installed ' % (simuVer(), simuRev()) +
        'but simuPOP revision >= %d is required. ' % simuOptions['Revision'] +
        'Please consider upgrading your simuPOP installation.')

_dbgCode = {
	'DBG_ALL': DBG_ALL,
	'DBG_GENERAL': DBG_GENERAL,
	'DBG_UTILITY': DBG_UTILITY,
	'DBG_OPERATOR': DBG_OPERATOR,
	'DBG_SIMULATOR': DBG_SIMULATOR,
	'DBG_INDIVIDUAL': DBG_INDIVIDUAL,
	'DBG_OUTPUTER': DBG_OUTPUTER,
	'DBG_MUTATOR': DBG_MUTATOR,
	'DBG_TRANSMITTER': DBG_TRANSMITTER,
	'DBG_INITIALIZER': DBG_INITIALIZER,
	'DBG_POPULATION': DBG_POPULATION,
	'DBG_STATOR': DBG_STATOR,
	'DBG_TERMINATOR': DBG_TERMINATOR,
	'DBG_TAGGER': DBG_TAGGER,
	'DBG_VISUALIZER': DBG_VISUALIZER,
	'DBG_SELECTOR': DBG_SELECTOR,
	'DBG_MATING': DBG_MATING,
	'DBG_MIGRATOR': DBG_MIGRATOR,
	'DBG_PROFILE': DBG_PROFILE,
	'DBG_DEVEL': DBG_DEVEL,
}

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
    TurnOnDebug(_dbgCode['DBG_GENERAL'])

if simuOptions['Debug'] != []:
    for g in simuOptions['Debug']:
        if g not in ['', None]:
            try:
                print "Turn on debug '%s'" % g
                TurnOnDebug(_dbgCode[g])
            except:
                print 'Incorrect debug code. Please use one of'
                print _dbgCode.keys()

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


def DebugCodes():
    'Return names of all debug codes'
    return _dbgCode.keys()

#
# functions to corresponding operators
def Dump(pop, *args, **kwargs):
    'Apply operator ``dumper`` to population *pop*.'
    dumper(*args, **kwargs).apply(pop)

def InitSex(pop, *args, **kwargs):
    'Apply operator ``initSex`` to population *pop*.'
    initSex(*args, **kwargs).apply(pop)

def InitByFreq(pop, *args, **kwargs):
    'Apply operator ``initByFreq`` to population *pop*.'
    initByFreq(*args, **kwargs).apply(pop)

def InitByValue(pop, *args, **kwargs):
    'Apply operator ``initByValue`` to population *pop*.'
    initByValue(*args, **kwargs).apply(pop)

def PyEval(pop, *args, **kwargs):
    '''Evaluate statements *stmts* (optional) and expression *expr* in
    population *pop*\ 's local namespace and return the result of *expr*.
    If *exposePop* is given, population *pop* will be exposed in its local
    namespace as a variable with a name specified by *exposePop*.

    .. note::

       Unlike its operator counterpart, this function returns the result of
       *expr* rather than writting it to an output.
    '''
    return pyEval(*args, **kwargs).evaluate(pop)

def PyExec(pop, *args, **kwargs):
    '''Execute *stmts* in population *pop*\ 's local namespace.'''
    pyExec(*args, **kwargs).apply(pop)

def InfoEval(pop, *args, **kwargs):
    '''Evaluate *expr* for each individual, using information fields as variables.
    Please refer to operator ``infoEval`` for details.
    '''
    infoEval(*args, **kwargs).apply(pop)

def InfoExec(pop, *args, **kwargs):
    '''Execute *stmts* for each individual, using information fields as variables.
    Please refer to operator ``infoExec`` for details.
    '''
    infoExec(*args, **kwargs).apply(pop)

def Migrate(pop, *args, **kwargs):
    'Function form of operator ``migrator``.'
    migrator(*args, **kwargs).apply(pop)

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

def MapSelect(pop, loci, fitness, phase = False, *args, **kwargs):
    mapSelector(loci, fitness, phase, PostMating, *args, **kwargs).apply(pop)

if mapSelector.__init__.__doc__ is not None:
    MapSelect.__doc__ = "Function version of operator mapSelect whose __init__ function is \n" + mapSelector.__init__.__doc__

def MaSelect(pop, loci, fitness, wildtype, *args, **kwargs):
    maSelector(loci, fitness, wildtype, PostMating, *args, **kwargs).apply(pop)

if maSelector.__init__.__doc__ is not None:
    MaSelect.__doc__ = "Function version of operator maSelect whose __init__ function is \n" + maSelector.__init__.__doc__

def MlSelect(pop, selectors, mode = Multiplicative, *args, **kwargs):
    mlSelector(selectors, mode, PostMating, *args, **kwargs).apply(pop)

if mlSelector.__init__.__doc__ is not None:
    MlSelect.__doc__ = "Function version of operator mlSelect whose __init__ function is \n" + mlSelector.__init__.__doc__

def PySelect(pop, loci, func, *args, **kwargs):
    pySelector(loci, func, PostMating, *args, **kwargs).apply(pop)

if pySelector.__init__.__doc__ is not None:
    PySelect.__doc__ = "Function version of operator pySelect whose __init__ function is \n" + pySelector.__init__.__doc__

def MapPenetrance(pop, loci, penetrance, phase = False, ancGen = -1, *args, **kwargs):
    mapPenetrance(loci, penetrance, phase, ancGen, PostMating, *args, **kwargs).apply(pop)

if mapPenetrance.__init__.__doc__ is not None:
    MapPenetrance.__doc__ = "Function version of operator mapPenetrance whose __init__ function is \n" + mapPenetrance.__init__.__doc__

def MaPenetrance(pop, loci, penetrance, wildtype = 0, ancGen = -1, *args, **kwargs):
    maPenetrance(loci, penetrance, wildtype, ancGen, PostMating, *args, **kwargs).apply(pop)

if maPenetrance.__init__.__doc__ is not None:
    MaPenetrance.__doc__ = "Function version of operator maPenetrance whose __init__ function is \n" + maPenetrance.__init__.__doc__

def MlPenetrance(pop, peneOps, mode = Multiplicative, ancGen = -1, *args, **kwargs):
    mlPenetrance(peneOps, mode, ancGen, PostMating, *args, **kwargs).apply(pop)

if mlPenetrance.__init__.__doc__ is not None:
    MlPenetrance.__doc__ = "Function version of operator mlPenetrance whose __init__ function is \n" + mlPenetrance.__init__.__doc__

def PyPenetrance(pop, loci, func, ancGen = -1, *args, **kwargs):
    pyPenetrance(loci, func, ancGen, PostMating, *args, **kwargs).apply(pop)

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


# mating schemes

def cloneOffspringGenerator(ops=[], *args, **kwargs):
    '''An offspring generator that uses ``cloneGenoTransmitter()`` as a
    default genotype transmitter. Additional during mating operators can be
    specified using the *ops* parameter. Other parameters are passed directly
    to ``offspringGenerator``.
    '''
    return offspringGenerator([cloneGenoTransmitter()] + ops, *args, **kwargs)

def mendelianOffspringGenerator(ops=[], *args, **kwargs):
    '''An offspring generator that uses ``mendelianGenoTransmitter()`` as a
    default genotype transmitter. Additional during mating operators can be
    specified using the *ops* parameter. Other parameters are passed directly
    to ``offspringGenerator``.
    '''
    return offspringGenerator([mendelianGenoTransmitter()] + ops, *args, **kwargs)

def haplodiploidOffspringGenerator(ops=[], *args, **kwargs):
    '''An offspring generator that uses ``haplodiploidGenoTransmitter()`` as a
    default genotype transmitter. Additional during mating operators can be
    specified using the *ops* parameter. Other parameters are passed directly
    to ``offspringGenerator``.
    '''
    return offspringGenerator([haplodiploidGenoTransmitter()] + ops, *args, **kwargs)

def selfingOffspringGenerator(ops=[], *args, **kwargs):
    '''An offspring generator that uses ``selfingGenoTransmitter()`` as a
    default genotype transmitter. Additional during mating operators can be
    specified using the *ops* parameter. Other parameters are passed directly
    to ``offspringGenerator``.
    '''
    return offspringGenerator([selfingGenoTransmitter()] + ops, *args, **kwargs)


def cloneMating(numOffspring = 1, sexMode = None, ops = [], subPopSize = [],
        subPop = (), weight = 0, selectionField = None):
    '''A homogeneous mating scheme that uses a sequential parent chooser and
    a clone offspring generator. Please refer to class ``offspringGenerator``
    for parameters *ops* and *numOffspring*, and to class ``homoMating`` for
    parameters  *subPopSize*, *subPop* and *weight*. Parameters *sexMode* and
    *selectionField* are ignored because this mating scheme does not support
    natural selection, and ``cloneOffspringGenerator`` copies sex from parents
    to offspring.
    '''
    return homoMating(
        chooser = sequentialParentChooser(),
        generator = cloneOffspringGenerator(ops, numOffspring, RandomSex),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def randomSelection(numOffspring = 1, sexMode = None, ops = [], subPopSize = [],
        subPop = (), weight = 0, selectionField = 'fitness'):
    '''A homogeneous mating scheme that uses a random single-parent parent
    chooser with replacement, and a clone offspring generator. This mating
    scheme is usually used to simulate the basic haploid Wright-Fisher model
    but it can also be applied to diploid populations. Please refer to class
    ``randomParentChooser`` for parameter *selectionField*, to class
    ``offspringGenerator`` for parameters *ops* and *numOffspring*, and to
    class ``homoMating`` for parameters *subPopSize*, *subPop* and *weight*.
    Parameter *sexMode* is ignored because ``cloneOffspringGenerator`` copies
    sex from parents to offspring.
    '''
    return homoMating(
        chooser = randomParentChooser(True, selectionField),
        generator = cloneOffspringGenerator(ops, numOffspring, RandomSex),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def randomMating(numOffspring = 1, sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0, selectionField = 'fitness'):
    '''A homogeneous mating scheme that uses a random parents chooser with
    replacement and a Mendelian offspring generator. This mating scheme is
    widely used to simulate diploid sexual Wright-Fisher random mating.
    Please refer to class ``randomParentsChooser`` for parameter
    *selectionField*, to class ``offspringGenerator`` for parameters *ops*,
    *sexMode* and *numOffspring*, and to class ``homoMating`` for parameters
    *subPopSize*, *subPop* and *weight*.
    '''
    return homoMating(
        chooser = randomParentsChooser(True, selectionField),
        generator = mendelianOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def monogamousMating(numOffspring = 1, sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0, selectionField = None):
    '''A homogeneous mating scheme that uses a random parents chooser without
    replacement and a Mendelian offspring generator. It differs from the basic
    random mating scheme in that each parent can mate only once so there is no
    half-sibling in the population. Please refer to class ``offspringGenerator``
    for parameters *ops*, *sexMode* and *numOffspring*, and to class
    ``homoMating`` for parameters *subPopSize*, *subPop* and *weight*.
    Parameter *selectionField* is ignored because this mating scheme does not
    support natural selection.
    '''
    return homoMating(
        chooser = randomParentsChooser(replacement=False),
        generator = mendelianOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def polygamousMating(polySex=Male, polyNum=1, numOffspring = 1,
        sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0, selectionField = 'fitness'):
    '''A homogeneous mating scheme that uses a multi-spouse parents chooser
    and a Mendelian offspring generator. It differs from the basic random
    mating scheme in that each parent of sex *polySex* will have *polyNum*
    spouses. Please refer to class ``polyParentsChooser`` for parameters
    *polySex*, *polyNum* and *selectionField*, to class ``offspringGenerator``
    for parameters *ops*,  *sexMode* and *numOffspring*, and to class
    ``homoMating`` for parameters *subPopSize*, *subPop* and *weight*.
    '''
    return homoMating(
        chooser = polyParentsChooser(polySex, polyNum),
        generator = mendelianOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def alphaMating(alphaSex=Male, alphaNum=0, alphaField='', numOffspring = 1, 
    	sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0, selectionField = 'fitness'):
    '''A homogeneous mating scheme that uses a alpha-individual parents chooser
    and a Mendelian offspring generator. It differs from the basic random
    mating scheme in that selection of parents of sex *alphaSex* is limited to
    certain alpha individuals, which are chosen either randomly (parameter
    *alphaNum*) or from an information field (parameter *alphaField*). This
    mating scheme is usually used to simulate animal population where only a
    few alpha individuals have the right to mate. Please refer to class
    ``alphaParentsChooser`` for parameters *alphaSex*, *alphaNum*, *alphaField*
    and *selectionField*, to class ``offspringGenerator`` for parameters *ops*,
    *sexMode* and *numOffspring*, and to class ``homoMating`` for parameters
    *subPopSize*, *subPop* and *weight*.
    '''
    return homoMating(
        chooser = alphaParentsChooser(alphaSex, alphaNum, alphaField, selectionField),
        generator = mendelianOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def haplodiploidMating(numOffspring = 1., sexMode = RandomSex, ops = [],
        subPopSize = [], subPop = (), weight = 0, selectionField = 'fitness'):
    '''A homogeneous mating scheme that uses a random parents chooser with
    replacement and a haplodiploid offspring generator. It should be used
    in a haplodiploid population where male individuals only have one set
    of homologous chromosomes. Please refer to class ``randomParentsChooser``
    for parameter *selectionField*, to class ``offspringGenerator`` for
    parameters *ops*, *sexMode* and *numOffspring*, and to class ``homoMating``
    for parameters *subPopSize*, *subPop* and *weight*.
    '''
    return homoMating(
        chooser = randomParentsChooser(True, selectionField),
        generator = haplodiploidOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def selfMating(replacement=True, numOffspring = 1, sexMode = RandomSex,
        ops = [], subPopSize = [], subPop = (), weight = 0,
        selectionField = 'fitness'):
    '''A homogeneous mating scheme that uses a random single-parent parent
    chooser with or without replacement (parameter *replacement*) and a
    selfing offspring generator. It is used to mimic self-fertilization
    in certain plant populations. Please refer to class ``randomParentChooser``
    for parameter *replacement* and  *selectionField*, to class
    ``offspringGenerator`` for parameters *ops*, *sexMode* and *numOffspring*,
    and to class ``homoMating`` for parameters *subPopSize*, *subPop* and
    *weight*.
    '''
    return homoMating(
        chooser = randomParentChooser(replacement, selectionField),
        generator = selfingOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def consanguineousMating(infoFields = [], func = None, param = None,
        replacement = False, numOffspring = 1.,	 sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0, selectionField = 'fitness'):
    '''A homogeneous mating scheme that uses an information parents chooser and
    a Mendelian offspring generator. A function *func* should be defined to
    locate certain types of relative to each individual and save their indexes
    to information fields *infoFields*. This mating scheme will then choose a
    parent randomly and then another parent from his/her relatives using their
    saved indexes. Please refer to class ``infoParentsChooser`` for parameters
    *infoFields*, *func*, *param* and  *selectionField*, to class
    ``offspringGenerator`` for parameters *ops*, *sexMode* and *numOffspring*,
    and to class ``homoMating`` for parameters *subPopSize*, *subPop* and
    *weight*.
    '''
    return homoMating(
        chooser = infoParentsChooser(infoFields, func, param, selectionField),
        generator = mendelianOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


def controlledRandomMating(loci=[], alleles=[], freqFunc=None,
        numOffspring = 1, sexMode = RandomSex, ops = [], subPopSize = [],
		subPop = (), weight = 0, selectionField = 'fitness'):
    '''A homogeneous mating scheme that uses a random sexual parents chooser
    with replacement and a controlled offspring generator using Mendelian
    genotype transmitter. At each generation, function *freqFunc* will be
    called to obtain intended frequencies of alleles *alleles* at loci
    *loci*. The controlled offspring generator will control the acceptance of
    offspring so that the generation reaches desired allele frequencies at
    these loci. Rationals and applications of this mating scheme is described
    in details in a paper *Peng et al, 2007 (PLoS Genetics)*. Please refer to
    class ``randomParentsChooser`` for parameters *selectionField*, to class
    ``controlledOffspringGenerator`` for parameters *loci*, *alleles*,
    *freqFunc*, to class ``offspringGenerator`` for parameters *ops*, *sexMode*
    and *numOffspring*, and to class ``homoMating`` for parameters *subPopSize*,
    *subPop* and *weight*.
    '''
    return homoMating(chooser = randomParentsChooser(True, selectionField),
        generator = controlledOffspringGenerator(loci, alleles, freqFunc,
            [mendelianGenoTransmitter()] + ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)


# Ascertainment operators and functions

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
    '''This operator draws random individuals from a population repeatedly and
    forms a number of random samples. These samples can be put in the
    population's local namespace, or save to disk files. The function form
    of this operator returns a list of samples directly.
    '''
    def __init__(self, size, *args, **kwargs):
        '''Draw *size* random samples from a population *times* times. *size* can
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
        self.pedigree = pedigree(pop, fatherField='', motherField='')
        self.pedigree.addInfoField('sample', -1)
        return True

    def drawSample(self, pop):
        import random
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
 
RandomSample.__doc__ = "Function version of operator ``randomSample``."


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
        self.pedigree = pedigree(pop, fatherField='', motherField='')
        self.pedigree.addInfoField('sample', -1)
        self.pedigree.setVirtualSplitter(affectionSplitter())
        return True

    def drawSample(self, pop):
        import random
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
        self.pedigree.locateRelatives(Offspring, ['offspring0', 'offspring1'])
        self.pedigree.locateRelatives(Spouse, ['spouse'])
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
        import random
        pedindex = self.pedigree.infoIdx('pedindex')
        sample = self.pedigree.infoIdx('sample')
        #
        # clear information sample in case this operator is applied twice
        self.pedigree.setIndInfo([0], sample)
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
                else:
                    ind.setInfo(-1, sample)
        sample = pop.extract(field='sample', ancGen=1, ped=self.pedigree)
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





