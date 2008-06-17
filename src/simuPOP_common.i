///////////////////////////////////////////////////////////////////
//     Copyright (C) 2004 by Bo Peng
//     bpeng@rice.edu
//
//     $LastChangedDate$
//     $Rev$
//
//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with this program; if not, write to the
//     Free Software Foundation, Inc.,
//     59 Temple Place - Suite 330, Boston, MA    02111-1307, USA.
//
///////////////////////////////////////////////////////////////////

// for uint16_t
%include "stdint.i"

////////////////////////// INCLUDE FILES //////////////////////////
%{

#include "../config.h"
#include "simuPOP_cfg.h"
#include "genoStru.h"
#include "individual.h"
#include "population.h"
#include "virtualSubPop.h"
#include "operator.h"
#include "simulator.h"

#include "utility.h"
#include "misc.h"
#include "pedigree.h"
#include "initializer.h"
#include "outputer.h"
#include "terminator.h"
#include "mating.h"
#include "tagger.h"
#include "stator.h"
#include "migrator.h"
#include "mutator.h"
#include "recombinator.h"
#include "selector.h"
#include "qtrait.h"
#include "penetrance.h"
#include "sampler.h"

%}

////////////////////////// DEFINE CARRAY //////////////////////////
%{
extern "C"
{
#include "arraymodule.c"
}
%}

////////////////////////// CLEAN EXTRA SYMBOLS //////////////////////////

// do not load these constants in ../config.h
%ignore HAVE_DECL_ACOSH;
%ignore HAVE_DECL_ATANH;
%ignore HAVE_DECL_EXPM1;
%ignore HAVE_DECL_FINITE;
%ignore HAVE_DECL_FREXP;
%ignore HAVE_DECL_HYPOT;
%ignore HAVE_DECL_ISFINITE;
%ignore HAVE_DECL_ISINF;
%ignore HAVE_DECL_ISNAN;
%ignore HAVE_DECL_LDEXP;
%ignore HAVE_DEV_NULL;
%ignore HAVE_DLFCN_H;
%ignore HAVE_FLOAT_H;
%ignore HAVE_INTTYPES_H;
%ignore HAVE_ISWPRINT;
%ignore HAVE_LIMITS_H;
%ignore HAVE_MALLOC;
%ignore HAVE_MEMORY_H;
%ignore HAVE_MEMSET;
%ignore HAVE_PTRDIFF_T;
%ignore HAVE_SNPRINTF;
%ignore HAVE_STDBOOL_H;
%ignore HAVE_STDDEF_H;
%ignore HAVE_STDINT_H;
%ignore HAVE_STDLIB_H;
%ignore HAVE_STRINGS_H;
%ignore HAVE_STRING_H;
%ignore HAVE_STRTOUL;
%ignore HAVE_SYS_STAT_H;
%ignore HAVE_SYS_TYPES_H;
%ignore HAVE_UNISTD_H;
%ignore HAVE_BOOL;
%ignore PACKAGE_NAME;
%ignore PACKAGE_STRING;
%ignore PACKAGE_TARNAME;
%ignore PACKAGE_VERSION;
%ignore PACKAGE_BUGREPORT;
%ignore STDC_HEADERS;
%ignore GSL_RANGE_CHECK;


%ignore DBG_CODE_LENGTH;
%ignore SIMUPOP_VAR_NAME;
%ignore SIMUPOP_VER;
%ignore PopSWIGType;

////////////////////////// VECTOR/MATRIX TYPES //////////////////////////

%include "std_vector.i"
%include "std_string.i"

// stl.i for std::pair
%include "stl.i"
%include "std_map.i"

%include "../config.h"
%include "simuPOP_cfg.h"

namespace std
{
    // used in invidiaul.h
    %template()            pair<UINT, UINT>;
    // used in populaiton.h
    %template()            pair<UINT, ULONG>;
    %template()            vector<Allele>;

//#ifndef LONGALLELE
    %template()         vector<UINT>;
//#endif

    %template()         vector<int>;
    %template()         vector<LONG>;
    %template()         vector<ULONG>;
    %template()         vector<InfoType>;
    %template()         vector<double>;
    %template()         vector<string>;

    %template()         pair<string, double>;
    %template()         map<string, double>;
    %template()         map<int, int>;

    %template()         pair<ULONG, ULONG>;
    %template()         vector<pair<ULONG, ULONG> >;

    %template()         vector< vector<int> >;
    %template()         vector< vector<double> >;
}


////////////////////////// SWIG_INIT FUNCTION //////////////////////////
%init
%{
    simuPOP::initialize();
%}


////////////////////////// C++=>PYTHON EXCEPTIONS //////////////////////////

%include exception.i

%exception
{
    try
    {
        $function
    }
    catch(simuPOP::StopIteration e)
    {
        SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
        SWIG_fail;
    }
    catch(simuPOP::OutOfMemory e)
    {
        SWIG_exception(SWIG_MemoryError, e.message());
    }
    catch(simuPOP::IOError e)
    {
        SWIG_exception(SWIG_IOError, e.message());
    }
    catch(simuPOP::IndexError e)
    {
        SWIG_exception(SWIG_IndexError, e.message());
    }
    catch(simuPOP::TypeError e)
    {
        SWIG_exception(SWIG_TypeError, e.message());
    }
    catch(simuPOP::ValueError e)
    {
        SWIG_exception(SWIG_ValueError, e.message());
    }
    catch(simuPOP::SystemError e)
    {
        SWIG_exception(SWIG_SystemError, e.message());
    }
    catch(...)
    {
        SWIG_exception(SWIG_RuntimeError, "Unknown runtime error happened.");
    }
}


////////////////////////// SIMUPOP TYPES //////////////////////////

%ignore std::operator<<(ostream&, const strDict&);
%ignore std::operator<<(ostream&, const intDict&);
%ignore simuPOP::IndAlleleIterator;
%ignore simuPOP::IndInfoIterator;

// individual and population are type names, and can not be used
// as function name. ind and pop are used instead.
// at the python level, individual and population are better.
%rename(individual) ind(ULONG, UINT);
%rename(population) pop(UINT) const;

%newobject LoadPopulation;
%newobject LoadSimulator;

%newobject simuPOP::population::newPopByIndID;
%newobject simuPOP::population::newPopWithPartialLoci;
%newobject simuPOP::population::clone;
%newobject simuPOP::simulator::getPopulation;
%newobject simuPOP::simulator::clone;
%newobject simuPOP::baseOperator::clone;
%newobject simuPOP::mating::clone;
%newobject simuPOP::stat::clone;

// the following load a docstring file extracted from doxgen output.
// there will also be a bunch of %ignore directives as well
//
%include "simuPOP_doc.i"

%include "utility.h"
%include "misc.h"
%include "pedigree.h"
%include "genoStru.h"
%include "individual.h"

namespace std
{
    %template()    vector<simuPOP::vsp >;
    %template()    vector<simuPOP::vspSplitter * >;
}

%include "virtualSubPop.h"
%include "population.h"

%include "operator.h"

%extend simuPOP::population
{
    %template(setIndInfo) setIndInfo<vectori, UINT>;
    %template(setIndInfo) setIndInfo<vectorf, UINT>;
    %template(setIndInfo) setIndInfo<vectori>;
    %template(setIndInfo) setIndInfo<vectorf>;
}


namespace std
{
    %template()    vector<PyObject*>;
    %template()    vector<simuPOP::baseOperator * >;
    %template()    vector<simuPOP::mating * >;
}

////////////////////////// SIMUPOP CLASSES //////////////////////////

%include "mating.h"
%include "simulator.h"
%include "stator.h"
%include "outputer.h"
%include "initializer.h"
%include "terminator.h"
%include "tagger.h"
%include "migrator.h"
%include "mutator.h"
%include "recombinator.h"
%include "selector.h"
%include "qtrait.h"
%include "penetrance.h"
%include "sampler.h"

////////////////////////// Provide a iterator as pop.individuals() ////////////////
// This part is now implemented at the C++ level (see population.h)
// Performance comparison shows 12.18s vs. 4.10s, the individual(i) method
// takes 4.7 to 5.1s
//
//%pythoncode %{
//
//class populationIterator:
//    def __init__(self, pop):
//        self.pop = pop
//        self.index = -1
//        self.endIndex = pop.popSize()
//
//    def __iter__(self):
//        return self
//
//    def next(self):
//        self.index += 1
//        if self.index == self.endIndex:
//            raise StopIteration
//        else:
//            return self.pop.individual(self.index)
//
//def inds(self):
//    '''
//    Return an iterator iterate over all individuals
//    '''
//    return populationIterator(self)
//
//population.inds = inds
//%}


////////////////////////// SIMUPOP PYTHON UTILITY FUNCTIONS //////////////////////////

%pythoncode %{
import exceptions, types
# population.__init__ needs to check simuOptions
from simuOpt import simuOptions

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

def MergePopulations(pops, newSubPopSizes=[], keepAncestralPops=-1):
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
    if len(newSubPopSizes) != 0:
        if sum(newSubPopSizes) != res.popSize():
            raise exceptions.ValueError("MergePopulations: can not change total population size")
        res.setSubPopStru(newSubPopSizes)
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


#### /////////////////// FUNCTION COUNTERPART OF OPERATORS ////////////////////////////


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

def PyInit(pop, *args, **kwargs):
    pyInit(*args, **kwargs).apply(pop)

if pyInit.__init__.__doc__ is not None:
	PyInit.__doc__ = "Function version of operator pyInit whose __init__ function is \n" + pyInit.__init__.__doc__

def Spread(pop,    *args, **kwargs):
    spread(*args, **kwargs).apply(pop)

if spread.__init__.__doc__ is not None:
	Spread.__doc__ = "Function version of operator spread whose __init__ function is \n" + spread.__init__.__doc__

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

def RemoveEmptySubPops(pop, *args, **kwargs):
    pop.removeEmptySubPops(*args, **kwargs)

if population.removeEmptySubPops.__doc__ is not None:
	RemoveEmptySubPops.__doc__ = "Function versionof member function population::removeEmptySubPops with help info:\n" + population.removeEmptySubPops.__doc__

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

def Sample(pop, *args, **kwargs):
    s = sample(*args, **kwargs)
    s.apply(pop)
    return s.sample(pop)

if sample.__init__.__doc__ is not None:
	Sample.__doc__ = "Function version of operator sample whose __init__function is \n" + sample.__init__.__doc__

def RandomSample(pop, *args, **kwargs):
    s = randomSample(*args, **kwargs)
    s.apply(pop)
    return s.samples(pop)

if randomSample.__init__.__doc__ is not None:
	RandomSample.__doc__ = "Function version of operator randomSample whose __init__function is \n" + randomSample.__init__.__doc__

def CaseControlSample(pop, *args, **kwargs):
    s = caseControlSample(*args, **kwargs)
    s.apply(pop)
    return s.samples(pop)

if caseControlSample.__init__.__doc__ is not None:
	CaseControlSample.__doc__ = "Function version of operator caseControlSample whose __init__function is \n" + caseControlSample.__init__.__doc__

def PySample(pop, *args, **kwargs):
    s = pySample(*args, **kwargs)
    s.apply(pop)
    return s.samples(pop)

if pySample.__init__.__doc__ is not None:
	PySample.__doc__ = "Function version of operator pySample whose __init__function is \n" + pySample.__init__.__doc__

def AffectedSibpairSample(pop, *args, **kwargs):
    s = affectedSibpairSample(*args, **kwargs)
    s.apply(pop)
    return s.samples(pop)

if affectedSibpairSample.__init__.__doc__ is not None:
	AffectedSibpairSample.__doc__ = "Function version of operator affectedSibpairSample whose __init__function is \n" + affectedSibpairSample.__init__.__doc__

def LargePedigreeSample(pop, *args, **kwargs):
    s = largePedigreeSample(*args, **kwargs)
    s.apply(pop)
    return s.samples(pop)

if largePedigreeSample.__init__.__doc__ is not None:
	LargePedigreeSample.__doc__ = "Function version of operator largePedigreeSample whose __init__function is \n" + largePedigreeSample.__init__.__doc__

def NuclearFamilySample(pop, *args, **kwargs):
    s = nuclearFamilySample(*args, **kwargs)
    s.apply(pop)
    return s.samples(pop)

if nuclearFamilySample.__init__.__doc__ is not None:
	NuclearFamilySample.__doc__ = "Function version of operator nuclearFamilySample whose __init__function is \n" + nuclearFamilySample.__init__.__doc__

def PySubset(pop, *args, **kwargs):
    s = pySubset(*args, **kwargs)
    s.apply(pop)

if pySubset.__init__.__doc__ is not None:
	PySubset.__doc__ = "Function version of operator pySubset whose __init__function is \n" + pySubset.__init__.__doc__

def InfoEval(pop, *args, **kwargs):
    infoEval(*args, **kwargs).apply(pop)

if infoEval.__init__.__doc__ is not None:
	InfoEval.__doc__ = "Function version of operator infoEval whose __init__function is \n" + infoEval.__init__.__doc__

def InfoExec(pop, *args, **kwargs):
    infoExec(*args, **kwargs).apply(pop)

if infoExec.__init__.__doc__ is not None:
	InfoExec.__doc__ = "Function version of operator infoExec whose __init__function is \n" + infoExec.__init__.__doc__

def SavePopulation(pop, *args, **kwargs):
    pop.savePopulation(*args, **kwargs)

if population.savePopulation.__doc__ is not None:
	SavePopulation.__doc__ = "Function versionof member function population::savePopulation with help info:\n" + population.savePopulation.__doc__

def SaveSimulator(simu, *args, **kwargs):
    simu.saveSimulator(*args, **kwargs)

if simulator.saveSimulator.__doc__ is not None:
	SaveSimulator.__doc__ = "Function versionof member function simulator::saveSimulator with help info:\n" + simulator.saveSimulator.__doc__

#### /////////////////// SIMUPOP PYTHON REDEFINITION FUNCTIONS ////////////////////////
def new_population(self, size=[], ploidy=2, loci=[], sexChrom=False,
    lociPos=[], subPop=[], ancestralDepth=0, chromNames=[], alleleNames=[], lociNames=[],
    maxAllele=ModuleMaxAllele, infoFields=[]):
    if subPop != []:
        print 'Parameter subPop is obsolete. Please use size instead'
        if size != 0:
            print 'In addition, you can not specify both size and subPop'
        size = subPop
    if type(size) in [type(0), type(0L)]:
        sp = [size]
    elif type(size) in [types.TupleType, types.ListType]:
        sp = size
    if type(infoFields) not in (type([]), type(())):
        raise exceptions.ValueError('infoFields needs to be an array')
    ld = lociPos
    if len(lociPos) > 0 and type(lociPos[0]) in [types.TupleType, types.ListType]:
        ld = []
        for i in range(0, len(lociPos)):
            if len(lociPos[i]) != loci[i]:
                raise exceptions.ValueError("Loci position specification for chromosome %d is wrong" % i)
            ld.extend( lociPos[i])
    ln = lociNames
    if len(lociNames) > 0 and type(lociNames[0]) in [types.TupleType, types.ListType]:
        ln = []
        for i in range(0, len(lociNames)):
            ln.extend(lociNames[i])
    cppModule.population_swiginit(self,
        cppModule.new_population(sp, ploidy, loci, sexChrom, ld,
            ancestralDepth, chromNames, alleleNames, ln, maxAllele, infoFields))

new_population.__doc__ = population.__init__.__doc__
del population.__init__
population.__init__ = new_population


def new_dumper(self, chrom=[], subPop=[], indRange=[], *args, **kwargs):
    # param chrom
    if type(chrom) == types.IntType:
        ch = [chrom]
    else:
        ch = chrom
    # parameter subPop
    if type(subPop) == types.IntType:
        sp = [subPop]
    else:
        sp = subPop
    # parameter indRange
    ir = indRange
    if len(indRange) > 0 and type(indRange[0]) in [types.TupleType, types.ListType]:
        ir = []
        for i in indRange:
            ir.extend(i)
    cppModule.dumper_swiginit(self,
        cppModule.new_dumper(chrom=ch, subPop=sp, indRange=ir, *args, **kwargs))

new_dumper.__doc__ = dumper.__init__.__doc__
del dumper.__init__
dumper.__init__ = new_dumper


def new_controlledMating(self, matingScheme, freqFunc, locus=-1, loci=[],
    allele=-1, alleles=[], *args, **kwargs):
    # parameter locus
    if locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        raise exceptions.TypeError('Please specify locus or loci')
    if allele != -1 and type(allele) in [types.IntType, types.LongType]:
        al = [allele]
    elif type(alleles) in [types.IntType, types.LongType]:
        al = [alleles]
    elif type(alleles) in [types.TupleType, types.ListType] and len(alleles)>0:
        al = alleles
    else:
        raise exceptions.TypeError('Please specify allele or alleles')
    cppModule.controlledMating_swiginit(self,
        cppModule.new_controlledMating(matingScheme=matingScheme,
            loci=loc, alleles=al, freqFunc=freqFunc, *args, **kwargs))

new_controlledMating.__doc__ = controlledMating.__init__.__doc__
del controlledMating.__init__
controlledMating.__init__ = new_controlledMating


def new_controlledRandomMating(self, freqFunc, locus=-1, loci=[],
    allele=-1, alleles=[], *args, **kwargs):
    # parameter locus
    if locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        raise exceptions.TypeError('Please specify locus or loci')
    if allele != -1 and type(allele) in [types.IntType, types.LongType]:
        al = [allele]
    elif type(alleles) in [types.IntType, types.LongType]:
        al = [alleles]
    elif type(alleles) in [types.TupleType, types.ListType] and len(alleles)>0:
        al = alleles
    else:
        raise exceptions.TypeError('Please specify allele or alleles')
    cppModule.controlledRandomMating_swiginit(self,
        cppModule.new_controlledRandomMating( loci=loc, alleles=al,
        freqFunc=freqFunc, *args, **kwargs))

new_controlledRandomMating.__doc__ = controlledRandomMating.__init__.__doc__
del controlledRandomMating.__init__
controlledRandomMating.__init__ = new_controlledRandomMating


def mutator_setRate(self, rate, loci=[], atLoci=[], *args, **kwargs):
    # rate -> [rate] if needed
    if type(rate)    in [types.IntType, types.FloatType]:
        r = [rate]
    else:
        r = rate
    #
    if atLoci != []:
        print 'Parameter atLoci is obsolete. Please use loci'
        loc = atLoci
    else:
        loc = loci
    return cppModule.mutator_setRate(self, r, loc, *args, **kwargs)

del mutator.setRate
mutator.setRate = mutator_setRate

def new_kamMutator(self, rate=[], loci=[], atLoci=[], *args, **kwargs):
    # parameter rate
    if type(rate) in [types.IntType, types.FloatType]:
        r = [rate]
    else:
        r = rate
    if atLoci != []:
        print 'Parameter atLoci is obsolete. Please use loci'
        loc = atLoci
    else:
        loc = loci
    cppModule.kamMutator_swiginit(self,
        cppModule.new_kamMutator(rate=r, loci=loc, *args, **kwargs))

new_kamMutator.__doc__ = kamMutator.__init__.__doc__
del kamMutator.__init__
kamMutator.__init__ = new_kamMutator



def new_smmMutator(self, rate=[], loci=[], atLoci=[], *args, **kwargs):
    # parameter rate
    if type(rate) in [types.IntType, types.FloatType]:
        r = [rate]
    else:
        r = rate
    if atLoci != []:
        print 'Parameter atLoci is obsolete. Please use loci'
        loc = atLoci
    else:
        loc = loci
    cppModule.smmMutator_swiginit(self,
        cppModule.new_smmMutator(rate=r, loci=loc, *args, **kwargs))

new_smmMutator.__doc__ = smmMutator.__init__.__doc__
del smmMutator.__init__
smmMutator.__init__ = new_smmMutator



def new_gsmMutator(self, rate=[], loci=[], atLoci=[], *args, **kwargs):
    # parameter rate
    if type(rate) in [types.IntType, types.FloatType]:
        r = [rate]
    else:
        r = rate
    if atLoci != []:
        print 'Parameter atLoci is obsolete. Please use loci'
        loc = atLoci
    else:
        loc = loci
    cppModule.gsmMutator_swiginit(self,
        cppModule.new_gsmMutator(rate=r, loci=loc, *args, **kwargs))

new_gsmMutator.__doc__ = gsmMutator.__init__.__doc__
del gsmMutator.__init__
gsmMutator.__init__ = new_gsmMutator



def new_pyMutator(self, rate=[], loci=[], atLoci=[], *args, **kwargs):
    # parameter rate
    if type(rate) in [types.IntType, types.FloatType]:
        r = [rate]
    else:
        r = rate
    if atLoci != []:
        print 'Parameter atLoci is obsolete. Please use loci'
        loc = atLoci
    else:
        loc = loci
    cppModule.pyMutator_swiginit(self,
        cppModule.new_pyMutator(rate=r, loci=loc, *args, **kwargs))

new_pyMutator.__doc__ = pyMutator.__init__.__doc__
del pyMutator.__init__
pyMutator.__init__ = new_pyMutator

def new_pointMutator(self, locus=None, loci=[], atLoci=[], *args, **kwargs):
    if atLoci != []:
        print 'Parameter atLoci is obsolete. Please use locus or loci'
        loc = atLoci
    elif locus != None:
        loc = [locus]
    else:
        loc = loci
    cppModule.pointMutator_swiginit(self,
        cppModule.new_pointMutator(loci=loc, *args, **kwargs))

new_pointMutator.__doc__ = pointMutator.__init__.__doc__
del pointMutator.__init__
pointMutator.__init__ = new_pointMutator


def new_migrator(self, rate, fromSubPop=[], toSubPop=[], *args, **kwargs):
    # parameter rate
    r = rate
    if type(rate) in [types.IntType, types.LongType, types.FloatType]:
        r = [[rate]]
    # if a single vector, [a,b] ==> [[a,b]]
    if type(rate) in [types.ListType, types.TupleType, types.FloatType]:
        if len(rate) == 0:
            raise exceptions.ValueError('Migration rate can not be empty')
        elif type(rate[0]) in [types.IntType, types.LongType, types.FloatType]:
            r = [rate]
    # parameter fromSubPop
    if type(fromSubPop) in [types.IntType, types.LongType]:
        fs = [fromSubPop]
    else:
        fs = [vsp(x) for x in fromSubPop]
    print 'from subpopulations ', fromSubPop, fs
    # parameter toSubPop
    if type(toSubPop) in [types.IntType, types.LongType]:
        ts = [toSubPop]
    else:
        ts = toSubPop
    cppModule.migrator_swiginit(self,
        cppModule.new_migrator(rate=r, fromSubPop=fs, toSubPop=ts, *args, **kwargs))

new_migrator.__doc__ = migrator.__init__.__doc__
del migrator.__init__
migrator.__init__ = new_migrator


def new_pyMigrator(self, fromSubPop=[], toSubPop=[], *args, **kwargs):
    # parameter fromSubPop
    if type(fromSubPop) in [types.IntType, types.LongType]:
        fs = [fromSubPop]
    else:
        fs = fromSubPop
    # parameter toSubPop
    if type(toSubPop) in [types.IntType, types.LongType]:
        ts = [toSubPop]
    else:
        ts = toSubPop
    cppModule.pyMigrator_swiginit(self,
        cppModule.new_pyMigrator(fromSubPop=fs, toSubPop=ts, *args, **kwargs))

new_pyMigrator.__doc__ = pyMigrator.__init__.__doc__
del pyMigrator.__init__
pyMigrator.__init__ = new_pyMigrator


def new_recombinator(self, intensity=-1, rate=[],
    locus=None, loci=[], afterLoci=[],
    maleIntensity=-1, maleRate=[],
    maleLocus=None, maleLoci=[], maleAfterLoci=[], *args, **kwargs):
    # parameter loci
    if len(afterLoci) != 0 and len(loci) == 0:
        print 'Warning: parameter afterLoci of recombinator is renamed to loci'
        loci = afterLoci
    # parameter maleLoci
    if len(maleAfterLoci) != 0 and len(maleLoci) == 0:
        print 'Warning: parameter maleAfterLoci of recombinator is renamed to maleLoci'
        maleLoci = maleAfterLoci
    # paramter locus
    if locus is not None and len(loci) == 0:
        loci = [locus]
    if maleLocus is not None and len(maleLoci) == 0:
        maleLoci = [maleLocus]
    # parameter rate
    if type(rate) in [types.IntType, types.FloatType]:
        if len(afterLoci) > 0:
            r = [rate]*len(afterLoci)
        else:
            r = [rate]
    else:
        r = rate
    # parameter maleRate
    if type(maleRate) in [types.IntType, types.FloatType]:
        if len(maleAfterLoci) > 0:
            mr = [maleRate]*len(maleAfterLoci)
        elif len(afterLoci) > 0:
            mr = [maleRate]*len(afterLoci)
        else:
            mr = [maleRate]
    else:
        mr = maleRate
    cppModule.recombinator_swiginit(self,
        cppModule.new_recombinator(intensity=intensity,
        rate=r, afterLoci=loci,
        maleIntensity=maleIntensity,
        maleRate=mr, maleAfterLoci=maleLoci,
        *args, **kwargs))

new_recombinator.__doc__ = recombinator.__init__.__doc__
del recombinator.__init__
recombinator.__init__ = new_recombinator



def new_initByFreq(self, alleleFreq=[], locus=-1, loci=[], atLoci=[], indRange=[], subPop=-1, subPops=[], *args, **kwargs):
    # parameter alleleFreq
    if len(alleleFreq) > 0 and type(alleleFreq[0]) in [types.IntType, types.LongType, types.FloatType]:
        af = [alleleFreq]
    else:
        af = alleleFreq
    # parameter loci
    if type(atLoci) in [types.IntType, types.LongType]:
        loc = [atLoci]
    elif type(atLoci) in [types.TupleType, types.ListType] and len(atLoci) > 0:
        loc = atLoci
    elif locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif locus != -1 and type(locus) in [types.ListType, types.TupleType] and len(locus) > 0:
        loc = locus
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        loc = []
    # parameter indRange
    if len(indRange) > 0 and type(indRange[0]) in    [types.IntType, types.LongType]:
        ir = [indRange]
    else:
        ir = indRange
    # parameter subpop
    if subPop != -1 and type(subPop) in [types.IntType, types.LongType]:
        sp = [subPop]
    elif subPop != -1 and type(subPop) in [types.ListType, types.TupleType] and len(subPop) > 0:
        sp = subPop
    elif type(subPops) in [types.IntType, types.LongType]:
        sp = [subPops]
    elif type(subPops) in [types.TupleType, types.ListType] and len(subPops)>0:
        sp = subPops
    else:
        sp = []
    cppModule.initByFreq_swiginit(self,
        cppModule.new_initByFreq(alleleFreq=af, loci=loc, indRange=ir, subPop=sp,
            *args, **kwargs))

new_initByFreq.__doc__ = initByFreq.__init__.__doc__
del initByFreq.__init__
initByFreq.__init__ = new_initByFreq


def new_initByValue(self, value=[], indRange=[], locus=-1, loci=[], atLoci=[], subPop=-1, subPops=[], *args, **kwargs):
    # parameter value
    if len(value) > 0 and type(value[0]) in [types.IntType, types.LongType]:
        val = [value]
    else:
        val = value
    # parameter indRange
    if len(indRange) > 0 and type(indRange[0]) in    [types.IntType, types.LongType]:
        ir = [indRange]
    else:
        ir = indRange
        # parameter loci
    if type(atLoci) in [types.IntType, types.LongType]:
        loc = [atLoci]
    elif type(atLoci) in [types.TupleType, types.ListType] and len(atLoci) > 0:
        loc = atLoci
    elif locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif locus != -1 and type(locus) in [types.ListType, types.TupleType] and len(locus) > 0:
        loc = locus
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        loc = []
    if subPop != -1 and type(subPop) in [types.IntType, types.LongType]:
        sp = [subPop]
    elif subPop != -1 and type(subPop) in [types.ListType, types.TupleType] and len(subPop) > 0:
        sp = subPop
    elif type(subPops) in [types.IntType, types.LongType]:
        sp = [subPops]
    elif type(subPops) in [types.TupleType, types.ListType] and len(subPops)>0:
        sp = subPops
    else:
        sp = []
    cppModule.initByValue_swiginit(self,
        cppModule.new_initByValue(value=val, loci=loc, indRange=ir, subPop=sp, *args, **kwargs))

new_initByValue.__doc__ = initByValue.__init__.__doc__
del initByValue.__init__
initByValue.__init__ = new_initByValue


def new_pyInit(self, locus=-1, loci=[], atLoci=[], indRange=[], *args, **kwargs):
    # parameter indRange
    if len(indRange) > 0 and type(indRange[0]) in    [types.IntType, types.LongType]:
        ir = [indRange]
    else:
        ir = indRange
    # parameter loci
    if type(atLoci) in [types.IntType, types.LongType]:
        loc = [atLoci]
    elif type(atLoci) in [types.TupleType, types.ListType] and len(atLoci) > 0:
        loc = atLoci
    elif locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif locus != -1 and type(locus) in [types.ListType, types.TupleType] and len(locus) > 0:
        loc = locus
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        loc = []
    cppModule.pyInit_swiginit(self,
        cppModule.new_pyInit(indRange=ir, loci=loc, *args, **kwargs))

new_pyInit.__doc__ = pyInit.__init__.__doc__
del pyInit.__init__
pyInit.__init__ = new_pyInit


def new_stat(self, haploFreq=[], LD=[], LD_param={}, association=[], association_param={},
    relGroups=[], relMethod=[], midValues=None, *args, **kwargs):
    # midValues is now obsolete
    if midValues is not None:
        print 'Parameter midValues is now obsolete. Please use the _param parameter of corresponding statistics'
    # parameter haploFreq
    if len(haploFreq) > 0 and type(haploFreq[0]) in [types.IntType, types.LongType]:
        hf = [haploFreq]
    else:
        hf = haploFreq
    # parameter LD
    if len(LD) > 0 and type(LD[0]) in    [types.IntType, types.LongType]:
        ld = [LD]
    else:
        ld = LD
    # parameters of LD, convert 'stat':['LD', 'LD_prime'] etc to 'LD':True, 'LD_prime':True
    ldp = {}
    for key in LD_param.keys():
        if 'stat' == key and type(LD_param['stat']) in [types.TupleType, types.ListType]:
            for stat in LD_param['stat']:
                ldp[stat] = True
        else:
            ldp[key] = LD_param[key]
    # parameters of association, convert 'stat':['ChiSq', 'UCU'] etc to 'ChiSq':True, 'UCU':True
    assp = {}
    for key in association_param.keys():
        if 'stat' == key and type(association_param['stat']) in [types.TupleType, types.ListType]:
            for stat in association_param['stat']:
                assp[stat] = True
        else:
            assp[key] = association_param[key]
    # parameter association
    if len(association) > 0 and type(association[0]) in [types.IntType, types.LongType]:
        Association = [association]
    else:
        Association = association
    # parameter relGroups
    if relGroups == []:
        rg = [[]]
        useSubPop = True
    elif len(relGroups) > 0 and type(relGroups[0]) in [types.IntType, types.LongType]:
        rg = [relGroups]
        useSubPop = True
    else:
        rg = relGroups
        useSubPop = False
    # parameter relMethod
    if type(relMethod) in [types.IntType, types.LongType]:
        rm = [relMethod]
    else:
        rm = relMethod
    cppModule.stat_swiginit(self,
        cppModule.new_stat(haploFreq=hf, LD=ld, LD_param=ldp,
            association=Association, association_param=assp,
            relGroups=rg, relBySubPop=useSubPop,
            relMethod =rm, *args, **kwargs))

new_stat.__doc__ = stat.__init__.__doc__
del stat.__init__
stat.__init__ = new_stat



def new_randomSample(self, size=[], *args, **kwargs):
    if type(size) in [types.IntType, types.LongType]:
        sz=[size]
    else:
        sz = size
    cppModule.randomSample_swiginit(self,
        cppModule.new_randomSample(size=sz, *args, **kwargs))

new_randomSample.__doc__ = randomSample.__init__.__doc__
del randomSample.__init__
randomSample.__init__ = new_randomSample


def new_caseControlSample(self, cases=[], controls=[], *args, **kwargs):
    if type(cases) in [types.IntType, types.LongType]:
        ca = [cases]
        spSample = False
    else:
        ca = cases
        spSample = True
    if type(controls) in [types.IntType, types.LongType]:
        ct = [controls]
        spSample = False
    else:
        ct = controls
        spSample = True
    cppModule.caseControlSample_swiginit(self,
        cppModule.new_caseControlSample(cases=ca, controls=ct,
            spSample=spSample, *args, **kwargs))

new_caseControlSample.__doc__ = caseControlSample.__init__.__doc__
del caseControlSample.__init__
caseControlSample.__init__ = new_caseControlSample


def new_affectedSibpairSample(self,size=[], *args, **kwargs):
    if type(size) in [types.IntType, types.LongType]:
        sz=[size]
    else:
        sz = size
    cppModule.affectedSibpairSample_swiginit(self,
        cppModule.new_affectedSibpairSample(size=sz, *args, **kwargs))

new_affectedSibpairSample.__doc__ = affectedSibpairSample.__init__.__doc__
del affectedSibpairSample.__init__
affectedSibpairSample.__init__ = new_affectedSibpairSample


def new_largePedigreeSample(self, size=[], *args, **kwargs):
    if type(size) in [types.IntType, types.LongType]:
        sz= [size]
    else:
        sz = size
    cppModule.largePedigreeSample_swiginit(self,
        cppModule.new_largePedigreeSample(size=sz, *args, **kwargs))

new_largePedigreeSample.__doc__ = largePedigreeSample.__init__.__doc__
del largePedigreeSample.__init__
largePedigreeSample.__init__ = new_largePedigreeSample


def new_nuclearFamilySample(self, size=[], *args, **kwargs):
    if type(size) in [types.IntType, types.LongType]:
        sz= [size]
    else:
        sz = size
    cppModule.nuclearFamilySample_swiginit(self,
        cppModule.new_nuclearFamilySample(size=sz, *args, **kwargs))

new_nuclearFamilySample.__doc__ = nuclearFamilySample.__init__.__doc__
del nuclearFamilySample.__init__
nuclearFamilySample.__init__ = new_nuclearFamilySample


def new_mapSelector(self, locus=-1, loci=[], subPop=-1, subPops=[], *args, **kwargs):
    if locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        raise exceptions.TypeError('Please specify locus or loci')
    if subPop != -1 and type(subPop) in [types.IntType, types.LongType]:
        sp = [subPop]
    elif subPop != -1 and type(subPop) in [types.ListType, types.TupleType] and len(subPop) > 0:
        sp = subPop
    elif type(subPops) in [types.IntType, types.LongType]:
        sp = [subPops]
    elif type(subPops) in [types.TupleType, types.ListType] and len(subPops)>0:
        sp = subPops
    else:
        sp = []
    cppModule.mapSelector_swiginit(self,
        cppModule.new_mapSelector(loci=loc, subPops=sp, *args, **kwargs))

new_mapSelector.__doc__ = mapSelector.__init__.__doc__
del mapSelector.__init__
mapSelector.__init__ = new_mapSelector

def new_maSelector(self, locus=-1, loci=[], wildtype=[0], subPop=-1, subPops=[], *args, **kwargs):
    if locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        raise exceptions.TypeError('Please specify locus or loci')
    if type(wildtype) in [types.IntType, types.LongType]:
        wt = [wildtype]
    else:
        wt = wildtype
    if subPop != -1 and type(subPop) in [types.IntType, types.LongType]:
        sp = [subPop]
    elif subPop != -1 and type(subPop) in [types.ListType, types.TupleType] and len(subPop) > 0:
        sp = subPop
    elif type(subPops) in [types.IntType, types.LongType]:
        sp = [subPops]
    elif type(subPops) in [types.TupleType, types.ListType] and len(subPops)>0:
        sp = subPops
    else:
        sp = []
    cppModule.maSelector_swiginit(self,
        cppModule.new_maSelector(loci=loc, wildtype=wt, subPops=sp, *args, **kwargs))

new_maSelector.__doc__ = maSelector.__init__.__doc__
del maSelector.__init__
maSelector.__init__ = new_maSelector

def new_mlSelector(self, selectors=[], subPop=-1, subPops=[], *args, **kwargs):
    if subPop != -1 and type(subPop) in [types.IntType, types.LongType]:
        sp = [subPop]
    elif subPop != -1 and type(subPop) in [types.ListType, types.TupleType] and len(subPop) > 0:
        sp = subPop
    elif type(subPops) in [types.IntType, types.LongType]:
        sp = [subPops]
    elif type(subPops) in [types.TupleType, types.ListType] and len(subPops)>0:
        sp = subPops
    else:
        sp = []
    cppModule.mlSelector_swiginit(self,
        cppModule.new_mlSelector(selectors, subPops=sp, *args, **kwargs))


new_mlSelector.__doc__ = mlSelector.__init__.__doc__
del mlSelector.__init__
mlSelector.__init__ = new_mlSelector

def new_pySelector(self, locus=-1, loci=[], subPop=-1, subPops=[], *args, **kwargs):
    if locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        raise exceptions.TypeError('Please specify locus or loci')
    if subPop != -1 and type(subPop) in [types.IntType, types.LongType]:
        sp = [subPop]
    elif subPop != -1 and type(subPop) in [types.ListType, types.TupleType] and len(subPop) > 0:
        sp = subPop
    elif type(subPops) in [types.IntType, types.LongType]:
        sp = [subPops]
    elif type(subPops) in [types.TupleType, types.ListType] and len(subPops)>0:
        sp = subPops
    else:
        sp = []
    cppModule.pySelector_swiginit(self,
        cppModule.new_pySelector(loci=loc, subPops=sp, *args, **kwargs))

new_pySelector.__doc__ = pySelector.__init__.__doc__
del pySelector.__init__
pySelector.__init__ = new_pySelector

# for backward compatibility, keep long penetrance parameter
def new_mapPenetrance(self, locus=-1, loci=[], penetrance={}, *args, **kwargs):
    if locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        raise exceptions.TypeError('Please specify locus or loci')
    cppModule.mapPenetrance_swiginit(self,
        cppModule.new_mapPenetrance(loci=loc, penet=penetrance, *args, **kwargs))

new_mapPenetrance.__doc__ = mapPenetrance.__init__.__doc__
del mapPenetrance.__init__
mapPenetrance.__init__ = new_mapPenetrance

def new_maPenetrance(self, locus=-1, loci=[], wildtype=[0], penetrance=[], *args, **kwargs):
    if locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        raise exceptions.TypeError('Please specify locus or loci')
    if type(wildtype) in [types.IntType, types.LongType]:
        wt = [wildtype]
    else:
        wt = wildtype
    cppModule.maPenetrance_swiginit(self,
        cppModule.new_maPenetrance(loci=loc, wildtype=wt, penet=penetrance, *args, **kwargs))

new_maPenetrance.__doc__ = maPenetrance.__init__.__doc__
del maPenetrance.__init__
maPenetrance.__init__ = new_maPenetrance

def new_pyPenetrance(self, locus=-1, loci=[], *args, **kwargs):
    if locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        raise exceptions.TypeError('Please specify locus or loci')
    cppModule.pyPenetrance_swiginit(self,
        cppModule.new_pyPenetrance(loci=loc, *args, **kwargs))

new_pyPenetrance.__doc__ = pyPenetrance.__init__.__doc__
del pyPenetrance.__init__
pyPenetrance.__init__ = new_pyPenetrance

def new_mapQuanTrait(self, locus=-1, loci=[], *args, **kwargs):
    if locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        raise exceptions.TypeError('Please specify locus or loci')
    cppModule.mapQuanTrait_swiginit(self,
        cppModule.new_mapQuanTrait(loci=loc, *args, **kwargs))

new_mapQuanTrait.__doc__ = mapQuanTrait.__init__.__doc__
del mapQuanTrait.__init__
mapQuanTrait.__init__ = new_mapQuanTrait

def new_maQuanTrait(self, locus=-1, loci=[], wildtype=[0], sigma=[0], qtrait=[], *args, **kwargs):
    if locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        raise exceptions.TypeError('Please specify locus or loci')
    if type(wildtype) in [types.IntType, types.LongType]:
        wt = [wildtype]
    else:
        wt = wildtype
    if type(sigma) in [types.IntType, types.LongType, types.FloatType]:
        s = [sigma]*len(qtrait)
    elif len(sigma) == 1:
        s = sigma*len(qtrait)
    else:
        s = sigma
    cppModule.maQuanTrait_swiginit(self,
        cppModule.new_maQuanTrait(loci=loc, wildtype=wt, sigma=s, qtrait=qtrait, *args, **kwargs))

new_maQuanTrait.__doc__ = maQuanTrait.__init__.__doc__
del maQuanTrait.__init__
maQuanTrait.__init__ = new_maQuanTrait

def new_pyQuanTrait(self, locus=-1, loci=[], *args, **kwargs):
    if locus != -1 and type(locus) in [types.IntType, types.LongType]:
        loc = [locus]
    elif type(loci) in [types.IntType, types.LongType]:
        loc = [loci]
    elif type(loci) in [types.TupleType, types.ListType] and len(loci)>0:
        loc = loci
    else:
        raise exceptions.TypeError('Please specify locus or loci')
    cppModule.pyQuanTrait_swiginit(self,
        cppModule.new_pyQuanTrait(loci=loc, *args, **kwargs))

new_pyQuanTrait.__doc__ = pyQuanTrait.__init__.__doc__
del pyQuanTrait.__init__
pyQuanTrait.__init__ = new_pyQuanTrait


def new_rangeSplitter(self, range=[], ranges=[], *args, **kwargs):
    if len(range) == 0:
        rg = ranges
    else:
        rg = [range]
    for r in rg:
        if len(r) != 2:
            raise exceptions.TypeError("Expecting a range: %s" % r)
    cppModule.rangeSplitter_swiginit(self,
        cppModule.new_rangeSplitter(ranges=rg, *args, **kwargs))

new_rangeSplitter.__doc__ = rangeSplitter.__init__.__doc__
del rangeSplitter.__init__
rangeSplitter.__init__ = new_rangeSplitter


def new_genotypeSplitter(self, locus=None, loci=[],
    alleles=[], *args, **kwargs):
    if locus is not None:
        loc = [locus]
    else:
        loc = loci
    if len(alleles) == 0:
        raise exceptions.ValueError("Please specify alleles at each locus")
    if type(alleles[0]) in [type(0), type(0L)]:
        als = [alleles]
    else:
        als = alleles
    cppModule.genotypeSplitter_swiginit(self,
        cppModule.new_genotypeSplitter(loci=loc, alleles=als, *args, **kwargs))

new_genotypeSplitter.__doc__ = genotypeSplitter.__init__.__doc__
del genotypeSplitter.__init__
genotypeSplitter.__init__ = new_genotypeSplitter


def new_pedigreeMating(self, generator=None, pedigree=None, *args, **kwargs):
    if generator is None:
        generator = mendelianOffspringGenerator()
    cppModule.pedigreeMating_swiginit(self,
        cppModule.new_pedigreeMating(generator=generator, ped=pedigree, *args, **kwargs))

new_pedigreeMating.__doc__ = pedigreeMating.__init__.__doc__
del pedigreeMating.__init__
pedigreeMating.__init__ = new_pedigreeMating
%}
