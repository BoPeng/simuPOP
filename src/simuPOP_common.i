//
// $File: simuPOP_common.i $
// $LastChangedDate$
// $Rev$
//
// This file is part of simuPOP, a forward-time population genetics
// simulation environment. Please visit http://simupop.sourceforge.net
// for details.
//
// Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//

// for uint16_t
%include "stdint.i"

////////////////////////// INCLUDE FILES //////////////////////////
%{

#include "../config.h"
#include "simuPOP_cfg.h"
#include "genoStru.h"
#include "individual.h"
#include "population.h"
#include "pedigree.h"
#include "virtualSubPop.h"
#include "operator.h"
#include "simulator.h"

#include "utility.h"
#include "pedigree.h"
#include "initializer.h"
#include "outputer.h"
#include "mating.h"
#include "tagger.h"
#include "stator.h"
#include "migrator.h"
#include "mutator.h"
#include "transmitter.h"
#include "selector.h"
#include "qtrait.h"
#include "penetrance.h"

%}

////////////////////////// DEFINE CARRAY //////////////////////////
%{
extern "C"
{
#include "customizedTypes.c"
}
%}

%pythoncode %{
#redefine __repr__ to make it shorter.
def _swig_repr(self):
    if hasattr(self, 'describe'):
        return self.describe()
    else:
        return "<%s.%s>" % (self.__class__.__module__.split('.')[-1].split('_')[0], self.__class__.__name__)
%}



////////////////////////// CLEAN EXTRA SYMBOLS //////////////////////////

// do not load these constants in ../config.h
%ignore HAVE__BOOL;
%ignore HAVE_DECL_ACOSH;
%ignore HAVE_DECL_ASINH;
%ignore HAVE_DECL_LOG1P;
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
%ignore DBG_ALL;
%ignore DBG_GENERAL;
%ignore DBG_UTILITY;
%ignore DBG_POPULATION;
%ignore DBG_OPERATOR;
%ignore DBG_SIMULATOR;
%ignore DBG_INDIVIDUAL;
%ignore DBG_MUTATOR;
%ignore DBG_TRANSMITTER;
%ignore DBG_INITIALIZER;
%ignore DBG_STATOR;
%ignore DBG_TAGGER;
%ignore DBG_SELECTOR;
%ignore DBG_MATING;
%ignore DBG_MIGRATOR;
%ignore DBG_PROFILE;
%ignore DBG_BATCHTESTING;
%ignore DBG_INTEROPERABILITY;
%ignore DBG_COMPATIBILITY;
%ignore DBG_DEVEL;
%ignore DBG_CODE_LENGTH;
%ignore SIMUPOP_VAR_NAME;
%ignore SIMUPOP_VER;
%ignore PopSWIGType;
%ignore IndSWIGType;

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
    %template()         pair<UINT, UINT>;
    %template()         pair<UINT, ULONG>;
    %template()         vector<Allele>;
    %template()         vector<long int>;
    %template()         vector<ULONG>;
    %template()         vector<InfoType>;
    %template()         vector<double>;
    %template()         vector<string>;
    %template()         pair<string, double>;
    %template()         map<string, double>;
    %template()         map<int, int>;
    %template()         map<vector<long int>, double>;
    %template()         pair<ULONG, ULONG>;
    %template()         vector<pair<ULONG, ULONG> >;
    %template()         vector< vector<long int> >;
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
    catch(simuPOP::IndexError e)
    {
        SWIG_exception(SWIG_IndexError, e.message());
    }
    catch(simuPOP::ValueError e)
    {
        SWIG_exception(SWIG_ValueError, e.message());
    }
    catch(simuPOP::SystemError e)
    {
        SWIG_exception(SWIG_SystemError, e.message());
    }
    catch(simuPOP::RuntimeError e)
    {
        SWIG_exception(SWIG_RuntimeError, e.message());
    }
    catch(...)
    {
        SWIG_exception(SWIG_UnknownError, "Unknown runtime error happened.");
    }
}


////////////////////////// SIMUPOP TYPES //////////////////////////

%ignore std::operator<<(ostream&, const strDict&);
%ignore std::operator<<(ostream&, const intDict&);
%ignore simuPOP::IndAlleleIterator;
%ignore simuPOP::IndInfoIterator;

%newobject loadPopulation;

// %newobject simuPOP::Population::extract;
%newobject simuPOP::Population::extractSubPops;
%newobject simuPOP::Population::extractIndividuals;
%newobject simuPOP::Population::clone;
%newobject simuPOP::Simulator::extract;
%newobject simuPOP::Simulator::clone;
%newobject simuPOP::BaseOperator::clone;
%newobject simuPOP::MatingScheme::clone;
%newobject simuPOP::Stat::clone;

// the following load a docstring file extracted from doxgen output.
// there will also be a bunch of %ignore directives as well
//
%include "simuPOP_doc.i"

%implicitconv floatList;
%implicitconv repList;

%implicitconv intList;
%implicitconv uintString;
%implicitconv uintList;
%implicitconv uintListFunc;
%implicitconv floatListFunc;
%implicitconv stringList;
%implicitconv stringMatrix;
%implicitconv stringFunc;
%implicitconv lociList;
%implicitconv opList;

%include "utility.h"
%include "genoStru.h"
%include "individual.h"

%implicitconv vspID;

namespace std
{
    %template()    vector<simuPOP::vspID >;
    %template()    vector<simuPOP::BaseVspSplitter * >;
}

%implicitconv subPopList;

%include "virtualSubPop.h"
%include "population.h"

namespace std {
    // this place should be vector<const simuPOP::BaseOperator *> but
    // swig right now does not handle this well......q
    %template()    vector<simuPOP::BaseOperator * >;
}

%include "operator.h"

namespace std
{
    %template()    vector<PyObject*>;
    %template()    vector<simuPOP::HomoMating * >;
}


////////////////////////// SIMUPOP CLASSES //////////////////////////

%include "mating.h"
%include "simulator.h"
%include "stator.h"
%include "outputer.h"
%include "initializer.h"
%include "tagger.h"
%include "migrator.h"
%include "mutator.h"
%include "transmitter.h"
%include "selector.h"
%include "qtrait.h"
%include "penetrance.h"
%include "pedigree.h"

////////////////////////// SIMUPOP PYTHON UTILITY FUNCTIONS //////////////////////////

%pythoncode %{
import exceptions, types

# This constant is used by parameters loci, reps and subPops to 
# input a 'all available' case. Although both None and True represents
# this case, a devoted constant is easier to remember.
#
ALL_AVAIL = True
UNSPECIFIED = False

def evolve_pop(self, initOps=[], preOps=[], matingScheme=MatingScheme(), postOps=[],
    finalOps=[], gen=-1, dryrun=False):
    '''Evolve the current population <em>gen</em> generations using mating
    scheme <em>matingScheme</em> and operators <em>initOps</em> (applied before
    evolution), <em>preOps</em> (applied to the parental population at the
    beginning of each life cycle), <em>postOps</em> (applied to the offspring
    population at the end of each life cycle) and <em>finalOps</em> (applied at
    the end of evolution). More specifically, this function creates a
    <em>Simulator</em> using the current population, call its <em>evolve</em>
    function using passed parameters and then replace the current population
    with the evolved population. Please refer to function 
    <tt>Simulator.evolve</tt> for more details about each parameter.'''
    if dryrun:
        print describeEvolProcess(initOps, preOps, matingScheme, postOps, finalOps, gen, 1)
        return (0,)
    if isinstance(self, Pedigree):
        raise ValueError("Evolving a pedigree object directly is not allowed.")
    # create a simulator with self
    simu = Simulator(self)
    # evolve
    gen = simu.evolve(initOps, preOps, matingScheme, postOps, finalOps, gen)
    # get the evolved population
    self.swap(simu.population(0))
    return gen[0]

Population.evolve = evolve_pop

def all_individuals(self, subPops=ALL_AVAIL, ancGens=ALL_AVAIL):
    '''Return an iterator that iterat through all (virtual) subpopulations
    in all ancestral generations. A list of (virtual) subpopulations
    (<em>subPops</em>) and a list of ancestral generations (<em>ancGens</em>,
    can be a single number) could be specified to iterate through only
    selected subpopulation and generations.
    '''
    if ancGens is ALL_AVAIL:
        gens = range(self.ancestralGens() + 1)
    elif hasattr(ancGens, '__iter__'):
        gens = ancGens
    else:
        gens = [ancGens]
    #
    curGen = self.curAncestralGen()
    for gen in gens:
        self.useAncestralGen(gen)
        if subPops is ALL_AVAIL:
            for ind in self.individuals():
                yield ind
        else:
            for subPop in subPops:
                for ind in self.individuals(subPop):
                    yield ind
    self.useAncestralGen(curGen)

Population.allIndividuals = all_individuals

def as_pedigree(self, idField='ind_id', fatherField='father_id', motherField='mother_id'):
    '''Convert the existing population object to a pedigree. After this function
    pedigree function should magically be usable for this function.
    '''
    if isinstance(self, Pedigree):
        return
    ped = Pedigree(self, loci=ALL_AVAIL, infoFields=ALL_AVAIL, ancGens=ALL_AVAIL,
        idField=idField, fatherField=fatherField, motherField=motherField,
        stealPop=True)
    # swap ped and this object. (I do not know if this is the right thing to do)
    self.__class__, ped.__class__ = ped.__class__, self.__class__
    self.this, ped.this = ped.this, self.this

Population.asPedigree = as_pedigree

def as_population(self):
    '''Convert the existing pedigree object to a population. This function will
    behave like a regular population after this function call.'''
    if isinstance(self, Population):
        return
    pop = population(0)
    # the pedigree data has been swapped to pop
    pop.swap(self)
    # swap ped and this object. (I do not know if this is the right thing to do)
    self.__class__, pop.__class__ = pop.__class__, self.__class__
    self.this, pop.this = pop.this, self.this

Pedigree.asPopulation = as_population

def _new_Migrator(self, rate=[], *args, **kwargs):
    # parameter rate
    r = rate
    if type(rate) in [types.IntType, types.LongType, types.FloatType]:
        r = [[rate]]
    # if a single vector, [a,b] ==> [[a,b]]
    if type(rate) in [types.ListType, types.TupleType, types.FloatType] and \
        (len(rate) == 0 or type(rate[0]) in [types.IntType, types.LongType, types.FloatType]):
        r = [rate]
    cppModule.Migrator_swiginit(self,
        cppModule.new_Migrator(r, *args, **kwargs))

_new_Migrator.__doc__ = Migrator.__init__.__doc__
del Migrator.__init__
Migrator.__init__ = _new_Migrator


def _new_Stat(self, haploFreq=[], LD=[], *args, **kwargs):
    # parameter haploFreq
    if len(haploFreq) > 0 and type(haploFreq[0]) in [types.IntType, types.LongType]:
        hf = [haploFreq]
    else:
        hf = haploFreq
    # parameter LD
    if len(LD) > 0 and type(LD[0]) in [types.IntType, types.LongType]:
        ld = [LD]
    else:
        ld = LD
    cppModule.Stat_swiginit(self,
        cppModule.new_Stat(haploFreq=hf, LD=ld, *args, **kwargs))

_new_Stat.__doc__ = Stat.__init__.__doc__
del Stat.__init__
Stat.__init__ = _new_Stat

def _new_GenotypeSplitter(self, loci=[], alleles=[], *args, **kwargs):
    if len(alleles) == 0:
        raise exceptions.ValueError("Please specify alleles at each locus")
    if type(alleles[0]) in [type(0), type(0L)]:
        als = [alleles]
    else:
        als = alleles
    cppModule.GenotypeSplitter_swiginit(self,
        cppModule.new_GenotypeSplitter(loci, als, *args, **kwargs))

_new_GenotypeSplitter.__doc__ = GenotypeSplitter.__init__.__doc__
del GenotypeSplitter.__init__
GenotypeSplitter.__init__ = _new_GenotypeSplitter

%}
