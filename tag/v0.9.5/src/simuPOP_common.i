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
#include "terminator.h"
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
	%template()         vector< vector<string> >;
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

// individual and population are type names, and can not be used
// as function name. ind and pop are used instead.
// at the python level, individual and population are better.
%rename(individual) ind(ULONG);
%rename(individual) ind(ULONG, UINT);
%rename(population) pop(UINT) const;

%newobject LoadPopulation;
%newobject LoadSimulator;

%newobject simuPOP::population::extract;
%newobject simuPOP::population::clone;
%newobject simuPOP::simulator::extract;
%newobject simuPOP::simulator::clone;
%newobject simuPOP::baseOperator::clone;
%newobject simuPOP::mating::clone;
%newobject simuPOP::stat::clone;

// the following load a docstring file extracted from doxgen output.
// there will also be a bunch of %ignore directives as well
//
%include "simuPOP_doc.i"

%implicitconv floatList;
%implicitconv repList;

%implicitconv intList;
%implicitconv uintList;
%implicitconv uintListFunc;
%implicitconv floatListFunc;
%implicitconv stringList;
%implicitconv stringFunc;
%implicitconv opList;

%include "utility.h"
%include "genoStru.h"
%include "individual.h"

%implicitconv vspID;

namespace std
{
    %template()    vector<simuPOP::vspID >;
    %template()    vector<simuPOP::vspSplitter * >;
}

%implicitconv subPopList;

%include "virtualSubPop.h"
%include "population.h"

namespace std {
    %template()    vector<simuPOP::baseOperator * >;
}

%include "operator.h"

namespace std
{
    %template()    vector<PyObject*>;
    %template()    vector<simuPOP::homoMating * >;
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
%include "transmitter.h"
%include "selector.h"
%include "qtrait.h"
%include "penetrance.h"
%include "pedigree.h"

////////////////////////// SIMUPOP PYTHON UTILITY FUNCTIONS //////////////////////////

%pythoncode %{
import exceptions, types
# population.__init__ needs to check simuOptions
from simuOpt import simuOptions

# This wrapper function converts all the following inputs to
# a list of vspID so that the real constructor can correctly
# recognize them.
#
# subPops = 0
# subPops = [1, 2, 3]
# subPops = [(0, 1), (2, 1), 3]
#
def new_subPopList(self, subPops=[]):
    sp = []
    # in the case of subPops=0
    if type(subPops) not in [type([]), type(())]:
        sp = [vspID(subPops)]
    else:
        for s in subPops:
            sp.append(vspID(s))
    cppModule.subPopList_swiginit(self,
        cppModule.new_subPopList(sp))

new_subPopList.__doc__ = subPopList.__init__.__doc__
del subPopList.__init__
subPopList.__init__ = new_subPopList
    
def new_population(self, size=[], ploidy=2, loci=[], chromTypes=[],
    lociPos=[], ancGen=0, chromNames=[], alleleNames=[], lociNames=[],
    subPopNames=[], infoFields=[]):
    if type(size) in [type(0), type(0L)]:
        sp = [size]
    elif type(size) in [types.TupleType, types.ListType]:
        sp = size
    else:
        raise exceptions.TypeError('Parameter size should be an integer or a list of integers, %s given.' % size)
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
        cppModule.new_population(sp, ploidy, loci, chromTypes, ld,
            ancGen, chromNames, alleleNames, ln, subPopNames, infoFields))

new_population.__doc__ = population.__init__.__doc__
del population.__init__
population.__init__ = new_population


def new_extract(self, field=None, loci=None, infoFields=None, ancGen=-1, ped=None):
    removeInd = field is not None
    if field is None:
        field = ''
    removeLoci = loci is not None
    if loci is None:
        loci = []
    removeInfo = infoFields is not None
    if infoFields is None:
        infoFields = []
    if ped is None:
        return cppModule.population_extract(self, removeInd, field, removeLoci, loci, removeInfo,
            infoFields, ancGen)
    else:
        return cppModule.population_extract(self, removeInd, field, removeLoci, loci, removeInfo,
            infoFields, ancGen, ped)

if population.extract.__doc__ is not None:
    new_extract.__doc__ = population.extract.__doc__
population.extract = new_extract


def new_migrator(self, rate=[], *args, **kwargs):
    # parameter rate
    r = rate
    if type(rate) in [types.IntType, types.LongType, types.FloatType]:
        r = [[rate]]
    # if a single vector, [a,b] ==> [[a,b]]
    if type(rate) in [types.ListType, types.TupleType, types.FloatType] and \
        (len(rate) == 0 or type(rate[0]) in [types.IntType, types.LongType, types.FloatType]):
        r = [rate]
    cppModule.migrator_swiginit(self,
        cppModule.new_migrator(r, *args, **kwargs))

new_migrator.__doc__ = migrator.__init__.__doc__
del migrator.__init__
migrator.__init__ = new_migrator


def new_initByFreq(self, alleleFreq=[], *args, **kwargs):
    # parameter alleleFreq
    if len(alleleFreq) > 0 and type(alleleFreq[0]) in [types.IntType, types.LongType, types.FloatType]:
        af = [alleleFreq]
    else:
        af = alleleFreq
    cppModule.initByFreq_swiginit(self,
        cppModule.new_initByFreq(af, *args, **kwargs))

new_initByFreq.__doc__ = initByFreq.__init__.__doc__
del initByFreq.__init__
initByFreq.__init__ = new_initByFreq


def new_initByValue(self, value=[], *args, **kwargs):
    # parameter value
    if len(value) > 0 and type(value[0]) in [types.IntType, types.LongType]:
        val = [value]
    else:
        val = value
    cppModule.initByValue_swiginit(self,
        cppModule.new_initByValue(val, *args, **kwargs))

new_initByValue.__doc__ = initByValue.__init__.__doc__
del initByValue.__init__
initByValue.__init__ = new_initByValue

def new_stat(self, haploFreq=[], LD=[], LD_param={}, *args, **kwargs):
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
    # parameters of LD, convert 'stat':['LD', 'LD_prime'] etc to 'LD':True, 'LD_prime':True
    ldp = {}
    for key in LD_param.keys():
        if 'stat' == key and type(LD_param['stat']) in [types.TupleType, types.ListType]:
            for stat in LD_param['stat']:
                ldp[stat] = True
        else:
            ldp[key] = LD_param[key]
    cppModule.stat_swiginit(self,
        cppModule.new_stat(haploFreq=hf, LD=ld, LD_param=ldp,
            *args, **kwargs))

new_stat.__doc__ = stat.__init__.__doc__
del stat.__init__
stat.__init__ = new_stat



def new_genotypeSplitter(self, loci=[], alleles=[], *args, **kwargs):
    if len(alleles) == 0:
        raise exceptions.ValueError("Please specify alleles at each locus")
    if type(alleles[0]) in [type(0), type(0L)]:
        als = [alleles]
    else:
        als = alleles
    cppModule.genotypeSplitter_swiginit(self,
        cppModule.new_genotypeSplitter(loci, als, *args, **kwargs))

new_genotypeSplitter.__doc__ = genotypeSplitter.__init__.__doc__
del genotypeSplitter.__init__
genotypeSplitter.__init__ = new_genotypeSplitter


# def new_pedigreeMating(self, pedigree=None, generator=None, *args, **kwargs):
#     if generator is None:
#         generator = mendelianOffspringGenerator()
#     cppModule.pedigreeMating_swiginit(self,
#         cppModule.new_pedigreeMating(ped=pedigree, generator=generator, *args, **kwargs))
# 
# new_pedigreeMating.__doc__ = pedigreeMating.__init__.__doc__
# del pedigreeMating.__init__
# pedigreeMating.__init__ = new_pedigreeMating
%}
