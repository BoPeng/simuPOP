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

#include "config.h"
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

#include "sandbox.h"

%}

////////////////////////// DEFINE CARRAY //////////////////////////
%{
// templates are used to define array types for genotype and lineage
// because python types uses C functions, they need to be redefined with
// different names when templates are instantiated for these types.
#include "customizedTemplates.cpp"
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

// do not load these constants in config.h
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

%include "config.h"
//%include "mutant_vector.h"
%include "simuPOP_cfg.h"

typedef long ssize_t;

namespace std
{
    %template()         pair<size_t, size_t>; /* e.g. chromLocusPair */
    %template()         vector<string>;     /* e.g. infoFields */
    %template()         vector<size_t>;      /* e.g. subPopSizes */
    %template()         vector<double>;     /* e.g. lociPos */
    %template()         vector<long>;       /* e.g. vspID(vectori) */
    %template()         map<vector<long>, double>; /* e.g. MapSelector */
    %template()         vector<pair<size_t, size_t> >; /* e.g. mutant list */
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
    catch(std::bad_alloc)
    {
        SWIG_exception(SWIG_MemoryError, "Memory allocation error");
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

//%newobject simuPOP::Population::extract;
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

%implicitconv vspID;
%implicitconv intList;
%implicitconv uintString;
%implicitconv uintList;
%implicitconv lociList;
%implicitconv uintListFunc;
%implicitconv floatListFunc;
%implicitconv stringList;
%implicitconv intMatrix;
%implicitconv floatMatrix;
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

%rename(sb_RevertFixedSites)      simuPOP::sandbox::RevertFixedSites;
%rename(sb_MutSpaceSelector)      simuPOP::sandbox::MutSpaceSelector;
%rename(sb_MutSpaceMutator)       simuPOP::sandbox::MutSpaceMutator;
%rename(sb_MutSpaceRecombinator)  simuPOP::sandbox::MutSpaceRecombinator;

%include "sandbox.h"

