/**
 *  $File: simuPOP_cfg.h $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
   \file
   \brief system config file (typedef and macros)
 */

/// global configuration
#ifndef _SIMUPOP_CONFIG_H
#define _SIMUPOP_CONFIG_H

// include the sytem wide config.h
// Note: there is some strange #include problem if cmath is included later.
#include <cmath>

// for max allele etc.
#include <limits>

// under parent directory. Included with -I.. option.
#include "config.h"

// For the handling of binary modules and for the use of unordered_map in tr1
#ifdef _MSC_VER

// Visual C++ does not support threadprivate in DLL
#  define THREADPRIVATE_SUPPORT 0

#  define WORDBIT (8 * sizeof(unsigned))
#  define WORDTYPE unsigned
#  define BITPTR(ref) ref._Myptr
#  define BITOFF(ref) ref._Myoff
#  if _MSC_VER >= 1400
//   unordered_map in ''
#    define TR1_SUPPORT 1
#  else
#    define TR1_SUPPORT 0
#  endif
#else
#  define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#  if GCC_VERSION > 30400
#    define WORDBIT std::_S_word_bit
#  else
// previous version uses _M_word_bit
#    define WORDBIT std::_M_word_bit
#  endif
#  if GCC_VERSION > 40000
//   unordered_map in 'tr1/'
#    if defined(__INTEL_COMPILER)
//     intel c++ does not yet support tr1 because it cannot handle gcc tr1 headers such as type_traits
#      define TR1_SUPPORT 0
#    else
#      define TR1_SUPPORT 2
#    endif
#  else
#    define TR1_SUPPORT 0
#  endif
#  define WORDTYPE std::_Bit_type
#  define BITPTR(ref) ref._M_p
#  define BITOFF(ref) ref._M_offset

// Gcc 4.2 does not support threadprivate
#  if GCC_VERSION > 43000 || defined(__INTEL_COMPILER)
#    define THREADPRIVATE_SUPPORT 1
#  else
#    define THREAFPRIVATE_SUPPORT 0
#  endif
#endif

#ifdef _WIN64
/* Inclusion of windows.h is needed because InterlockedIncrement only accept LONGLONG
 * under windows 64 bit, which is only defined in this header file.
 * windows.h has its own definition of min and max, inclusion of windows.h will
 * cause problem with the use of std::min and std::max in the source code.
 * Definition of NOMINMAX before the inclusion of windows.h addresses this problem.
 */
#  define NOMINMAX
#  include <windows.h>
#  define ATOMICLONG LONGLONG
#elif defined(_WIN32)
#  define NOMINMAX
#  include <windows.h>
#  define ATOMICLONG LONG
#else
#  define ATOMICLONG unsigned long
#endif

/*
 * size_t string format of MSVC, GCC, and Intel compiler is different. The format of MSVC
 * is %lu, and others are %zu.
 */
#if defined (_WIN32) || defined (_WIN64)
#  define SIZE_T_FORMAT "%lu"
#else
#  define SIZE_T_FORMAT "%zu"
#endif

#include <string>
using std::string;

#define toID(val)  (static_cast<size_t>((val) + 0.5))

/// needed by the following typedefs
#include <vector>
using std::vector;


/// UINT should not be changed to unsigned long
/// since python extension use it as int.
typedef unsigned int UINT;

// for msvc, I have a portable stdint.h under win32 directory
#ifdef HAVE_STDINT_H
#  include <stdint.h>
#else
// for solaris, I have to use inttypes.h because there is no stdint.h
#  include <inttypes.h>
#endif

// NOTE: the change of allele type here may need similar changes
// in the wrapper file simuPOP_common.i
#ifdef LONGALLELE
//#  ifdef _MSC_VER
// according to MSVC manual, unsigned int is 32 bit, regardless of
// platform. Because uint32_t is not handled correctly by swig under
// windows, I am using unsigned int here.
typedef unsigned int Allele;
typedef unsigned int & AlleleRef;
//#  else
// under a unix platform, uint32_t seems to work fine.
//typedef uint32_t Allele;
//typedef uint32_t & AlleleRef;
//#  endif

#  define AlleleInc(a)  ++ (a)
#  define AlleleDec(a)  -- (a)
#  define AlleleAdd(a, b) (a) += (b)
#  define AlleleMinus(a, b) (a) -= (b)
#  define AlleleUnsigned(a) (a)
#  define ToAllele(a)    static_cast<Allele>(a)

#else

#  ifdef BINARYALLELE

typedef bool Allele;
typedef vector<Allele>::reference AlleleRef;
// bool type, inc go to 1, dec go to 0
#    define AlleleInc(a)  (a) = true
#    define AlleleDec(a)  (a) = false
#    define AlleleAdd(a, b) (a) = ((b) == 0 ? (a) : ((b) > 0 ? true : false))
#    define AlleleMinus(a, b) (a) = ((b) == 0 ? (a) : ((b) > 0 ? false : true))
#    define AlleleUnsigned(a) (static_cast<unsigned char>(a))
// this is used to avoid a VC++ C4800 warning when converting int to bool
#    define ToAllele(a)   ((a) != 0)

#  else
#    ifndef MUTANTALLELE
typedef unsigned char Allele;
typedef unsigned char & AlleleRef;
#    endif
#    define AlleleInc(a)  ++ (a)
#    define AlleleDec(a)  -- (a)
#    define AlleleAdd(a, b) (a) += (b)
#    define AlleleMinus(a, b) (a) -= (b)
#    define AlleleUnsigned(a) (a)
#    define ToAllele(a)   static_cast<Allele>(a)
#  endif
#endif

// for mutant vector -- the class wrapper for compressed vector
#include "mutant_vector.h"

#ifdef MUTANTALLELE
typedef simuPOP::mutant_vector<Allele>::iterator GenoIterator;
typedef simuPOP::mutant_vector<Allele>::const_iterator ConstGenoIterator;
#else
typedef std::vector<Allele>::iterator GenoIterator;
typedef std::vector<Allele>::const_iterator ConstGenoIterator;
#endif

// max allowed allele state
extern const unsigned long ModuleMaxAllele;
extern const unsigned long MaxRandomNumber;

#define PopSWIGType "simuPOP::Population *"
#define IndSWIGType "simuPOP::Individual *"

// For genotypic structure
enum Sex {
	MALE = 1,
	FEMALE = 2
};

// For genotypic structure
enum ChromType {
	CUSTOMIZED = 11,
	AUTOSOME = 12,
	CHROMOSOME_X = 13,
	CHROMOSOME_Y = 14
};

// For numOffspring and gene conversion
enum Distribution {
	CONSTANT = 21,
	BINOMIAL_DISTRIBUTION = 22,
	EXPONENTIAL_DISTRIBUTION = 23,
	GEOMETRIC_DISTRIBUTION = 24,
	POISSON_DISTRIBUTION = 25,
	UNIFORM_DISTRIBUTION = 26,
	NORMAL_DISTRIBUTION = 27,
	GAMMA_DISTRIBUTION = 28,
};

// For sexMode
enum SexMode {
	NO_SEX = 30,
	RANDOM_SEX = 31,
	PROB_OF_MALES = 32,
	NUM_OF_MALES = 33,
	NUM_OF_FEMALES = 34,
	SEQUENCE_OF_SEX = 35,
	GLOBAL_SEQUENCE_OF_SEX = 36,
};

// For gene conversion
enum ConversionMode {
	NO_CONVERSION = 41,
	NUM_MARKERS = 42,
	TRACT_LENGTH = 43
};

// For pedigree tracing
enum RelativeType {
	OFFSPRING = 50,                 // All offspring with all spouses (if there are more than one spouse)
	COMMON_OFFSPRING = 51,          // One spouse and their common offspring.
	SPOUSE = 52,                    // All spouses (with at least one offspring)
	OUTBRED_SPOUSE = 53,            // Spouse who is not a sibling.
	SIBLING = 54,                   // Siblings who share at least one parent
	FULLSIBLING = 55,               // Siblings who share two parents
};

// For pedigree tracing
enum SexChoice {
	ANY_SEX = 60,
	MALE_ONLY = 61,
	FEMALE_ONLY = 62,
	SAME_SEX = 63,
	OPPOSITE_SEX = 64
};

// this one does not use 71, 72 to allow users to specify affection
// status using True and False if they accidentally do so.
enum AffectionStatus {
	UNAFFECTED = 0,
	AFFECTED = 1,
	ANY_AFFECTION_STATUS = 2,
};

enum MultiLociMode {
	MULTIPLICATIVE = 81,
	ADDITIVE = 82,
	HETEROGENEITY = 83,
	EXPONENTIAL = 84,
};

enum MigrMode {
	BY_IND_INFO = 91,
	BY_PROBABILITY = 92,
	BY_PROPORTION = 93,
	BY_COUNTS = 94
};


//
enum InheritanceType {
	PATERNAL = 101,
	MATERNAL = 102,
	MEAN = 103,
	MAXIMUM = 104,
	MINIMUM = 105,
	SUMMATION = 106,
	MULTIPLICATION = 107,
};


#define DBG_CODE_LENGTH 21

/// CPPONLY
enum DBG_CODE {
	DBG_ALL = 0,
	DBG_GENERAL = 1,
	DBG_UTILITY = 2,
	DBG_POPULATION = 3,
	DBG_OPERATOR = 4,
	DBG_SIMULATOR = 5,
	DBG_INDIVIDUAL = 6,
	DBG_MUTATOR = 7,
	DBG_TRANSMITTER = 8,
	DBG_INITIALIZER = 9,
	DBG_STATOR = 10,
	DBG_TAGGER = 11,
	DBG_SELECTOR = 12,
	DBG_MATING = 13,
	DBG_MIGRATOR = 14,
	DBG_PROFILE = 15,
	DBG_BATCHTESTING = 16,
	DBG_INTEROPERABILITY = 17,
	DBG_COMPATIBILITY = 18,
	DBG_DEVEL = 19,
	DBG_WARNING = 20,
};

typedef unsigned char TraitIndexType;
extern const unsigned char MaxTraitIndex;

// info is usually used for subpopulation index.
// signed short should be enough.
// if this is changed Info_Var_As_Numarray in utility.cpp also needs to be changed.
typedef std::vector<double>::iterator InfoIterator;
typedef std::vector<double>::const_iterator ConstInfoIterator;
extern const size_t InvalidValue;

// FIXME: I need a type that is 32 or 64 bit long depending on platform
typedef unsigned long ULONG;
extern const size_t MaxIndexSize;
typedef long LONG;

typedef std::vector<long>                                vectori;
typedef std::vector<double>                              vectorf;
typedef std::vector<Allele>                              vectora;
typedef compressed_vector<unsigned int> 		 compressed_vectora;
typedef std::vector<size_t>                              vectoru;
typedef std::vector<std::string>                         vectorstr;
typedef std::pair<size_t, size_t>                        pairu;
typedef std::vector<std::vector<long> >                  matrixi;
typedef std::vector<std::vector<std::string > >          matrixstr;
typedef std::vector<std::vector<double > >               matrixf;

#include <map>
using std::map;
typedef std::map<string, double>           strDict;
// indDict is currently not defined in simuPOP_common.i
typedef std::map<int, double>              intDict;
typedef std::map<size_t, double>           uintDict;
typedef std::map<vectori, double>          tupleDict;

#define ValidPyObject(obj)   (obj != NULL && obj != Py_None)
#define InvalidPyObject(obj) (obj == NULL || obj == Py_None)

namespace simuPOP {


/// exception handler. Exceptions will be passed to Python.
class Exception
{
public:
	/// constructor
	/// \param msg error message
	Exception(const string & msg) : m_msg(msg)
	{
	}


	/// return error message
	const char * message()
	{
		return m_msg.c_str();
	}


	virtual ~Exception()
	{
	};

private:
	/// error message
	string m_msg;
};

/// exception, thrown if out of memory
class StopIteration : public Exception
{
public:
	StopIteration(const string msg) : Exception(msg)
	{
	};
};


/// exception, thrown if index out of range
class IndexError : public Exception
{
public:
	IndexError(const string msg) : Exception(msg)
	{
	};
};

/// exception, thrown if value of range etc
class ValueError : public Exception
{
public:
	ValueError(const string msg) : Exception(msg)
	{
	};
};

/// exception, thrown if system error occurs
class SystemError : public Exception
{
public:
	SystemError(const string msg) : Exception(msg)
	{
	};
};

/// exception, thrown if a runtime error occurs
class RuntimeError : public Exception
{
public:
	RuntimeError(const string msg) : Exception(msg)
	{
	};
};

/// exception, throw if an operator would like to stop
/// all replicates.
class StopEvolution : public Exception
{
public:
	StopEvolution(const string msg) : Exception(msg)
	{
	};
};


#define UnnamedSubPop        ""
}

/// This is not the most effecient method, but it is convenient to use.
// the pragma will suppress a lot of warnings when using toStr, generated
// from -Wconversion
// warning: conversion to char from int may alter its value
#pragma GCC system_header
#include "boost/lexical_cast.hpp"
#define toStr(val) (boost::lexical_cast < string > (val))

namespace simuPOP {

// standard library
#ifndef OPTIMIZED

// toStr(__LINE__) would trigger a icc warning about using temporary variate in function.
// a variable is therefore defined to hold __LINE__
#  define DBG_ASSERT(cond, exception, message) \
    if (!(cond)) \
	{ \
		int line = __LINE__; \
		const char * filename = strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : \
		                        (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__); \
		throw exception(filename + string(":") + toStr(line) + string(" ") + message); \
	}

#  define DBG_FAILIF(cond, exception, message) \
    if (cond) \
	{ \
		int line = __LINE__; \
		const char * filename = strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : \
		                        (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__); \
		throw exception(filename + string(":") + toStr(line) + string(" ") + message); \
	}

#  define PARAM_ASSERT DBG_ASSERT
#  define PARAM_FAILIF DBG_FAILIF

#  define DBG_WARNIF(cond, message) \
    if (debug(DBG_WARNING) && cond && !repeatedWarning(message)) \
	{ \
		cerr << "WARNING: " << message << endl; \
	}

#  define DBG_DO(dbgCode, expr) \
    if (debug(dbgCode)) { expr; }

#  define DBG_DO_(expr) expr

#else                                                                             // optimized mode

#  define DBG_ASSERT(cond, exception, message)
#  define DBG_FAILIF(cond, exception, message)

#  define PARAM_ASSERT(cond, exception, message) \
    if (!(cond)) \
	{ \
		int line = __LINE__; \
		const char * filename = strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : \
		                        (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__); \
		throw exception(filename + string(":") + toStr(line) + string(" ") + message); \
	}

#  define PARAM_FAILIF(cond, exception, message) \
    if (cond) \
	{ \
		int line = __LINE__; \
		const char * filename = strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : \
		                        (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__); \
		throw exception(filename + string(":") + toStr(line) + string(" ") + message); \
	}

#  define DBG_WARNIF(cond, message)
#  define DBG_DO(dbgCode, expr)
#  define DBG_DO_(expr)
#endif

// definition for all mode

// epsilon when during floating point comparison
#define cmp_epsilon (1.e-9)

#define CLEARFLAG(var) (var = 0)
#define SETFLAG(var, flag) (var |= flag)
#define RESETFLAG(var, flag) (var &= ~flag)
#define ISSETFLAG(var, flag) (!!(var & flag))

// check range.
#define CHECKRANGEPLOIDY(p)  DBG_FAILIF(p >= ploidy(), IndexError, "index (" + toStr(p) + ") out of range of ploidy of 0 ~ " + toStr(ploidy() - 1))
#define CHECKRANGESEX(sex) DBG_FAILIF(sex != MALE && sex != FEMALE, IndexError, "Wrong sex info. Only MALE or FEMALE is allowed.")
#define CHECKRANGESUBPOP(subPop) DBG_FAILIF(static_cast<UINT>(subPop) >= numSubPop(), IndexError, "Subpop index (" + toStr(subPop) + ") out of range of 0  ~ " + toStr(numSubPop() - 1))
#define CHECKRANGEVIRTUALSUBPOP(subPop) DBG_FAILIF(subPop != InvalidValue && static_cast<UINT>(subPop) >= numVirtualSubPop(), IndexError, "No virtual subpopulation is defined, or subpop index (" + toStr(subPop) + ") out of range of 0  ~ " + toStr(numVirtualSubPop() - 1))
#define CHECKRANGECHROM(chrom)   DBG_FAILIF(chrom >= numChrom(), IndexError, "chromosome index (" + toStr(chrom) + ") out of range of 0 ~ " + toStr(numChrom() - 1))
#define CHECKRANGELOCUS(chrom, locus) DBG_FAILIF(locus >= numLoci(chrom), IndexError, "locus index (" + toStr(locus) + ") out of range of 0 ~ " + toStr(numLoci(chrom) - 1))
#define CHECKRANGEABSLOCUS(locus) DBG_FAILIF(locus >= totNumLoci(), IndexError, "absolute locus index (" + toStr(locus) + ") out of range of 0 ~ " + toStr(totNumLoci() - 1))
#define CHECKRANGEGENOSIZE(p) DBG_FAILIF(p >= genoSize(), IndexError, "locus index  (" + toStr(p) + ") out of range of 0 ~ " + toStr(genoSize() - 1))
#define CHECKRANGESUBPOPMEMBER(ind, sp) DBG_FAILIF(subPopSize(sp) > 0 && ind >= subPopSize(sp), IndexError, "individual index (" + toStr(ind) + ") out of range 0 ~" + toStr(subPopSize(sp) - 1) + " in subpopulation " + toStr(sp))
#define CHECKRANGEIND(ind) DBG_FAILIF(ind >= popSize(), IndexError, "individual index (" + toStr(ind) + ") out of range of 0 ~ " + toStr(popSize() - 1))
#define CHECKRANGEINFO(ind) DBG_FAILIF(ind >= infoSize(), IndexError, "info index (" + toStr(ind) + ") out of range of 0 ~ " + toStr(infoSize() - 1))
}
#endif