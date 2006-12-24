/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu
 *                                                                         *
 *   $LastChangedDate$
 *   $Rev$                                                     *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

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

// for mac, or earlier version of gcc, _M_word_bit is used
// while _S_word_bit is used for newer versions
// to compile for mac, I define
//    #define _Bit_type_size _M_word_bit
// in config_mac.h. For other OS, we use _S_word_bit
#ifdef _MSC_VER

#define WORDBIT (8*sizeof(unsigned))
#define WORDTYPE unsigned
#define BITPTR(ref) ref._Myptr
#define BITOFF(ref) ref._Myoff

#define INT32 __int32

#else

#ifndef WORDBIT
#define WORDBIT std::_S_word_bit
#endif
#define WORDTYPE std::_Bit_type
#define BITPTR(ref) ref._M_p
#define BITOFF(ref) ref._M_offset

#define INT32 int32_t
#endif

#include <string>
using std::string;

/// This is not the most effecient method, but it is convenient to use.
#include "boost/lexical_cast.hpp"
#define toStr(val) (boost::lexical_cast<string>(val))

/// needed by the following typedefs
#include <vector>
using std::vector;

/// UINT should not be changed to unsigned long
/// since python extension use it as int.
typedef unsigned int UINT;

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif
// NOTE: the change of allele type here may need similar changes
// in the wrapper file simuPOP_common.i
#ifdef LONGALLELE
#ifdef _MSC_VER
// swig can not handle __int16 on windows yet, partly because
// stdint is not defined under windows
typedef unsigned  Allele;
typedef unsigned & AlleleRef;
#else
typedef uint16_t  Allele;
typedef uint16_t& AlleleRef;
#endif

#define AlleleInc(a)  ++(a)
#define AlleleDec(a)  --(a)
#define AlleleAdd(a, b) (a)+=(b)
#define AlleleMinus(a, b) (a)-=(b)
#define AlleleUnsigned(a) (a)

#else

#ifdef BINARYALLELE

typedef bool Allele;
typedef vector<Allele>::reference AlleleRef;
// bool type, inc go to 1, dec go to 0
#define AlleleInc(a)  (a)=true
#define AlleleDec(a)  (a)=false
#define AlleleAdd(a, b) (a)=((b)==0?(a):((b)>0?true:false))
#define AlleleMinus(a, b) (a)=((b)==0?(a):((b)>0?false:true))
#define AlleleUnsigned(a) (static_cast<unsigned char>(a))

#else

typedef unsigned char  Allele;
typedef unsigned char& AlleleRef;
#define AlleleInc(a)  ++(a)
#define AlleleDec(a)  --(a)
#define AlleleAdd(a, b) (a)+=(b)
#define AlleleMinus(a, b) (a)-=(b)
#define AlleleUnsigned(a) (a)
#endif
#endif

typedef std::vector<Allele>::iterator GenoIterator;
typedef std::vector<Allele>::const_iterator constGenoIterator;

// max allowed allele state
const unsigned long MaxAllele = std::numeric_limits<Allele>::max();
const unsigned long int MaxRandomNumber = std::numeric_limits<INT32>::max();

#define PopSWIGType "simuPOP::population *"
#define IndSWIGType "simuPOP::individual *"

enum Sex{ Male = 1, Female = 2};

// info is usually used for subpopulation index.
// signed short should be enough.
// if this is changed Info_Var_As_Numarray in utility.cpp also needs to be changed.
typedef double InfoType;
typedef std::vector<double>::iterator InfoIterator;
typedef std::vector<double>::const_iterator InfoConstIterator;
typedef signed short SubPopID;
const unsigned long MaxSubPopID = std::numeric_limits<SubPopID>::max();

typedef unsigned long ULONG;
typedef long LONG;

typedef std::vector<int>                   vectori;
typedef std::vector<double>                vectorf;
typedef std::vector<Allele>                vectora;
typedef std::vector<LONG>                  vectorl;
typedef std::vector<UINT>                  vectoru;
typedef std::vector<ULONG>                 vectorlu;
typedef std::vector<std::string>           vectorstr;
typedef std::vector<std::vector<int > >    intMatrix;
typedef std::vector<InfoType>              vectorinfo;
typedef std::vector<std::vector<double > > matrix;

#include <map>
using std::map;
typedef std::map<string, double>           strDict;
typedef std::map<int, double>              intDict;

namespace simuPOP
{
	/// exception handler. Exceptions will be passed to Python.
	class Exception
	{
		public:

			/// constructor
			/// \param msg error message
			Exception(string msg):m_msg(msg){}

			/// return error message
			const char* message(){ return m_msg.c_str(); }

			virtual ~Exception(){};

		private:
			/// error message
			string m_msg;
	};

	/// exception, thrown if out of memory
	class StopIteration: public Exception
	{
		public:
			StopIteration(string msg):Exception(msg){};
	};

	/// exception, thrown if out of memory
	class OutOfMemory: public Exception
	{
		public:
			OutOfMemory(string msg):Exception(msg){};
	};

	/// exception, thrown if file io failure
	class IOError: public Exception
	{
		public:
			IOError(string msg):Exception(msg){};
	};

	/// exception, thrown if index out of range
	class IndexError: public Exception
	{
		public:
			IndexError(string msg):Exception(msg){};
	};

	/// exception, thrown if type mismatch
	class TypeError: public Exception
	{
		public:
			TypeError(string msg):Exception(msg){};
	};

	/// exception, thrown if value of range etc
	class ValueError: public Exception
	{
		public:
			ValueError(string msg):Exception(msg){};
	};

	/// exception, thrown if system error occurs
	class SystemError: public Exception
	{
		public:
			SystemError(string msg):Exception(msg){};
	};

	// define DEBUG codes
	// DEbUG_CODE_LENGTH should be the number of debug codes
#define DBG_CODE_LENGTH 20

	enum DBG_CODE
	{
		DBG_ALL=0,
		DBG_GENERAL=1,
		DBG_UTILITY=2,
		DBG_OPERATOR=3,
		DBG_SIMULATOR=4,
		DBG_INDIVIDUAL=5,
		DBG_OUTPUTER=6,
		DBG_MUTATOR=7,
		DBG_RECOMBINATOR=8,
		DBG_INITIALIZER=9,
		DBG_POPULATION=10,
		DBG_STAT=11,
		DBG_TERMINATOR=12,
		DBG_TAGGER=13,
		DBG_VISUALIZER=14,
		DBG_SELECTOR=15,
		DBG_MATING=16,
		DBG_MIGRATOR=17,
		DBG_PROFILE=18,
		DBG_DEVEL=19
	};

	// DBG_NAMES are defined in utility.cpp
	// if you add any debug code here, you need to add its name
	// in utility.cpp as well. Otherwise, turnOnDebug, turnOffDebug
	// etc will be function well.

#define REP_ALL             -2
#define REP_CUR             -2
#define REP_LAST            -1
#define GRP_ALL             -1
#define PreMating            1
#define DuringMating         2
#define PostMating           4
	// combinations of mating scheme.
#define PrePostMating        (PreMating+PostMating)
#define PreDuringMating      (PreMating+DuringMating)
#define DuringPostMating     (DuringMating+PostMating)
#define PreDuringPostMating  (PreMating+DuringMating+PostMating)

	// used to set local variables
#define subPopVar_String(sp, var) (string("subPop[") + toStr(sp) + "]{\'" + var + "\'}")

	// standard library
#ifndef OPTIMIZED

#define DBG_ASSERT(cond, exception, message) \
if(!(cond)) \
{ \
	throw exception( \
	toStr(__FILE__)+toStr(":")+toStr(__LINE__)+toStr(" ")+message); \
}

#define DBG_FAILIF(cond, exception, message) \
if(cond) \
{ \
	throw exception( \
	toStr(__FILE__)+toStr(":")+toStr(__LINE__)+toStr(" ")+message); \
}

#define DBG_WARNING(cond, message) \
if(cond) \
{ \
	cout << "Warning (line " << __LINE__ << " in " << __FILE__ << "): " << message << endl; \
}

#define DBG_DO(dbgCode, expr) \
if(debug(dbgCode)){ expr; }

#define DBG_DO_( expr) expr

#else										  // optimized mode

#define DBG_ASSERT(cond, exception, message)
#define DBG_FAILIF(cond, exception, message)
#define DBG_WARNING(cond, message)
#define DBG_DO(dbgCode, expr)
#define DBG_DO_(expr)
#endif

				// definition for all mode

				// epsilon when during floating point comparison
#define cmp_epsilon (1.e-9)

#define CLEARFLAG(var) (var = 0)
#define SETFLAG(var, flag) (var |= flag)
#define RESETFLAG(var, flag) (var &= ~flag)
#define ISSETFLAG(var, flag) (!!(var & flag))

				// check range.
#define CHECKRANGEPLOIDY(p)  DBG_FAILIF( p>=ploidy(), IndexError, "index (" + toStr(p) + ") out of range of ploidy of 0 ~ " + toStr(ploidy()-1))
#define CHECKRANGESEX(sex) DBG_FAILIF( sex!=Male && sex!=Female, IndexError, "Wrong sex info. Male " + toStr(Male) + " or Fenamle "  + toStr(Female) + " only.")
#define CHECKRANGESUBPOP(subPop) DBG_FAILIF ( subPop >= numSubPop(), IndexError, "Subpop index (" + toStr(subPop) + ") out of range of 0  - " + toStr(numSubPop()-1))
#define CHECKRANGECHROM(chrom)   DBG_FAILIF ( chrom >= numChrom(), IndexError, "chromosome index (" + toStr(chrom) + ") out of range of 0 - " + toStr(numChrom()-1))
#define CHECKRANGELOCUS(chrom, locus) DBG_FAILIF( locus >= numLoci(chrom), IndexError, "locus index (" + toStr(chrom) + ") out of range of 0 - " + toStr(numLoci(chrom)-1))
#define CHECKRANGEABSLOCUS(locus) DBG_FAILIF(  locus >= totNumLoci(), IndexError, "absolute locus index (" + toStr(locus) + ") out of range of 0 - " + toStr(totNumLoci()-1))
#define CHECKRANGEGENOSIZE(p) DBG_FAILIF( p>=genoSize(),IndexError, "locus index  (" + toStr(p) + ") out of range of 0 - " + toStr(genoSize()-1))
#define CHECKRANGESUBPOPMEMBER(ind,sp) DBG_FAILIF( subPopSize(sp)>0 && ind >= subPopSize(sp), IndexError, "individual index (" + toStr(ind) + ") out of range 0 ~" + toStr(subPopSize(sp)-1) + " in subpopulation " + toStr(sp))
#define CHECKRANGEIND(ind) DBG_FAILIF(ind >= popSize(), IndexError, "individual index (" + toStr(ind) + ") is out of range of 0 ~ " + toStr(popSize()-1))
#define CHECKRANGEINFO(ind) DBG_FAILIF(ind >= infoSize(), IndexError, "indo index (" + toStr(ind) + ") is out of rage of 0 ~ " + toStr(infoSize()-1))

}
#endif
