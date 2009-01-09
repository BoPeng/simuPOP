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

#ifndef _UTILITY_H
#define _UTILITY_H

/**
   \file
   \brief head file for all utility classes/functions

   This file is very important since it handles most of the
   memory, file, IO, random number generator problems.

   The classes that are implemented here (also in utility.cpp)
   include:
 * SharedVariables
 * Expression
 * StreamElem
 * StreamManager
 * StreamProvider
 * RNG
 * BernulliTrials
 * parameter polymorphism
 */

extern "C" {
#include "Python.h"
}

#include "simuPOP_cfg.h"

#include <stdarg.h>
#include <iostream>
using std::ostream;
using std::iostream;
using std::cout;
using std::endl;

#include <utility>
using std::pair;

#include <algorithm>
using std::copy;
using std::find;

#include <iomanip>
using std::setw;

/// for ranr generator
#include "gsl/gsl_sys.h"                                           // for floating point comparison
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"

/// for bernulli trials.
// use vector<bool> instead of dynamic_bitset since I can manipulate
// bits directly in vector<bool>
typedef vector<bool> BitSet;

namespace simuPOP {

// ////////////////////////////////////////////////////////////
// / Debug and info functions
// ////////////////////////////////////////////////////////////

/// set debug codes. Default to turn on all debug codes. Only available in non-optimized modules.
void TurnOnDebug(DBG_CODE code = DBG_ALL);

void TurnOnDebug(string code);

/// turn off debug information. Default to turn off all debug codes. Only available in non-optimized modules.
void TurnOffDebug(DBG_CODE code = DBG_ALL);

#ifndef OPTIMIZED
/// test if one code is turned on, CPPONLY
bool debug(DBG_CODE code);

#endif

/// list all debug codes
void ListDebugCode();

/// dbg string for a code CPPONLY
string dbgString(DBG_CODE code);

#ifdef Py_REF_DEBUG
/// refcount debug
void saveRefCount();

/// check if refcount increase
void checkRefCount();

#endif

// ////////////////////////////////////////////////////////////
// / Some common functions/templates
// ////////////////////////////////////////////////////////////

/// a utility function to check keyboard stroke
/// CPPONLY
int simuPOP_kbhit();

/// CPPONLY
//  a utility function to get char unbuffered
int simuPOP_getch();

/// very important: re-define floating point comparison

#define fcmp_lt(a, b) (gsl_fcmp(a, b, cmp_epsilon) == -1)
#define fcmp_le(a, b) (gsl_fcmp(a, b, cmp_epsilon) <= 0)
#define fcmp_gt(a, b) (gsl_fcmp(a, b, cmp_epsilon) == 1)
#define fcmp_ge(a, b) (gsl_fcmp(a, b, cmp_epsilon) >= 0)
#define fcmp_eq(a, b) (gsl_fcmp(a, b, cmp_epsilon) == 0)
#define fcmp_ne(a, b) (gsl_fcmp(a, b, cmp_epsilon) != 0)
#define f_asBool(a)  (fcmp_ne(a, 0.))
}


namespace std {
/// CPPONLY how to output any std::pair
template<class T1, class T2>
ostream & operator<<(ostream & out, const pair<T1, T2> & pair)
{
	out << "(" << setw(3) << pair.first << "," << setw(3) << pair.second << ")";
	return out;
}


/// CPPONLY how to output any vector.
template<class T>
ostream & operator<<(ostream & out, const vector<T> & vec)
{
	if (!vec.empty() ) {
		typename vector<T>::const_iterator it = vec.begin();
		out << *it;
		for (++it; it != vec.end(); ++it)
			out << ", " << *it ;
	}
	return out;
}


/// CPPONLY how to output a dictionary
ostream & operator<<(ostream & out, const strDict & dict);

/// CPPONLY how to output a dictionary
ostream & operator<<(ostream & out, const intDict & dict);

/// CPPONLY can not use pow(3, n) because of overloading problem in msvc.
unsigned pow3(unsigned n);

}


namespace simuPOP {

/// CPPONLY operator to tell the affected status of an individual
template<class individual>
class isAffected
{
public:
	isAffected()
	{
	};

	bool operator()(const individual & ind)
	{
		return ind.affected();
	}


};

// ////////////////////////////////////////////////////////////
// / Shared variables
// ////////////////////////////////////////////////////////////

/// CPPONLY
void PyObj_As_Bool(PyObject * obj, bool & val);

/// CPPONLY
void PyObj_As_Int(PyObject * obj, int & val);

/// CPPONLY
void PyObj_As_Double(PyObject * obj, double & val);

/// CPPONLY
void PyObj_As_String(PyObject * obj, string & val);

/// CPPONLY
void PyObj_As_Array(PyObject * obj, vectorf & val);

/// CPPONLY
void PyObj_As_IntArray(PyObject * obj, vectori & val);

/// CPPONLY
void PyObj_As_Matrix(PyObject * obj, matrix & val);

/// CPPONLY
void PyObj_As_StrDict(PyObject * obj, strDict & val);

/// CPPONLY
void PyObj_As_IntDict(PyObject * obj, intDict & val);

/// CPPONLY
bool PyObj_Is_IntNumArray(PyObject * obj);

/// CPPONLY
bool PyObj_Is_DoubleNumArray(PyObject * obj);

/// CPPONLY
bool PyObj_Is_AlleleNumArray(PyObject * obj);

/// CPPONLY
PyObject * Double_Vec_As_NumArray(vectorf::iterator begin, vectorf::iterator end);

/// CPPONLY
PyObject * Allele_Vec_As_NumArray(GenoIterator begin, GenoIterator end);

/// CPPONLY
int NumArray_Size(PyObject * obj);

/// CPPONLY
char * NumArray_Data(PyObject * obj);

// ///////////////////////////////////////////////////////
/** CPPONLY shared variables.

   This class set and read Python variables using the given
   dictionary.

 */
class SharedVariables
{
public:
	/// CPPONLY
	SharedVariables() : m_dict(NULL), m_ownVars(false)
	{
	}


	/// CPPONLY
	SharedVariables(PyObject * dict, bool ownVars)
		: m_dict(dict), m_ownVars(ownVars)
	{
		// if not given a pre-existing dictionary, create one
		if (m_dict == NULL)
			m_dict = PyDict_New();

		DBG_ASSERT(m_dict != NULL, SystemError, "Can not create a new dictionary");

		if (!PyDict_Check(m_dict))
			throw SystemError("Invaid dictionary. The local namespace may have been cleared.");
	}


	/// CPPONLY, copy and increase ref count
	SharedVariables(const SharedVariables & rhs);

	/// CPPONLY swap two ShatedVariables
	void swap(SharedVariables & rhs)
	{
		std::swap(m_dict, rhs.m_dict);
		std::swap(m_ownVars, rhs.m_ownVars);
	}


	/// destructor
	/// I can not clear dict here since
	/// a resize of g_vars will copy this object and
	/// hence call this destructore.
	~SharedVariables();

	/// CPPONLY
	void clear()
	{
		PyDict_Clear(m_dict);
	}


	/// CPPONLY
	void setDict(PyObject * dict)
	{
		m_dict = dict;
	}


	/// CPPONLY set arbitrary value (type)
	PyObject * setVar(const string & name, const PyObject * val);

	/// CPPONLY get variable, NULL is returned if nothing is found.
	PyObject * getVar(const string & name, bool nameError = true);

	bool hasVar(const string & name)
	{
		// not null = has variable
		return PyDict_GetItemString(m_dict, const_cast<char *>(name.c_str())) != NULL;
	}


	/// remove variable
	void removeVar(const string & name);

	///CPPONLY
	PyObject * setBoolVar(const string & name, const bool val);

	///CPPONLY
	PyObject * setIntVar(const string & name, const int val);

	///CPPONLY
	PyObject * setDoubleVar(const string & name, const double val);

	///CPPONLY
	PyObject * setStringVar(const string & name, const string & val);

	///CPPONLY
	PyObject * setIntVectorVar(const string & name, const vectori & val);

	///CPPONLY
	PyObject * setDoubleVectorVar(const string & name, const vectorf & val);

	///CPPONLY
	PyObject * setStrDictVar(const string & name, const strDict & val);

	///CPPONLY
	PyObject * setIntDictVar(const string & name, const intDict & val);

	/// CPPONLY
	bool getVarAsBool(const string & name, bool nameError = true)
	{
		bool val;

		PyObj_As_Bool(getVar(name, nameError), val);
		return val;
	}


	/// CPPONLY
	int getVarAsInt(const string & name, bool nameError = true)
	{
		int val;

		PyObj_As_Int(getVar(name, nameError), val);
		return val;
	}


	/// CPPONLY
	double getVarAsDouble(const string & name, bool nameError = true)
	{
		double val;

		PyObj_As_Double(getVar(name, nameError), val);
		return val;
	}


	/// CPPONLY, obsolete
	// void numArrayOwnMem(const string& name);

	/// CPPONLY
	string getVarAsString(const string & name, bool nameError = true)
	{
		string val;

		PyObj_As_String(getVar(name, nameError), val);
		return val;
	}


	/// CPPONLY
	strDict getVarAsStrDict(const string & name, bool nameError = true)
	{
		strDict val;

		PyObj_As_StrDict(getVar(name, nameError), val);
		return val;
	}


	/// CPPONLY
	intDict getVarAsIntDict(const string & name, bool nameError = true)
	{
		intDict val;

		PyObj_As_IntDict(getVar(name, nameError), val);
		return val;
	}


	PyObject * & dict()
	{
		return m_dict;
	}


	/// CPPONLY, save m_dist as string
	string asString() const;

	void fromString(const string & vars);

private:
	/// the list
	PyObject * m_dict;

	/// whether or not the object owns the dictionary
	bool m_ownVars;

};

/** CPPONLY
 * get main dictionary (user namespace)
 */
SharedVariables & mainVars();

/** CPPONLY
 *  get module dictionary (it is different than mainDict!
 */
SharedVariables & moduleVars();

/// CPPONLY
PyObject * pyPopObj(void * p);

/// CPPONLY
PyObject * pyIndObj(void * p);

// ////////////////////////////////////////////////////////////
// Expression evaluation
// ////////////////////////////////////////////////////////////

/**  \brief expression

   CPPONLY

   This class evaluate python expressions.

 */
class Expression
{
public:
	Expression(const string & expr = "", const string & stmts = "",
		PyObject * locals = NULL)
		: m_expr(NULL), m_stmts(NULL), m_locals(locals)
	{
		if (m_locals == NULL)
			m_locals = mainVars().dict();

		// empty expression
		if (expr.empty() && stmts.empty())
			return;

		compileExpr(expr);
		compileStmts(stmts);
	}


	~Expression();

	/// Copy constructor, need to be defined because of ref count issue.
	/// CPPONLY
	Expression(const Expression & rhs);

	/** CPPONLY
	 * set local dictionary
	 */
	void setLocalDict(PyObject * dict)
	{
		m_locals = dict;
	}


	/// CPPONLY
	bool empty()
	{
		return m_expr == NULL && m_stmts == NULL;
	}


	/// CPPONLY
	void setExpr(const string & expr = "")
	{
		compileExpr(expr);
	}


	/// CPPONLY
	void setStmts(const string & stmts = "")
	{
		compileStmts(stmts);
	}


	/// CPPONLY evaluate with PyObject* output
	PyObject * evaluate();

	/// CPPONLY  return bool value
	bool valueAsBool();

	/// CPPONLY  return int value
	int valueAsInt();

	/// CPPONLY  return double value
	double valueAsDouble();

	/// CPPONLY  return string value
	string valueAsString();

	/// CPPONLY  return array value
	vectorf valueAsArray();

	/// CPPONLY  return dictionary value
	strDict valueAsStrDict();

	/// CPPONLY  return dictionary value
	intDict valueAsIntDict();

private:
	/// compile expression into byte code
	void compileExpr(const string & expr);

	/// compile statements into byte code
	void compileStmts(const string & stmts);

	/// compiled code
	PyObject * m_expr, * m_stmts;

	/// local namespace
	PyObject * m_locals;
};

// ////////////////////////////////////////////////////////////
// / Stream element, can be of different types
// ////////////////////////////////////////////////////////////

/// CPPONLY element of persistent stream
class StreamElem
{
public:
	/// type of persistent string:
	/// 1. output only file stream
	/// 2. read/write file stream
	/// 3. string stream for maximum performance.
	enum streamType { OFSTREAM, FSTREAM, SSTREAM };

	/** CPPONLY create a stream
	 * \param name filename
	 * \param readable iostream or just ostream
	 * \param realAppend whether or not keep old content when open an existing file
	 * \param useString use a stringstream rather than a file.
	 */
	StreamElem(const string & name, bool readable, bool realAppend, bool useString);

	/// CPPONLY copy constructor, called by map generator
	StreamElem(const StreamElem & rhs);

	/// destructor. Close stream and delete m_stream pointer.
	~StreamElem();

	/// CPPONLY if the file is opened as write-only, chage to read-write mode using this function
	void makeReadable();

	/// CPPONLY set realAppend status. clear existing file if necessary
	void makeAppend(bool append);

	/// CPPONLY return the stream
	ostream * stream()
	{
		return m_stream;
	}


	/// CPPONLY return type of the stream
	streamType type()
	{
		return m_type;
	}


	/** CPPONLY output info (filename, type and file status)
	 *  mostly for debug purposes
	 */
	string info();

	/// CPPONLY if this file is in append mode
	bool append()
	{
		return m_append;
	}


private:
	/// type of the stream
	streamType m_type;

	/// persistent across evolution etc.
	bool m_append;

	/// ostream pointer. Actual type might be ofstream, fstream or sstream
	ostream * m_stream;

	/// file name.
	string m_filename;
};

// ////////////////////////////////////////////////////////////
// / Stream Manager
// ////////////////////////////////////////////////////////////

/** CPPONLY stream manager:
   This class keeps track of persistent ostreams, aceesible by
   their names.
 */
class OstreamManager
{
public:
	/// CPPONLY
	OstreamManager();

	~OstreamManager();

	/** CPPONLY get an ostream pointer from a name.
	 * if the stream does not exist, create one and return.
	 */
	ostream * getOstream(const string & name, bool readable,  bool realAppend, bool useString);

	/** CPPONLY if persistant ostream exist for a filename
	 * this is mostly for debug purposes
	 */
	bool hasOstream(const string & filename);

	/// CPPONLY list all registered ostreams
	void listAll();

	/// CPPONLY close all files and clean the map
	void closeAll();

private:
	typedef map<string, StreamElem> ostreamMap;
	typedef map<string, StreamElem>::iterator ostreamMapIterator;
	typedef map<string, StreamElem>::value_type ostreamMapValue;

	/// collection of ostreams
	ostreamMap m_ostreams;
};

/// CPPONLY return a global Ostream Manager.
OstreamManager & ostreamManager();

// ////////////////////////////////////////////////////////////
// / OStream wrapper with our name extensions: >, >> etc
// ////////////////////////////////////////////////////////////

/** CPPONLY this class provides a unified IO interface for operators.

   1. the file stream is specified with the following syntax:

   file       the same as >file
   >file      default: output once
   >>file     append during this evolution session
   >>>file    append even acbdfhmstvross evolution
 | pipename in memory pipe (no output)
   ''         no output

   filename can have variables like %gen, %var so the files
   might be changing. This is the reason we can not use a
   fixed file strategy.

   2. operators use a simple interface to access files
   OStreamProvider.getStream()
   OStreamProdiver.closeStream()
   and m_flag operator
   OStreamProvider.noOutput()
   which means no input ("")

   3. ostream can be of the following types, but this is hidden to end users:
   ofstream     w   to a file, default
   fstream      r/w to a file, specified by >> and getStream(readable=true)
   sstream      r/w to a string, specifed by |

   This object will work with ostreamManager to keep track of persistent files.
 */
class StreamProvider
{
public:
	/// CPPONLY constructor. set the name parser
	StreamProvider(const string & output, const string & outputExpr);

	~StreamProvider()
	{
	}


	/// CPPONLY reset format string
	void setOutput(const string & output, const string & outputExpr);

	/// CPPONLY m_flag: if type ''
	bool noOutput()
	{
		return ISSETFLAG(m_flags, m_flagNoOutput);
	}


	/// CPPONLY get output stream. This function is not exposed to user.
	/**
	   get ostream.
	   - if this Operator uses cout return cout.
	   - if it does not use output, return something similar to /dev/null.
	   - if it uses a use-and-close file, create one and return its handle.
	   The file will be closed by \c closeOstream.
	   - if it uses a persistent file (>> or >>>), get from a global
	   repository of file handles. (If the file is not created yet,
	   the repository will create one.)
	   .
	   Note that if a use-and-close file is being opened in the repository,
	   the one from repository will be returned. This means that only
	   the first Operator that uses the file need to specify its
	   persistancy using >> or >>>

	   \param readable if the file need to be readable (other than writable).
	   The file has to be created by >>| or >>>| specifier.
	   \param gen current generation
	   \param rep calling replicate

	 */
	ostream & getOstream(PyObject * dict = NULL, bool readable = false);

	/// CPPONLY close ostream and delete ostream pointer.. if it is a ofstream.
	void closeOstream();

private:
	/// we need to set
	///  m_flagNoOutput
	///  m_flagUseDefault (when no file is specified.)
	///  m_flagAppend  (>> )
	///  m_flagRealAppend (>>>)
	///  m_flagUseString (|)
	void analyzeOutputString(const string & output);

	/// m_flag: if use and close type
	bool closeAfterUse()
	{
		return ISSETFLAG(m_flags, m_flagCloseAfterUse);
	}


private:
	/// plain text filename
	string m_filename;

	/// output filename parser
	Expression m_filenameExpr;

	/// internal m_flags of the operator. They are set during initialization for
	/// performance considerations.
	static const size_t m_flagNoOutput = 1;
	static const size_t m_flagUseDefault = 2;
	static const size_t m_flagAppend = 4;
	static const size_t m_flagRealAppend = 8;
	static const size_t m_flagCloseAfterUse = 16;
	static const size_t m_flagUseString = 32;
	static const size_t m_flagReadable = 64;

	/// m_flags
	unsigned char m_flags;

	/// file stream pointer, used ONLY if we create a new file stream
	ostream * m_filePtr;
};

// ////////////////////////////////////////////////////////////
// / Random number generator
// ////////////////////////////////////////////////////////////

/** \brief random number generator

   This random number generator class wraps around a number of
   random number generators from GNU Scientific Library. You can obtain
   and change system random number generator through the \c rng()
   function. Or create a separate random number generator and use
   it in your script.
 */
class RNG
{

public:
	/// Create a RNG object. You can also use \c rng() function to get the
	/// RNG used by simuPOP.
	RNG(const char * rng = NULL, unsigned long seed = 0);

	///
	~RNG();

	/// choose an random number generator, or set seed to the current RNG
	/**
	   \param rng name of the RNG. If rng is not given, environmental variable
	   GSL_RNG_TYPE will be used if it is available. Otherwise, RNG \c mt19937
	   will be used.
	   \param seed random seed. If not given, <tt>/dev/urandom</tt>,
	   <tt>/dev/random</tt>, system time will be used, depending on availability,
	   in that order. Note that windows system does not have \c /dev so system
	   time is used.
	 */
	void setRNG(const char * rng = NULL, unsigned long seed = 0);

	/// return RNG name
	const char * name()
	{
		return gsl_rng_name(m_RNG);
	}


	/// return the seed of this RNG
	unsigned long seed()
	{
		return m_seed;
	}


	/// return the maximum allowed seed value
	unsigned long maxSeed()
	{
		return std::numeric_limits<unsigned long>::max();
	}


	/// set random seed for this random number generator
	/// if seed is 0, method described in \c setRNG is used.
	void setSeed(unsigned long seed)
	{
		m_seed = seed;
		setRNG(name(), m_seed);
	}


	/// CPPONLY
	unsigned long generateRandomSeed();

	/// Maximum value of this RNG
	unsigned long max()
	{
		return gsl_rng_max(m_RNG);
	}


	string __repr__()
	{
		return toStr("<simuPOP::RNG ") + gsl_rng_name(m_RNG) + ">";
	}


	/// return a random number in the range of [0, 2, ... max()-1]
	unsigned long int randGet()
	{
		return gsl_rng_get(m_RNG);
	}


	bool randBit();

	/// return a random number in the range of [0, 1, 2, ... n-1]
	unsigned long int randInt(unsigned long int n)
	{
		return gsl_rng_uniform_int(m_RNG, n);
	}


	/// CPPONLY
	void randIntArray(ULONG n, ULONG size, ULONG * vec)
	{

		DBG_FAILIF(n <= 0, ValueError, "RandInt: n should be positive.");

		for (size_t i = 0; i < size; ++size)
			vec[i] = gsl_rng_uniform_int(m_RNG, n);
	}


	/// Geometric distribution
	int randGeometric(double p)
	{
		return gsl_ran_geometric(m_RNG, p);
	}


	/// Uniform distribution  [0,1)
	double randUniform01()
	{
		return gsl_rng_uniform(m_RNG);
	}


	/// Normal distribution
	double randNormal(double m, double v)
	{
		return gsl_ran_gaussian(m_RNG, v) + m;
	}


	double randExponential(double v)
	{
		return gsl_ran_exponential(m_RNG, v);
	}


	/// CPPONLY
	void randUniform01Array(ULONG size, double * vec)
	{
		for (size_t i = 0; i < size; ++size)
			vec[i] = gsl_rng_uniform(m_RNG);
	}


	/// Binomial distribution B(n, p)
	UINT randBinomial(UINT n, double p)
	{
		DBG_FAILIF(n <= 0, ValueError, "RandBinomial: n should be positive.");

		return gsl_ran_binomial(m_RNG, p, n);
	}


	/// Multinomial distribution
	void randMultinomial(unsigned int N, const vectorf & p, vectoru::iterator n)
	{
		// if sum p_i != 1, it will be normalized.
		// the size of n is not checked!
		gsl_ran_multinomial(m_RNG, p.size(), N, &p[0], & * n);
	}


	vectoru randMultinomialVal(unsigned int N, const vectorf & p)
	{
		// if sum p_i != 1, it will be normalized.
		// the size of n is not checked!
		vectoru val(p.size());

		gsl_ran_multinomial(m_RNG, p.size(), N, &p[0], &val[0]);
		return val;
	}


	/// Poisson distribution
	UINT randPoisson(double p)
	{
		return gsl_ran_poisson(m_RNG, p);
	}


	/// right hand side (single side) p-value for ChiSq value
	double pvalChiSq(double chisq, unsigned int df)
	{
		return 1 - gsl_cdf_chisq_P(chisq, df);
	}


private:
	/// global random number generator
	gsl_rng * m_RNG;

	/// seed used
	unsigned long m_seed;
};


// weighted sampling using Walker's alias algorithm
class weightedSampler
{
public:
	// set up AliasMethod table
	weightedSampler(RNG & rng, const vectorf & weight = vectorf(), bool fast = true)
		: m_RNG(&rng), m_q(0), m_a(0)                     // , m_fast(fast)
	{
		set(weight);
	};

	~weightedSampler()
	{
	};

	void set(const vectorf & weight);

	ULONG biSearch(double a)
	{
		unsigned long left = 0, right = m_N - 1;

		// a reasonable initial guess will improve the performance a little bit
		unsigned long mid = static_cast<ULONG>(a * (m_N - 1));

		// at most loop 4 times for n=31
		while (left <= right) {

			if (m_q[mid] > a) {
				if (mid == 0 || m_q[mid - 1] <= a)
					return mid;
				else
					right = mid - 1;
			} else {
				if (m_q[mid + 1] > a)
					return mid + 1;
				else
					left = mid + 1;
			}
			mid = (right + left) / 2;
		}
		// this should not be reached.
		return 0;
	}


	// sample without replacement from 0,...,n-1,
	// with weight freq
	ULONG get()
	{
		double rN = m_RNG->randUniform01() * m_N;

		size_t K = static_cast<size_t>(rN);

		rN -= K;

		if (rN < m_q[K])
			return K;
		else
			return m_a[K];
	}


	// sample without replacement from 0,...,n-1,
	// with weight freq
	ULONG get(vectorlu & res, ULONG shift = 0)
	{
		double rN;
		size_t K;

		for (vectorlu::iterator it = res.begin(); it != res.end(); ++it) {
			rN = m_RNG->randUniform01() * m_N;
			K = static_cast<size_t>(rN);

			rN -= K;

			if (rN < m_q[K])
				*it = K + shift;
			else
				*it = m_a[K] + shift;
		}
		return 0;
	}


	template<class Iterator>
	ULONG get(Iterator beg, Iterator end,  ULONG shift = 0)
	{
		double rN;
		size_t K;

		for (Iterator it = beg; it != end; ++it) {
			rN = m_RNG->randUniform01() * m_N;
			K = static_cast<size_t>(rN);

			rN -= K;

			if (rN < m_q[K])
				*it = static_cast<typename Iterator::value_type>(K + shift);
			else
				*it = static_cast<typename Iterator::value_type>(m_a[K] + shift);
		}
		return 0;
	}


	// print internal table
#ifndef OPTIMIZED
	vectorf q()
	{
		return m_q;
	}


	vectorlu a()
	{
		return m_a;
	}


#endif

private:
	/// pointer to a RNG
	RNG * m_RNG;

	/// length of weight.
	size_t m_N;

	/// internal table.
	vectorf m_q;

	vectorlu m_a;

	/// whether or not use fast (but use more RAM) method
	// bool m_fast;
};

/** this class encapsulate behavior of a sequence of Bernulli trial.
 *  the main idea is that when doing a sequence of Bernulli trials
 *  of the same probability, we can use much quicker algorithms
 *  instead of doing n Bernulli trials
 *
 *  For example, when N=10000, p=0.001. The usual way to do N Bin(p)
 *  trials is to do N randUnif(0,1)<p comparison.
 *
 *  using the new method, we can use geometric distrubution to find
 *  the next true event.
 *
 *  Also, for the cases of p=0.5, random bits are generated.
 *
 *  This class maintain a two dimensional table:
 *  a vector of probabilities cross expected number of trials
 *
 *            p1 p2 p3 p4 p5
 *  trial 1
 *  trial 2
 *  ...
 *  trial N
 *
 *  We expect that N is big (usually populaiton size) and p_i are small
 *
 *  using fast BernulliTrial method for fix p,
 *  we can fill up this table very quickly column by column
 *
 *  This class will provide easy access to row (each trial) or column
 *  (called each prob) of this table.
 *
 *  if this table is accessed row by row (each trial), a internal index
 *  is used.
 *
 *  if index exceeds N, trials will be generated all again.
 *  if trial will be called, e.g., N+2 times all the time,
 *  this treatment might not be very efficient.
 */
class BernulliTrials
{
public:
	/// CPPONLY
	BernulliTrials(RNG & rng);

	///
	BernulliTrials(RNG & rng, const vectorf & prob, ULONG trials);

	///
	~BernulliTrials();

	/** CPPONLY
	 * return size of trial
	 */
	ULONG trialSize() const
	{
		return m_N;
	}


	size_t probSize() const
	{
		return m_prob.size();
	}


	/// CPPONLY
	void setParameter(const vectorf & prob, ULONG trials);

	/// generate the trial table, reset m_cur
	void doTrial();

	/// CPPONLY
	UINT curTrial();

	/// if necessary, do trail again.
	void trial();

	bool trialSucc(size_t idx) const;

	bool trialSucc(size_t idx, size_t cur) const;

	// first and next succ across prob
	size_t probFirstSucc() const;

	size_t probNextSucc(size_t pos) const;

	// first and next succ across trial
	size_t trialFirstSucc(size_t idx) const;

	size_t trialNextSucc(size_t idx, size_t pos) const;

	void setTrialSucc(size_t idx, bool succ);

	/// return the succ rate for one index, used for verification pruposes
	double trialSuccRate(UINT index) const;

	/// return the succ rate for current trial, used for verification pruposes
	double probSuccRate() const;

	/// CPPONLY
	vectorf probabilities()
	{
		return m_prob;
	}


public:
	static const size_t npos = static_cast<size_t>(-1);

private:
	void setAll(size_t idx, bool v);

private:
	/// pointer to a random number generator.
	/// this is in preparation for multiple thread/RNG.
	RNG * m_RNG;

	/// number of trials.
	ULONG m_N;

	/// vector of probabilities
	vectorf m_prob;

	/// vectors to save result
	/// note that the result will be stored in the reversed order
	/// i.e.
	/// p1 last_succ_index .... first_succ_index
	/// p2
	/// p3 last_succ_index .... first_succ_index
	///
	/// then for each trial, we only need to compare m_cur with the last element
	/// and removing them if match.
	vector< BitSet > m_table;

	/// cache the actual point m_table[i].begin()._M_p and
	/// access bits through this pointer. This is much faster
	/// than using the reference interface.
	vector<WORDTYPE *> m_pointer;

	/// current trial. Used when user want to access the table row by row
	size_t m_cur;
};

/// return the currently used random number generator
RNG & rng();

/// set random number generator. If <tt>seed=0</tt> (default), a random seed will be given. If <tt>rng=""</tt>, seed will be set to the current random number generator.
void SetRNG(const string rng = "", unsigned long seed = 0);

/// list the names of all available random number generators
vectorstr AvailableRNGs();

// ////////////////////////////////////////////////////////////
// /  Global debug and initialization related functions
// ////////////////////////////////////////////////////////////

/// return the revision number of this simuPOP module. Can be used to test if a feature is available.
int simuRev();

/// return the version of this simuPOP module
string simuVer();

/// return the compiler used to compile this simuPOP module
string ModuleCompiler();

/// return the date when this simuPOP module is compiled
string ModuleDate();

/// return the Python version this simuPOP module is compiled for
string ModulePyVersion();

/// return the platform on which this simuPOP module is compiled
string ModulePlatForm();

#ifdef BINARYALLELE
// efficiently copy alleles (block by block, rather than 1 by 1)
/// CPPONLY
void copyGenotype(GenoIterator fr, GenoIterator to, size_t n);

void clearGenotype(GenoIterator to, size_t n);

#  ifndef OPTIMIZED
void testCopyGenotype();

#  endif
#endif

/// CPPONLY initialize module simuPOP when using "import simuPOP"
bool initialize();

/// return \c True if this simuPOP module is optimized
bool Optimized();

/// print out system limits
void Limits();

/// return the allele type of the current module. Can be \c binary, \c short, or \c long.
string AlleleType();

/** return the maximum allowed allele state of the current simuPOP module,
 *  which is \c 1 for binary modules, \c 255 for short modules and \c 65535
 *  for long modules.
 *  <group>allele</group>
 */
ULONG MaxAllele();

/// CPPONLY get a null stream that discard everything
ostream & cnull();

/// set the standard output (default to standard Python output)
void setLogOutput(const string filename = "");

// ////////////////////////////////////////////////////////////
// / Parameter polymorphism
// ////////////////////////////////////////////////////////////

/** A wrapper to a python function
 *  CPPONLY
 */
class pyObject
{
public:
	pyObject(PyObject * obj) : m_object(NULL)
	{
		if (obj != NULL && obj != Py_None)
			m_object = obj;

		Py_XINCREF(m_object);
	}


	~pyObject()
	{
		Py_XDECREF(m_object);
	}


	pyObject(const pyObject & rhs) : m_object(rhs.m_object)
	{
		Py_XINCREF(m_object);
	}


	PyObject * object() const
	{
		return m_object;
	}


	bool isValid() const
	{
		return m_object != NULL;
	}


private:
	PyObject * m_object;
};


/** A wrapper to a python function
 *  CPPONLY
 */
class pyFunc
{
public:
	pyFunc(PyObject * func) : m_func(func)
	{
		DBG_ASSERT(!m_func.isValid() || PyCallable_Check(m_func.object()),
			ValueError,
			"Passed parameter should be None or a Python function");
	}


	bool isValid() const
	{
		return m_func.isValid();
	}


	PyObject * func() const
	{
		return m_func.object();
	}


	// Note how ... are passed to Py_BuildValue using Py_VaBuildValue
	template <typename T>
	T operator()(void converter(PyObject *, T &), char * format, ...)
	{
		va_list argptr;

		va_start(argptr, format);
		PyObject * arglist = Py_VaBuildValue(format, argptr);
		va_end(argptr);
		PyObject * pyResult = PyEval_CallObject(m_func.object(), arglist);

		Py_XDECREF(arglist);
		if (pyResult == NULL) {
#ifndef OPTIMIZED
			if (debug(DBG_GENERAL)) {
				PyErr_Print();
				PyErr_Clear();
			}
#endif
			throw ValueError("Function call failed.\n");
		}
		T retValue;
		converter(pyResult, retValue);
		Py_DECREF(pyResult);
		return retValue;
	}


	PyObject * operator()(char * format, ...)
	{
		va_list argptr;

		va_start(argptr, format);
		PyObject * arglist = Py_VaBuildValue(format, argptr);
		va_end(argptr);
		PyObject * pyResult = PyEval_CallObject(m_func.object(), arglist);

		Py_XDECREF(arglist);
		if (pyResult == NULL) {
#ifndef OPTIMIZED
			if (debug(DBG_GENERAL)) {
				PyErr_Print();
				PyErr_Clear();
			}
#endif
			throw ValueError("Function call failed\n");
		}
		return pyResult;
	}


private:
	pyObject m_func;
};

/** This class defines an interface using which both a integer number and
 *  a list of numbers can be accpted.
 */
template <typename T>
class typeList
{
public:
	typeList(const vector<T> & values = vector<T>()) : m_elems(values)
	{
	}


	typeList(T value) : m_elems(1, value)
	{
	}


	T operator[](size_t i) const
	{
		DBG_FAILIF(i >= size(), IndexError, "Index out of range");
		return m_elems[i];
	}


	bool empty() const
	{
		return m_elems.empty();
	}


	size_t size() const
	{
		return m_elems.size();
	}


	vector<T> & elems()
	{
		return m_elems;
	}


protected:
	vector<T> m_elems;
};


typedef typeList<int> intList;

// I cannot use template here because otherwise SWIG does not
// handle the type correctly. I guess this can be my problem with using
// of template in simuPOP_common.i

class uintList
{
public:
	typedef vectorlu::iterator iterator;
	typedef vectorlu::const_iterator const_iterator;

public:
	uintList(const vectorlu & values = vectorlu()) : m_elems(values)
	{
	}


	uintList(ULONG value) : m_elems(1, value)
	{
	}


	ULONG operator[](size_t i) const
	{
		DBG_FAILIF(i >= size(), IndexError, "Index out of range");
		return m_elems[i];
	}


	bool empty() const
	{
		return m_elems.empty();
	}


	size_t size() const
	{
		return m_elems.size();
	}


	vectorlu & elems()
	{
		return m_elems;
	}


	const_iterator begin() const
	{
		return m_elems.begin();
	}


	const_iterator end() const
	{
		return m_elems.end();
	}


protected:
	vectorlu m_elems;
};


class floatList
{
public:
	floatList(const vectorf & values = vectorf()) : m_elems(values)
	{
	}


	floatList(double value) : m_elems(1, value)
	{
	}


	double operator[](size_t i) const
	{
		DBG_FAILIF(i >= size(), IndexError,
			"Index " + toStr(i) + " out of range of 0 ~ " + toStr(size() - 1));
		return m_elems[i];
	}


	bool empty() const
	{
		return m_elems.empty();
	}


	size_t size() const
	{
		return m_elems.size();
	}


	vectorf & elems()
	{
		return m_elems;
	}


protected:
	vectorf m_elems;
};


class uintListFunc : public uintList
{
public:
	uintListFunc(const vectorlu & values = vectorlu()) :
		uintList(values), m_func(NULL)
	{
	}


	uintListFunc(ULONG value) : uintList(value), m_func(NULL)
	{
	}


	uintListFunc(PyObject * func) : uintList(), m_func(func)
	{
	}


	pyFunc func() const
	{
		return m_func;
	}


private:
	pyFunc m_func;
};


class floatListFunc : public floatList
{
public:
	floatListFunc(const vectorf & values = vectorf()) :
		floatList(values), m_func(NULL)
	{
	}


	floatListFunc(double value) : floatList(value), m_func(NULL)
	{
	}


	floatListFunc(PyObject * func) : floatList(), m_func(func)
	{
	}


	pyFunc func() const
	{
		return m_func;
	}


private:
	pyFunc m_func;
};


/** A class to specify replicate list. The reason why I cannot simple
 *  use vectori() is that users have got used to use a single number
 *  to specify a single replicate.
 */
class repList : public intList
{
public:
	repList(const vectori & reps = vectori()) :
		intList(reps)
	{
	}


	repList(int rep) : intList(rep)
	{
	}


	bool match(UINT rep, const vector<bool> & activeRep);

};


}
#endif
