/**
 *  $File: utility.h $
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
 * Bernullitrials
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
using std::cerr;
using std::endl;

#include <utility>
using std::pair;

#include <algorithm>
using std::copy;
using std::find;

#include <iomanip>
using std::setw;

#include <set>

/// for ranr generator
#include "gsl/gsl_sys.h"                                           // for floating point comparison
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
/// for openMP
#ifdef _OPENMP
#  include "omp.h"
#  if defined(__INTEL_COMPILER)
#    include "tbb/parallel_sort.h"
#    include "tbb/task_scheduler_init.h"
#  elif defined(GCC_VERSION) && GCC_VERSION >= 40300
#    include "parallel/algorithm"
#  endif
#endif
/// for bernulli trials.
// use vector<bool> instead of dynamic_bitset since I can manipulate
// bits directly in vector<bool>
typedef vector<bool> BitSet;

// this is used to compare loci positions when loci are provided by
// (chr, pos) pair.
#if defined (_WIN32) || defined (_WIN64)
// windows VC does not have round function
#define PRECISION(f) (double(size_t(f)) == f ? f : (size_t((f) * 100000. + 0.5) / 100000.))
#else
#define PRECISION(f)  (round(f*100000.)/100000.)
#endif


namespace simuPOP {

// ////////////////////////////////////////////////////////////
// / Debug and info functions
// ////////////////////////////////////////////////////////////

/** Set debug code \e code. More than one code could be specified using a comma
 *  separated string. Name of available codes are available from
 *  <tt>moduleInfo()['debug'].keys()</tt>.
 */
void turnOnDebug(const string & code = string());

/** Turn off debug code \e code. More than one code could be specified using a
 *  comma separated string. Default to turn off all debug codes.
 */
void turnOffDebug(const string & code = "DBG_ALL");

#ifndef OPTIMIZED
/// CPPONLY test if one code is turned on, CPPONLY
bool debug(DBG_CODE code);

/// CPPONLY return true if the message is repeated (so that it will not be printed)
bool repeatedWarning(const string & message);

/// CPPONLY
void initClock();

void elapsedTime(const string & name);

#else
#  define initClock();
#  define elapsedTime(name);
#endif

#ifdef Py_REF_DEBUG
/// refcount debug
void saveRefCount();

/// check if refcount increase
void checkRefCount();

#endif


// ////////////////////////////////////////////////////////////
// / Some common functions/templates
// ////////////////////////////////////////////////////////////

/** First argument is to set number of thread in openMP. The number of threads can be be positive,
 *  integer (number of threads) or 0, which implies all available cores, or
 *  a number set by environmental variable \c OMP_NUM_THREADS.
 *  Second and third argument is to set the type or seed of existing random number generator using RNG \e name
 *  with \e seed. If using openMP, it sets the type or seed of random number
 *  generator of each thread.
 */
void setOptions(const int numThreads = -1, const char * name = NULL, unsigned long seed = 0);

/// CPPONLY get number of thread in openMP
UINT numThreads();

/// CPPONLY return val and increase val by 1, ensuring thread safety
ATOMICLONG fetchAndIncrement(ATOMICLONG * val);

/// CPPONLY parallel sort by using tbb or gnu parallel
template<class T1, class T2>
void parallelSort(T1 start, T1 end, T2 cmp)
{
	if (numThreads() > 1) {
#ifdef _OPENMP
#  if defined(__INTEL_COMPILER)
		tbb::task_scheduler_init init(numThreads());
		tbb::parallel_sort(start, end, cmp);
#  elif defined(GCC_VERSION) && GCC_VERSION >= 40300
		__gnu_parallel::sort(start, end, cmp);
#  else
		std::sort(start, end, cmp);
#  endif
#endif
	} else {
		std::sort(start, end, cmp);
	}
}


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
	if (!vec.empty()) {
		typename vector<T>::const_iterator it = vec.begin();
		out << *it;
		for (++it; it != vec.end(); ++it)
			out << ", " << *it ;
	}
	return out;
}


/// CPPONLY how to output any vector.
template<typename T1, typename T2>
ostream & operator<<(ostream & out, const map<T1, T2> & dict)
{
	if (!dict.empty()) {
		typename map<T1, T2>::const_iterator it = dict.begin();
		typename map<T1, T2>::const_iterator it_end = dict.end();
		out << it->first << ": " << it->second;
		for (++it; it != it_end; ++it)
			out << ", " << it->first << ": " << it->second;
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

// ////////////////////////////////////////////////////////////
// / Parameter polymorphism
// ////////////////////////////////////////////////////////////

/** A wrapper to a python function
 *  CPPONLY
 */
class pyObject
{
public:
	pyObject(PyObject * obj, bool defToNone = false) : m_object(obj)
	{
		if (m_object == NULL && defToNone)
			m_object = Py_None;

		Py_XINCREF(m_object);
	}


	~pyObject()
	{
		Py_XDECREF(m_object);
	}


	pyObject & operator=(const pyObject & rhs)
	{
		// this should be the combination of ~pyObject and
		// pyObject(pyObject)
		Py_XDECREF(m_object);
		m_object = rhs.m_object;
		Py_XINCREF(m_object);
		return *this;
	}


	// when a object is copied, its ref increased. This is
	// a constructor so its own m_object does not exist
	// yet.
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
	pyFunc(PyObject * func);

	/// return number of arguments this function accepts.
	/// This function does not count tuple parameters.
	size_t numArgs() const
	{
		return m_numArgs;
	}


	bool hasArg(const string & arg) const
	{
		return std::find(m_args.begin(), m_args.end(), arg) != m_args.end();
	}

	string name() const
	{
		return m_name;
	}


	string arg(size_t arg) const
	{
		return m_args[arg];
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
	T operator()(void converter(PyObject *, T &), const char * format, ...) const
	{
		va_list argptr;

		va_start(argptr, format);
		PyObject * arglist = Py_VaBuildValue(const_cast<char *>(format), argptr);
		va_end(argptr);
		PyObject * pyResult = PyEval_CallObject(m_func.object(), arglist);

		Py_XDECREF(arglist);
		if (pyResult == NULL) {
			PyErr_Print();
			PyErr_Clear();
			throw ValueError("Function call failed.\n");
		}
		T retValue;
		converter(pyResult, retValue);
		Py_DECREF(pyResult);
		return retValue;
	}


	template <typename T>
	T operator()(void converter(PyObject *, T &), PyObject * arglist) const
	{
		PyObject * pyResult = PyEval_CallObject(m_func.object(), arglist);

		if (pyResult == NULL) {
			PyErr_Print();
			PyErr_Clear();
			throw ValueError("Function call failed.\n");
		}
		T retValue;
		converter(pyResult, retValue);
		Py_DECREF(pyResult);
		return retValue;
	}


	PyObject * operator()(const char * format, ...) const
	{
		va_list argptr;

		va_start(argptr, format);
		PyObject * arglist = Py_VaBuildValue(const_cast<char *>(format), argptr);
		va_end(argptr);
		PyObject * pyResult = PyEval_CallObject(m_func.object(), arglist);

		Py_XDECREF(arglist);
		if (pyResult == NULL) {
			PyErr_Print();
			PyErr_Clear();
			throw ValueError("Function call failed\n");
		}
		return pyResult;
	}


	PyObject * operator()(PyObject * args) const
	{
		PyObject * pyResult = PyEval_CallObject(m_func.object(), args);

		if (pyResult == NULL) {
			PyErr_Print();
			PyErr_Clear();
			throw ValueError("Function call failed\n");
		}
		return pyResult;
	}


private:
	pyObject m_func;

	string m_name;

	size_t m_numArgs;

	vectorstr m_args;

	PyObject * m_circular_self;
};


/// CPPONLY
/// Remove circular references to release memory of objects that 
/// are derived from PyObject
void cleanupCircularRefs();

/// CPPONLY
class pyGenerator
{
public:
	pyGenerator(PyObject * gen = NULL) : m_generator(NULL), m_iterator(NULL)
	{
		set(gen);
	}


	bool isValid()
	{
		return m_generator != NULL;
	}


	void set(PyObject * gen);

	~pyGenerator()
	{
		set(NULL);
	}


	PyObject * next();

private:
	PyObject * m_generator;
	PyObject * m_iterator;
};


/** A class to specify replicate list. The reason why I cannot simple
 *  use vectori() is that users have got used to use a single number
 *  to specify a single replicate.
 */
class intList
{
public:
	intList(PyObject * obj = NULL);

	/// CPPONLY
	intList(const vectori & reps) :
		m_elems(reps), m_allAvail(false)
	{
	}


	/// CPPONLY
	vectori elems() const
	{
		return m_elems;
	}


	/// CPPONLY
	bool allAvail()
	{
		return m_allAvail;
	}


	/// CPPONLY
	bool match(ssize_t rep, const vector<bool> * activeRep = NULL) const;

private:
	vectori m_elems;
	bool m_allAvail;
};


// At the CPP level, uintList() is ALL_AVAIL, uintList(vectoru()) is regular empty, uintList(NULL) is UNSPECIFIED
// At the Python level True is ALL_AVAIL, a list if regular, False is UNSPECIFIED
class uintList
{
public:
	typedef vectoru::iterator iterator;
	typedef vectoru::const_iterator const_iterator;

private:
	enum listStatus {
		REGULAR = 0,
		ALL_AVAIL = 1,
		UNSPECIFIED = 2
	};

public:
	uintList(PyObject * obj = Py_True);

	/// CPPONLY
	uintList(const vectoru & values) : m_elems(values), m_status(REGULAR)
	{
	}


	/// CPPONLY
	bool allAvail() const
	{
		return m_status == ALL_AVAIL;
	}


	bool unspecified() const
	{
		return m_status == UNSPECIFIED;
	}


	/// CPPONLY
	const vectoru & elems() const
	{
		return m_elems;
	}


protected:
	vectoru m_elems;

private:
	listStatus m_status;
};


// At the CPP level, lociList() is ALL_AVAIL, lociList(vectoru()) is regular empty,
// lociList(NULL) is UNSPECIFIED
// At the Python level True is ALL_AVAIL, a list if regular, False is UNSPECIFIED
class GenoStruTrait;
class Population;

class lociList
{
public:
	typedef vectoru::iterator iterator;
	typedef vectoru::const_iterator const_iterator;

private:
	enum listStatus {
		REGULAR = 0,
		ALL_AVAIL = 1,
		UNSPECIFIED = 2,
		FROM_NAME = 3,
		FROM_POSITION = 4,
		FROM_FUNC = 5,
	};

public:
	lociList(PyObject * obj = Py_True);

	/// CPPONLY
	lociList(const vectoru & values) : m_elems(values), m_names(), m_func(NULL), m_status(REGULAR), m_trait(MaxTraitIndex)
	{
	}


	bool empty() const
	{
		return m_status != ALL_AVAIL && m_elems.empty();
	}


	/// CPPONLY
	size_t size() const
	{
		return m_status == ALL_AVAIL ? 0 : m_elems.size();
	}


	/// CPPONLY
	string name(size_t i) const
	{
		DBG_FAILIF(i >= m_names.size(), ValueError, "Index out of range.");
		return m_names[i];
	}


	/// CPPONLY
	bool allAvail() const
	{
		return m_status == ALL_AVAIL;
	}


	/// CPPONLY
	bool unspecified() const
	{
		return m_status == UNSPECIFIED;
	}


	bool dynamic() const
	{
		return m_status == FROM_NAME;
	}


	/// CPPONLY
	/// return the index of loc in the index list
	size_t indexOf(size_t loc) const;

	/// CPPONLY
	const vectoru & elems(const GenoStruTrait * trait = NULL) const;

	/// CPPONLY
	const vectoru & elems(const Population * trait) const;

protected:
	mutable vectoru m_elems;
	vectorstr m_names;
	vectorpos m_positions;

	mutable pyFunc m_func;
	// to control m_func so that it is called only once for
	// each generation
	mutable size_t m_func_gen;

private:
	listStatus m_status;


	mutable TraitIndexType m_trait;
	// used by indexOf
	mutable std::map<size_t, size_t> m_lociMap;
};


class floatList
{
public:
	floatList(PyObject * obj = NULL);

	/// CPPONLY
	floatList(double val) : m_elems(1, val)
	{
	}


	/// CPPONLY
	floatList(const vectorf & val) : m_elems(val)
	{
	}


	/// CPPONLY
	const vectorf & elems() const
	{
		return m_elems;
	}


protected:
	vectorf m_elems;
};

/* stringList() is supposed to work just like intList and I should be
 * able to define two constructors, using string and vectorstr. However,
 * because string itself is a sequence, the vectorstr version accepts
 * input such as 'abs' and yield ('a', 'b', 'c'). I therefore have to
 * parse the input Python object by myself.
 * NOTE
 * - stringList() lead to m_allAvail = true
 * - stringList(vectorstr()) lead to m_allAvail = false
 * - stringList(str) lead to m_allAvail = false
 * This is why most operators accept either vectorstr() or stringList(str).
 */
class stringList
{
public:
	stringList(PyObject * str = NULL);

	/// CPPONLY
	stringList(const string & str) : m_elems(1, str), m_allAvail(false), m_trait(MaxTraitIndex)
	{
	}


	/// CPPONLY
	stringList(const string & str1, const string & str2) : m_elems(), m_allAvail(false), m_trait(MaxTraitIndex)
	{
		m_elems.push_back(str1);
		m_elems.push_back(str2);
	}


	/// CPPONLY
	stringList(const vectorstr & str) : m_elems(str), m_allAvail(false), m_trait(MaxTraitIndex)
	{
	}


	/// CPPONLY
	void obtainFrom(const stringList & items, const char * allowedItems[],
		const char * defaultItems[]);

	/// CPPONLY
	bool allAvail() const
	{
		return m_allAvail;
	}


	/// CPPONLY
	bool empty() const
	{
		return m_elems.empty();
	}


	/// CPPONLY
	bool contains(const string & str) const
	{
		return find(m_elems.begin(), m_elems.end(), str) != m_elems.end();
	}


	/// CPPONLY
	void push_back(const string & str)
	{
		m_elems.push_back(str);
	}


	/// CPPONLY
	/// Return lists. Population from information fields of a population
	/// if allAvail is true.
	const vectorstr & elems(const GenoStruTrait * trait = NULL) const;


private:
	void addString(PyObject * str);

protected:
	mutable vectorstr m_elems;
	bool m_allAvail;
	mutable TraitIndexType m_trait;
};


class intMatrix
{
public:
	intMatrix(PyObject * obj = NULL);

	/// CPPONLY
	bool empty() const
	{
		return m_elems.empty();
	}


	/// CPPONLY
	const matrixi & elems() const
	{
		return m_elems;
	}


protected:
	matrixi m_elems;
};


class floatMatrix
{
public:
	floatMatrix(PyObject * obj = NULL);

	/// CPPONLY
	bool empty() const
	{
		return m_elems.empty();
	}


	/// CPPONLY
	const matrixf & elems() const
	{
		return m_elems;
	}


protected:
	matrixf m_elems;
};


class stringMatrix
{
public:
	stringMatrix(PyObject * str = NULL);

	/// CPPONLY
	bool empty() const
	{
		return m_elems.empty();
	}


	/// CPPONLY
	const matrixstr & elems() const
	{
		return m_elems;
	}


protected:
	matrixstr m_elems;
};


class uintString
{
public:
	uintString(size_t value) : m_string(), m_int(value)
	{
	}


	uintString(const string & name) : m_string(name), m_int(0)
	{
	}


	/// CPPONLY
	bool empty() const
	{
		return m_string.empty();
	}


	/// CPPONLY
	const string & name() const
	{
		return m_string;
	}


	/// CPPONLY
	size_t value() const
	{
		return m_int;
	}


private:
	string m_string;
	size_t m_int;
};


/// CPPONLY
string PyObj_AsString(PyObject * str);

class stringFunc
{
public:
	stringFunc(const char * value) :
		m_value(value), m_func(NULL), m_mode()
	{
	}


	// string func is strictly used for BaseOperator and its
	// derived classes. It is therefore OK to make it accept
	// both func and filehandler
	stringFunc(PyObject * obj) : m_value(), m_func(NULL), m_mode()
	{
		PyObject * func = obj;

		// if obj is a WithMode wrapper class
		if (PyObject_HasAttrString(obj, "_with_mode") && PyObject_HasAttrString(obj, "_with_output")) {
			PyObject * mode = PyObject_GetAttrString(obj, "_with_mode");
			m_mode = PyObj_AsString(mode);
			Py_DECREF(mode);
			func = PyObject_GetAttrString(obj, "_with_output");
			Py_DECREF(func);
		}
		if (PyCallable_Check(func))
			m_func = pyFunc(func);
		// is this a file handler (or an object with write function)?
		else if (PyObject_HasAttrString(func, "write")) {
			PyObject * write_func = PyObject_GetAttrString(func, "write");
			m_func = pyFunc(write_func);
			Py_DECREF(write_func);
		} else
			throw ValueError("Passed python object is not a function or a file handler");
	}


	/// CPPONLY
	string value() const
	{
		return m_value;
	}


	/// CPPONLY
	pyFunc func() const
	{
		return m_func;
	}


	/// CPPONLY
	bool empty() const
	{
		return m_value.empty() && !m_func.isValid();
	}


	/// COPY
	string mode() const
	{
		return m_mode;
	}


private:
	string m_value;

	pyFunc m_func;

	string m_mode;
};


class uintListFunc : public uintList
{
public:
	uintListFunc(const vectoru & values = vectoru()) :
		uintList(values), m_func(NULL)
	{
	}


	uintListFunc(ULONG value) : uintList(vectoru(1, value)), m_func(NULL)
	{
	}


	uintListFunc(PyObject * func) : uintList(), m_func(func)
	{
	}


	/// CPPONLY
	pyFunc func() const
	{
		return m_func;
	}


	/// CPPONLY
	bool empty() const
	{
		return m_elems.empty();
	}


private:
	pyFunc m_func;
};


class floatListFunc : public floatList
{
public:
	floatListFunc(PyObject * func);

	/// CPPONLY
	floatListFunc(double val) : floatList(val), m_func(NULL)
	{
	}


	/// CPPONLY
	floatListFunc(const vectorf & val) : floatList(val), m_func(NULL)
	{
	}


	/// CPPONLY
	double operator[](size_t i) const
	{
		DBG_FAILIF(i >= size(), IndexError,
			(boost::format("Index %1% out of range of 0 ~ %2%") % i % (size() - 1)).str());
		return m_elems[i];
	}

	/// CPPONLY
	bool valid() const
	{
		return !m_elems.empty() || m_func.func() != NULL;
	}

	/// CPPONLY
	bool empty() const
	{
		return m_elems.empty();
	}


	/// CPPONLY
	size_t size() const
	{
		return m_elems.size();
	}


	/// CPPONLY
	pyFunc func() const
	{
		return m_func;
	}


private:
	pyFunc m_func;
};


// ////////////////////////////////////////////////////////////
// / Shared variables
// ////////////////////////////////////////////////////////////

/// CPPONLY
void PyObj_As_Bool(PyObject * obj, bool & val);

/// CPPONLY
void PyObj_As_Int(PyObject * obj, long & val);

/// CPPONLY
void PyObj_As_SizeT(PyObject * obj, size_t & val);

/// CPPONLY
void PyObj_As_Double(PyObject * obj, double & val);

/// CPPONLY
void PyObj_As_String(PyObject * obj, string & val);

/// CPPONLY
void PyObj_As_Array(PyObject * obj, vectorf & val);

/// CPPONLY
void PyObj_As_IntArray(PyObject * obj, vectori & val);

/// CPPONLY
void PyObj_As_SizeTArray(PyObject * obj, vectoru & val);

/// CPPONLY
PyObject * Allele_Vec_As_NumArray(GenoIterator begin, GenoIterator end);

/// CPPONLY
PyObject * Lineage_Vec_As_NumArray(LineageIterator begin, LineageIterator end);

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
	PyObject * getVar(const string & name, bool nameError = true) const;

	bool hasVar(const string & name)
	{
		// not null = has variable
		return PyDict_GetItemString(m_dict, const_cast<char *>(name.c_str())) != NULL;
	}


	/// remove variable
	void removeVar(const string & name);

	///CPPONLY
	PyObject * setVar(const string & name, const bool val);

	///CPPONLY
	PyObject * setVar(const string & name, const long val);

	///CPPONLY
	PyObject * setVar(const string & name, const size_t val);

	///CPPONLY
	PyObject * setVar(const string & name, const double val);

	///CPPONLY
	PyObject * setVar(const string & name, const string & val);

	///CPPONLY
	PyObject * setVar(const string & name, const vectori & val);

	///CPPONLY
	PyObject * setVar(const string & name, const vectoru & val);

	///CPPONLY
	PyObject * setVar(const string & name, const vectorf & val);

	///CPPONLY
	PyObject * setVar(const string & name, const strDict & val);

	///CPPONLY
	PyObject * setVar(const string & name, const intDict & val);

	///CPPONLY
	PyObject * setVar(const string & name, const uintDict & val);

	///CPPONLY
	PyObject * setVar(const string & name, const tupleDict & val);

	/// CPPONLY
	bool getVarAsBool(const string & name, bool nameError = true) const
	{
		bool val;

		PyObj_As_Bool(getVar(name, nameError), val);
		return val;
	}


	/// CPPONLY
	long getVarAsInt(const string & name, bool nameError = true) const
	{
		long val;

		PyObj_As_Int(getVar(name, nameError), val);
		return val;
	}


	/// CPPONLY
	double getVarAsDouble(const string & name, bool nameError = true) const
	{
		double val;

		PyObj_As_Double(getVar(name, nameError), val);
		return val;
	}


	/// CPPONLY, obsolete
	// void numArrayOwnMem(const string& name);

	/// CPPONLY
	string getVarAsString(const string & name, bool nameError = true) const
	{
		string val;

		PyObj_As_String(getVar(name, nameError), val);
		return val;
	}


	/// CPPONLY
	void getVarAsIntDict(const string & name, uintDict & res, bool nameError = true) const;

	/// CPPONLY
	void getVectorVarAsIntDict(const string & name, uintDict & res, bool nameError = true) const;

	PyObject * & dict()
	{
		return m_dict;
	}


	/// CPPONLY, save m_dist as string
	string asString() const;

	/// CPPONLY
	void fromString(const string & vars);

	/// CPPONLY, save m_dist to pickle
	string to_pickle() const;

	/// CPPONLY
	void from_pickle(const string & vars);

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

/// CPPONLY
void * pyIndPointer(PyObject * p);

/// CPPONLY
void * pyPopPointer(PyObject * p);

/// CPPONLY
void * pyOpPointer(PyObject * p);

/// CPPONLY
string shorten(const string & val, size_t length = 40);

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
	Expression(const string & expr = string(), const string & stmts = string(),
		PyObject * locals = NULL)
		: m_expr(NULL), m_stmts(NULL), m_locals(locals)
	{
		if (m_locals == NULL)
			m_locals = mainVars().dict();

		/* PyEval_GetBuildtins cannot be executed without a Python interpreter
		 * so it is commented out when simuPOP is built as a standalone executable.
		 */
		if (PyDict_GetItemString(m_locals, "__builtins__") == NULL)
			if (PyDict_SetItemString(m_locals, "__builtins__", PyEval_GetBuiltins()) != 0)
				throw RuntimeError("Cannot set __builtins__ for a dictionary.");

		// empty expression
		if (expr.empty() && stmts.empty())
			return;

		//detect leading spaces from python expressions
		DBG_FAILIF(!expr.empty() && (expr[0] == ' ' || expr[0] == '\t'), ValueError,
			"Can not include leading space in python expression '" + expr + "'");
		DBG_FAILIF(!stmts.empty() && (stmts[0] == ' ' || stmts[0] == '\t'), ValueError,
			"Can not include leading space in python statement '" + stmts + "'");

		compileExpr(expr);
		compileStmts(stmts);
	}


	~Expression();

	/// Copy constructor, need to be defined because of ref count issue.
	/// CPPONLY
	Expression(const Expression & rhs);

	/// CPPONLY
	string expr() const
	{
		return m_exprString;
	}


	/// CPPONLY
	string stmts() const
	{
		return m_stmtsString;
	}


	/** CPPONLY
	 * set local dictionary
	 */
	void setLocalDict(PyObject * dict) const
	{
		m_locals = dict;
		if (PyDict_GetItemString(m_locals, "__builtins__") == NULL)
			if (PyDict_SetItemString(m_locals, "__builtins__", PyEval_GetBuiltins()) != 0)
				throw RuntimeError("Cannot set __builtins__ for a dictionary.");
	}


	/// CPPONLY
	bool empty()
	{
		return m_expr == NULL && m_stmts == NULL;
	}


	/// CPPONLY
	void setExpr(const string & expr = string())
	{
		compileExpr(expr);
	}


	/// CPPONLY
	void setStmts(const string & stmts = string())
	{
		compileStmts(stmts);
	}


	/// CPPONLY evaluate with PyObject* output
	PyObject * evaluate() const;

	/// CPPONLY  return bool value
	bool valueAsBool() const;

	/// CPPONLY  return dictionary value
	long int valueAsInt() const;

	/// CPPONLY  return double value
	double valueAsDouble() const;

	/// CPPONLY  return string value
	string valueAsString() const;

	/// CPPONLY  return array value
	vectorf valueAsArray() const;

private:
	/// compile expression into byte code
	void compileExpr(const string & expr);

	/// compile statements into byte code
	void compileStmts(const string & stmts);

	string m_exprString;
	string m_stmtsString;

	/// compiled code
	PyObject * m_expr, * m_stmts;

	/// local namespace
	mutable PyObject * m_locals;
};


/** CPPONLY
 *  Define a class that decipher simple statements that
 *  can be executed faster for InfoExec. It does not
 *  execute the statement.
 */
class simpleStmt
{
public:
	enum OperationType {
		NoOperation = 100,
		Assignment = 101,
		Increment = 102,
		Decrement = 103,
		MultipliedBy = 104,
		SetSex = 105,
		SetAffected = 106,
		SetUnaffected = 107,
	};

public:
	simpleStmt(const string & stmt, const string & indVar);

	string var() const
	{
		return m_var;
	}


	OperationType operation() const
	{
		return m_operation;
	}


	double value() const
	{
		return m_value;
	}


private:
	string m_var;

	OperationType m_operation;

	double m_value;
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

	/// CPPONLY
	void closeOstream(const string & filename);

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
   StreamProvider.getStream()
   StreamProdiver.closeStream()
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
	StreamProvider(const string & output, const pyFunc & func, const string & mode);

	~StreamProvider()
	{
	}


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

	/// output to a function
	pyFunc m_func;

	/// internal m_flags of the operator. They are set during initialization for
	/// performance considerations.
	static const unsigned char m_flagNoOutput = 1;
	static const unsigned char m_flagUseDefault = 2;
	static const unsigned char m_flagAppend = 4;
	static const unsigned char m_flagRealAppend = 8;
	static const unsigned char m_flagCloseAfterUse = 16;
	static const unsigned char m_flagUseString = 32;
	static const unsigned char m_flagReadable = 64;
	static const unsigned char m_flagUseFunc = 128;

	/// m_flags
	unsigned char m_flags;

	string m_mode;

	/// file stream pointer, used ONLY if we create a new file stream
	ostream * m_filePtr;
};


/** Output files specified by \c '>' are closed immediately after they are
 *  written. Those specified by \c '>>' and \c '>>>' are closed by a
 *  simulator after <tt>Simulator.evolve()</tt>. However, these files will
 *  be kept open if the operators are applied directly to a population using
 *  the operators' function form. In this case, function \c closeOutput can be
 *  used to close a specific file \e output, and close all unclosed files if
 *  \e output is unspecified. An exception will be raised if \e output does
 *  not exist or it has already been closed.
 */
void closeOutput(const string & output = string());

// ////////////////////////////////////////////////////////////
// / Random number generator
// ////////////////////////////////////////////////////////////


/// CPPONLY
class RNG_func
{
public:
	RNG_func(gsl_rng * rng) : m_RNG(rng)
	{

	}


	unsigned long int operator()(unsigned long int N) const
	{
		return gsl_rng_uniform_int(m_RNG, N);
	}


private:
	gsl_rng * m_RNG;
};


/** This random number generator class wraps around a number of random number
 *  generators from GNU Scientific Library. You can obtain and change the
 *  RNG used by the current simuPOP module through the \c getRNG() function,
 *  or create a separate random number generator and use it in your script.
 */
class RNG
{

public:
	/** Create a RNG object using specified name and seed. If \e rng is not
	 *  given, environmental variable \c GSL_RNG_TYPE will be used if it is
	 *  available. Otherwise, generator \c mt19937 will be used. If \e seed is
	 *  not given, <tt>/dev/urandom</tt>, <tt>/dev/random</tt>, or other system
	 *  random number source will be used to guarantee that random seeds are
	 *  used even if more than one simuPOP sessions are started simultaneously.
	 *  Names of supported random number generators are available from
	 *  <tt>moduleInfo()['availableRNGs']</tt>.
	 */
	RNG(const char * name = NULL, unsigned long seed = 0);

	/// CPPONLY Copy contructor, needed because of m_RNG
	RNG(const RNG &);

	///
	~RNG();

	/** Replace the existing random number generator using RNG \e name with
	 *  seed \e seed. If \e seed is 0, a random seed will be used. If \e name
	 *  is empty, use the existing RNG but reset the seed.
	 *  <group>1-setup</group>
	 */
	void set(const char * name = NULL, unsigned long seed = 0);


	/** Return the name of the current random number generator.
	 *  <group>2-info</group>
	 */
	const char * name() const
	{
		return gsl_rng_name(m_RNG);
	}


	/** Return the seed used to initialize the RNG. This can be used to
	 *  repeat a previous session.
	 *  <group>2-info</group>
	 */
	unsigned long seed() const
	{
		return m_seed;
	}


	/// CPPONLY
	static unsigned long generateRandomSeed();


	/** Generate a random number following a rng_uniform [0, 1) distribution.
	 *  <group>3-rng</group>
	 */
	double randUniform()
	{
		return gsl_rng_uniform(m_RNG);
	}


	/** Return a random bit. This is not part of GSL.
	 *  HIDDEN
	 */
	bool randBit();

	/** return a random number in the range of <tt>[0, 1, 2, ... n-1]</tt>
	 *  <group>3-rng</group>
	 */
	unsigned long int randInt(unsigned long int n)
	{
		return gsl_rng_uniform_int(m_RNG, n);
	}


	/** Generate a random number following a normal distribution with mean
	 *  \e mu and standard deviation \e sigma.
	 *  <group>4-distribution</group>
	 */
	double randNormal(double mu, double sigma)
	{
		return gsl_ran_gaussian(m_RNG, sigma) + mu;
	}


	/** Generate a random number following a exponential distribution with
	 *  parameter \e mu.
	 *  <group>4-distribution</group>
	 */
	double randExponential(double mu)
	{
		return gsl_ran_exponential(m_RNG, mu);
	}


	/** Generate a random number following a gamma distribution with
	 *  a shape parameters \e a and scale parameter \e b.
	 *  <group>4-distribution</group>
	 */
	double randGamma(double a, double b)
	{
		return gsl_ran_gamma(m_RNG, a, b);
	}


	/** Generate a random number following a Chi-squared distribution with
	 *  \e nu degrees of freedom.
	 *  <group>4-distribution</group>
	 */
	double randChisq(double nu)
	{
		return gsl_ran_chisq(m_RNG, nu);
	}


	/** Generate a random number following a geometric distribution with
	 *  parameter \e p.
	 *  <group>4-distribution</group>
	 */
	long randGeometric(double p)
	{
		return gsl_ran_geometric(m_RNG, p);
	}


	/** Generate a random number following a binomial distribution with
	 *  parameters \e n and \e p.
	 *  <group>4-distribution</group>
	 */
	ULONG randBinomial(UINT n, double p)
	{
		DBG_FAILIF(n <= 0, ValueError, "RandBinomial: n should be positive.");

		return gsl_ran_binomial(m_RNG, p, n);
	}


	/** Generate a random number following a Poisson distribution with
	 *  parameter \e mu.
	 *  <group>4-distribution</group>
	 */
	ULONG randPoisson(double mu)
	{
		return gsl_ran_poisson(m_RNG, mu);
	}


	/** Generate a positive random number following a zero-truncated Poisson
	 *  distribution with parameter \e mu.
	 */
	ULONG randTruncatedPoisson(double mu);

	/** Generate a positive random number following a zero-truncated binomial
	 *  distribution with parameters \e n and \e p.
	 */
	ULONG randTruncatedBinomial(UINT n, double p);

	/** Generate a random number following a multinomial distribution with
	 *  parameters \e N and \e p (a list of probabilities).
	 *  <group>4-distribution</group>
	 */
	vectoru randMultinomial(unsigned int N, const vectorf & p)
	{
		// if sum p_i != 1, it will be normalized.
		// the size of n is not checked!
		vector<unsigned int> val(p.size());
		vectoru res(p.size());
		gsl_ran_multinomial(m_RNG, p.size(), N, &p[0], &val[0]);
		for (size_t i = 0; i < p.size(); ++i)
			res[i] = val[i];
		return res;
	}


	/** Randomly shuffle a sequence
	 *  CPPONLY
	 */
	template<typename T>
	void randomShuffle(T begin, T end) const
	{
		// use the STL random shuffle function but a GSL random number
		// generator.
		RNG_func rng(m_RNG);

		std::random_shuffle(begin, end, rng);
	}


private:
	ULONG search_poisson(UINT y, double * z, double p, double lambda);

	ULONG search_binomial(UINT y, double * z, double p, UINT n, double pr);

private:
	/// global random number generator
	gsl_rng * m_RNG;

	/// seed used
	unsigned long m_seed;

	/// used by RNG::rand_bit(). I was using static but this make it difficult
	/// to reset a RNG when a new seed is set.
	uint16_t m_bitByte;
	UINT m_bitIndex;
};

/// return the currently used random number generator
RNG & getRNG();

/// CPPONLY
void chisqTest(const vector<vectoru> & table, double & chisq, double & chisq_p);

/// CPPONLY
double armitageTrendTest(const vector<vectoru> & table, const vectorf & weight);

/// CPPONLY
double hweTest(const vectoru & cnt);

/// CPPONLY
template <typename IT>
void propToCount(IT first, IT last, size_t N, vectoru & count)
{
	size_t sz = last - first;

	count.resize(sz);
	size_t tot = 0;
	for (size_t i = 0; i < sz; ++i) {
		count[i] = static_cast<ULONG>(N * *(first + i) + 0.5);
		tot += count[i];
		if (tot > N) {
			count[i] -= tot - N;
			for (size_t j = i + 1; j < sz; ++j)
				count[j] = 0;
			return;
		}
	}
	if (N == tot)
		return;
	// if tot < N, spead the round offs to the first several counts
	for (size_t i = 0; tot < N && i < sz; ++i) {
		if (count[i] < *(first + i) * N) {
			count[i] += 1;
			++tot;
		}
	}
	if (N != tot)
		count.back() += N - tot;

}


/// CPPONLY
string formatDescription(const string & text);

/** A random number generator that returns \c 0, \c 1, ..., \c k-1 with
 *  probabilites that are proportional to their weights. For example, a
 *  weighted sampler with weights \c 4, \c 3, \c 2 and \c 1 will return numbers
 *  \c 0, \c 1, \c 2 and \c 3 with probabilities \c 0.4, \c 0.3, \c 0.2 and
 *  \c 0.1, respectively. If an additional parameter \c N is specified, the
 *  weighted sampler will return exact proportions of numbers if \c N numbers
 *  are returned. The version without additional parameter is similar to
 *  the <tt>sample(prob, replace=FALSE)</tt> function of the R statistical
 *  package.
 */
class WeightedSampler
{
public:
	/** Creates a weighted sampler that returns \c 0, \c 1, ... \c k-1 when a
	 *  list of \c k weights are specified (\e weights). \e weights do not have
	 *  to add up to 1. If a non-zero \e N is specified, exact proportions of
	 *  numbers will be returned in \e N returned numbers.
	 */
	WeightedSampler(const vectorf & weights = vectorf(), ULONG N = 0)
		: m_algorithm(0), m_q(0), m_a(0), m_param(0),
		m_sequence(0), m_index(0)
	{

		set(weights.begin(), weights.end(), N);
	}


	/// destructor
	~WeightedSampler()
	{
	}


	/** CPPONLY
	 *  Set parameters for the weighted sampler.
	 */
	template<typename IT>
	void set(IT first, IT last, size_t N = 0)
	{
		size_t sz = last - first;

		// this is the case with unknown number of outputs
		if (N == 0) {
			m_N = sz;
			// no weight (wrong case)
			if (m_N == 0) {
				// invalid
				m_algorithm = 0;
				return;
			}
			// only one weight?
			if (m_N == 1) {
				// return 0 all time
				m_algorithm = 1;
				m_param = 0;
				return;
			}
			// fixed value
			bool allEqual = true;

			for (size_t i = 1; i < sz; ++i) {
				if (*(first + i) != *(first + i - 1)) {
					allEqual = false;
					break;
				}
			}

			if (allEqual) {
				m_algorithm = 2;
				m_param = m_N;
				return;
			}
			// only one value
			bool fixed = true;
			ssize_t prevIndex = -1;
			for (size_t i = 0; i < sz; ++i) {
				if (*(first + i) != 0) {
					if (prevIndex == -1) {
						m_param = i;
						prevIndex = i;
					} else { // two non-zero index, not fixed.
						fixed = false;
						break;
					}
				}
			}
			if (fixed) {
				m_algorithm = 1;
				return;
			}
			// the mos difficult case
			m_algorithm = 3;
			// sum of weight
			double w = accumulate(first, last, 0.0);

			DBG_FAILIF(fcmp_le(w, 0), ValueError, "Sum of weight is <= 0.");

			w = m_N / w;

			// initialize p with N*p0,...N*p_k-1
			m_q.resize(m_N);

			for (size_t i = 0; i < m_N; ++i)
				m_q[i] = *(first + i) * w;

			// initialize Y with values
			m_a.resize(m_N);
			for (size_t i = 0; i < m_N; ++i)
				m_a[i] = i;
			// use two sets H and L
			// for efficiency purpose, use a single vector.
			size_t * HL = new size_t[m_N];
			size_t * L = HL;
			size_t * H = HL + m_N - 1;                                 // point to the end.

			for (size_t i = 0; i < m_N; ++i) {
				if (m_q[i] > 1)
					*H-- = i;
				else
					*L++ = i;
			}

			//
			size_t j, k;
			while (L != HL && H != HL + m_N - 1) {
				j = *(L - 1);
				k = *(H + 1);
				m_a[j] = k;
				m_q[k] += m_q[j] - 1;

				L--;                                                                    // remove j from L
				if (m_q[k] < 1.) {
					*L++ = k;                                                           // add k to L
					++H;                                                                // remove k from H
				}
			}
			delete[] HL;
		} else {
			m_algorithm = 4;
			for (size_t i = 0; i < sz; ++i) {
				DBG_FAILIF(*(first + i) < 0 || *(first + i) > 1, ValueError,
					"Proportions should be between 0 and 1");
			}
			// sum of weight
			double w = accumulate(first, last, 0.0);
			(void)w; // fix compiler warning.

			DBG_FAILIF(fcmp_eq(w, 0), ValueError, "Proportions sum up to 0");

			vectoru count(N);
			propToCount(first, last, N, count);

			m_sequence.resize(N);
			// turn weight into percentage
			for (size_t i = 0, j = 0; i < sz; ++i)
				for (size_t k = 0; k < count[i]; ++k, ++j)
					m_sequence[j] = i;

			// random shuffle
			getRNG().randomShuffle(m_sequence.begin(), m_sequence.end());
			m_index = 0;
		}

	}


	/** Returns a random number between \c 0 and \c k-1 with probabilities that
	 *  are proportional to specified weights.
	 */
	size_t draw();

	/** Returns a list of \e n random numbers
	 */
	vectoru drawSamples(ULONG n = 1);

private:
	/// which algorithm to use
	int m_algorithm;

	/// length of weight.
	size_t m_N;

	/// internal table.
	vectorf m_q;

	///
	vectoru m_a;

	///
	size_t m_param;

	///
	vectoru m_sequence;

	ATOMICLONG m_index;
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
 *  using fast bernulliTrial method for fix p,
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
class Bernullitrials
{
public:
	/// CPPONLY
	Bernullitrials(RNG & /* rng */);

	///
	Bernullitrials(RNG & /* rng */, const vectorf & prob, ULONG trials = 0);

	///
	~Bernullitrials();

	/** CPPONLY
	 * return size of trial
	 */
	size_t trialSize() const
	{
		return m_N;
	}


	size_t probSize() const
	{
		return m_prob.size();
	}


	/// CPPONLY
	void setParameter(const vectorf & prob, size_t trials = 0);

	/// generate the trial table, reset m_cur
	void doTrial();

	/// CPPONLY
	size_t curTrial();

	/// if necessary, do trail again.
	void trial();

#define getBit(ptr, i)    ((*((ptr) + (i) / WORDBIT) & (1UL << ((i) - ((i) / WORDBIT) * WORDBIT))) != 0)
	inline bool trialSucc(size_t idx) const
	{
		DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
		return getBit(m_pointer[idx], m_cur);
	}


	inline bool trialSucc(size_t idx, size_t cur) const
	{
		return getBit(m_pointer[idx], cur);
	}


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
	// We cannot cache m_RNG because differenct m_RNG will be used for
	// different threads
	//RNG * m_RNG;

	/// number of trials.
	size_t m_N;

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


class Bernullitrials_T
{
public:
	/// CPPONLY
	Bernullitrials_T(RNG & /* rng */);

	///
	Bernullitrials_T(RNG & /* rng */, const vectorf & prob, size_t N = 1024);

	///
	~Bernullitrials_T();

	/// CPPONLY
	void setParameter(const vectorf & prob, size_t N = 1024);

	/// generate the trial table, reset m_cur
	void doTrial();

	/// if necessary, do trail again.
	void trial();


#define getBit(ptr, i)    ((*((ptr) + (i) / WORDBIT) & (1UL << ((i) - ((i) / WORDBIT) * WORDBIT))) != 0)
	inline bool trialSucc(size_t idx) const
	{
		// DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
		return getBit(m_pointer[m_cur], idx);
	}


	/*
	    inline bool trialSucc(size_t idx, size_t cur) const
	    {
	        return getBit(m_pointer[cur], idx);
	    }
	 */

	// first and next succ across prob
	size_t probFirstSucc() const;

	size_t probNextSucc(size_t pos) const;

	void setTrialSucc(size_t idx, bool succ);

	/// return the succ rate for one index, used for verification pruposes
	double trialSuccRate(UINT index) const;

	/// return the succ rate for current trial, used for verification pruposes
	double probSuccRate() const;

	/*
	   /// CPPONLY
	   vectorf probabilities()
	   {
	    return m_prob;
	   }
	 */

public:
	static const size_t npos = static_cast<size_t>(-1);

private:
	void setAll(size_t idx, bool v);

private:
	// We cannot cache m_RNG because differenct m_RNG will be used for
	// different threads
	//RNG * m_RNG;

	size_t m_N;

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

	/// than using the reference interface.
	vector<WORDTYPE *> m_pointer;

	/// current trial. Used when user want to access the table row by row
	size_t m_cur;
};


// ////////////////////////////////////////////////////////////
// /  Global debug and initialization related functions
// ////////////////////////////////////////////////////////////

/** Return a dictionary with information regarding the currently loaded simuPOP
 *  module. This dictionary has the following keys:
 *  \li \c revision: revision number.
 *  \li \c version: simuPOP version string.
 *  \li \c optimized: Is this module optimized (\c True or \c False).
 *  \li \c alleleType: Allele type of the module (\c short, \c long or \c binary).
 *  \li \c maxAllele: the maximum allowed allele state, which is \c 1 for
 *       binary modules, \c 255 for short modules and \c 65535 for long modules.
 *  \li \c compiler: the compiler that compiles this module.
 *  \li \c date: date on which this module is compiled.
 *  \li \c python: version of python.
 *  \li \c platform: platform of the module.
 *  \li \c wordsize: size of word, can be either 32 or 64.
 *  \li \c alleleBits: the number of bits used to store an allele
 *  \li \c maxNumSubPop: maximum number of subpopulations.
 *  \li \c maxIndex: maximum index size (limits population size * total number of marker).
 *  \li \c debug: A dictionary with debugging codes as keys and the status of each
 *       debugging code (\c True or \c False) as their values.
 */
PyObject * moduleInfo();

#ifdef BINARYALLELE
// efficiently copy alleles (block by block, rather than 1 by 1)
/// CPPONLY
void copyGenotype(GenoIterator fr, GenoIterator to, size_t n);

/// CPPONLY
/// This function is not called fillGenotype because when A is not zero
/// fillGenotype should fill all genotypes with mutant, something we
/// will not do here.
void clearGenotype(GenoIterator to, size_t n);

#  ifndef OPTIMIZED
void testCopyGenotype();

#  endif
#endif

#ifdef MUTANTALLELE
/// COPY MUTANT ALLELES
/// CPPONLY
void copyGenotype(const ConstGenoIterator begin, const ConstGenoIterator end,
	GenoIterator it);

/// CPPONLY
void clearGenotype(GenoIterator begin, GenoIterator end);

#endif


/// CPPONLY initialize module simuPOP when using "import simuPOP"
bool initialize(PyObject * module);

/// CPPONLY get a null stream that discard everything
ostream & cnull();

}
#endif
