/**
 *  $File: utility.h $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
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

/** Set debug code \e code. Name of available codes are available from \c DebugCodes.
 */
void TurnOnDebug(DBG_CODE code);

/** Turn off debug code \e code. Default to turn off all debug codes.
 */
void TurnOffDebug(DBG_CODE code = DBG_ALL);

#ifndef OPTIMIZED
/// test if one code is turned on, CPPONLY
bool debug(DBG_CODE code);

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
	T operator()(void converter(PyObject *, T &), const char * format, ...)
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


	PyObject * operator()(const char * format, ...)
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


private:
	pyObject m_func;
};

/** This class defines an interface using which both a integer number and
 *  a list of numbers can be accpted.
 */
class intList
{
public:
	intList(const vectorl & values = vectorl()) : m_elems(values)
	{
	}


	intList(long value) : m_elems(1, value)
	{
	}


	/// CPPONLY
	const vectorl & elems() const
	{
		return m_elems;
	}


protected:
	vectorl m_elems;
};


// I cannot use template here because otherwise SWIG does not
// handle the type correctly. I guess this can be my problem with using
// of template in simuPOP_common.i

class uintList
{
public:
	typedef vectoru::iterator iterator;
	typedef vectoru::const_iterator const_iterator;

public:
	uintList(const vectoru & values = vectoru()) : m_elems(values)
	{
	}


	uintList(ULONG value) : m_elems(1, value)
	{
	}


	/// CPPONLY
	const vectoru & elems() const
	{
		return m_elems;
	}


protected:
	vectoru m_elems;
};


class lociList
{
public:
	lociList(PyObject * obj = NULL);

	/// CPPONLY
	lociList(const vectoru & values) : m_elems(values), m_allAvail(false)
	{
	}


	/// CPPONLY
	bool allAvail() const
	{
		return m_allAvail;
	}


	/// CPPONLY
	const vectoru & elems() const
	{
		return m_elems;
	}


private:
	vectoru m_elems;
	bool m_allAvail;
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
	stringList(const string & str) : m_elems(1, str), m_allAvail(false)
	{
	}


	/// CPPONLY
	stringList(const string & str1, const string & str2) : m_elems(), m_allAvail(false)
	{
		m_elems.push_back(str1);
		m_elems.push_back(str2);
	}


	/// CPPONLY
	stringList(const vectorstr & str) : m_elems(str), m_allAvail(false)
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
	const vectorstr & elems() const
	{
		return m_elems;
	}


private:
	void addString(PyObject * str);

protected:
	vectorstr m_elems;
	bool m_allAvail;
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
	uintString(UINT value) : m_string(), m_int(value)
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
	UINT value() const
	{
		return m_int;
	}


private:
	string m_string;
	UINT m_int;
};


class stringFunc
{
public:
	stringFunc(const char * value) :
		m_value(value), m_func(NULL)
	{
	}


	stringFunc(PyObject * func) : m_value(), m_func(func)
	{
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


private:
	string m_value;

	pyFunc m_func;
};


class uintListFunc : public uintList
{
public:
	uintListFunc(const vectoru & values = vectoru()) :
		uintList(values), m_func(NULL)
	{
	}


	uintListFunc(ULONG value) : uintList(value), m_func(NULL)
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


	/// CPPONLY
	double operator[](size_t i) const
	{
		DBG_FAILIF(i >= size(), IndexError,
			"Index " + toStr(i) + " out of range of 0 ~ " + toStr(size() - 1));
		return m_elems[i];
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


/** A class to specify replicate list. The reason why I cannot simple
 *  use vectori() is that users have got used to use a single number
 *  to specify a single replicate.
 */
class repList
{
public:
	repList(PyObject * obj = NULL);

	/// CPPONLY
	repList(const vectori & reps) :
		m_elems(reps), m_allAvail(false)
	{
	}


	/// CPPONLY
	vectori elems() const
	{
		return m_elems;
	}


	/// CPPONLY
	bool match(UINT rep, const vector<bool> & activeRep);

private:
	vectori m_elems;
	bool m_allAvail;
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
bool PyObj_Is_IntNumArray(PyObject * obj);

/// CPPONLY
bool PyObj_Is_DoubleNumArray(PyObject * obj);

/// CPPONLY
bool PyObj_Is_AlleleNumArray(PyObject * obj);

/// CPPONLY
PyObject * Double_Vec_As_NumArray(vectorf::iterator begin, vectorf::iterator end);

/// CPPONLY
PyObject * Int_Vec_As_NumArray(vectori::iterator begin, vectori::iterator end);

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

	///CPPONLY
	PyObject * setIntDefDictVar(const string & name, const intDict & val);

	///CPPONLY
	PyObject * setTupleDefDictVar(const string & name, const tupleDict & val);

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

/// CPPONLY
void * pyIndPointer(PyObject * p);

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
	PyObject * evaluate();

	/// CPPONLY  return bool value
	bool valueAsBool();

	/// CPPONLY  return dictionary value
	int valueAsInt();

	/// CPPONLY  return double value
	double valueAsDouble();

	/// CPPONLY  return string value
	string valueAsString();

	/// CPPONLY  return array value
	vectorf valueAsArray();

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


/** CPPONLY
 *  Define a class that decipher simple statements that
 *  can be executed faster for infoExec. It does not
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
	};

public:
	simpleStmt(const string & stmt);

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
	StreamProvider(const string & output, const pyFunc & func);

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
	static const size_t m_flagNoOutput = 1;
	static const size_t m_flagUseDefault = 2;
	static const size_t m_flagAppend = 4;
	static const size_t m_flagRealAppend = 8;
	static const size_t m_flagCloseAfterUse = 16;
	static const size_t m_flagUseString = 32;
	static const size_t m_flagReadable = 64;
	static const size_t m_flagUseFunc = 128;

	/// m_flags
	unsigned char m_flags;

	/// file stream pointer, used ONLY if we create a new file stream
	ostream * m_filePtr;
};


/** Output files specified by \c '>' are closed immediately after they are
 *  written. Those specified by \c '>>' and \c '>>>' are closed by a
 *  simulator after <tt>simulator.evolve()</tt>. However, these files will
 *  be kept open if the operators are applied directly to a population using
 *  the operators' function form. In this case, function \c closeOutput can be
 *  used to close a specific file \e output, and close all unclosed files if
 *  \e output is unspecified. An exception will be raised if \e output does
 *  not exist or it has already been closed.
 */
void CloseOutput(const string & output = string());

// ////////////////////////////////////////////////////////////
// / Random number generator
// ////////////////////////////////////////////////////////////

/** This random number generator class wraps around a number of random number
 *  generators from GNU Scientific Library. You can obtain and change the
 *  RNG used by the current simuPOP module through the \c GetRNG() function,
 *  or create a separate random number generator and use it in your script.
 */
class RNG
{

public:
	/** Create a RNG object using specified name and seed. If \e rng is not
	 *  given, environmental variable \c GSL_RNG_TYPE will be used if it is
	 *  available. Otherwise, RNG \c mt19937 will be used. If \e seed is not
	 *  given, <tt>/dev/urandom</tt>, <tt>/dev/random</tt>, or other system
	 *  random number source will be used to guarantee that random seeds
	 *  are used even if more than one simuPOP sessions are started
	 *  simultaneously.
	 */
	RNG(const char * name = NULL, unsigned long seed = 0);

	///
	~RNG();

	/** Use another underlying RNG for the current RNG object. The handling of
	 *  parameters \e rng and \e seed is the same as \c RNG::RNG(name, seed).
	 */
	void setRNG(const char * name = NULL, unsigned long seed = 0);

	/** Return the name of the current random number generator.
	 */
	const char * name()
	{
		return gsl_rng_name(m_RNG);
	}


	/** Return the seed used to initialize the RNG. This can be used to
	 *  repeat a previous session.
	 */
	unsigned long seed()
	{
		return m_seed;
	}


	/** Return the maximum allowed seed value
	 */
	unsigned long maxSeed()
	{
		return std::numeric_limits<unsigned long>::max();
	}


	/** Set random seed for this random number generator. If seed is 0, method
	 *  described in \c setRNG is used.
	 */
	void setSeed(unsigned long seed)
	{
		m_seed = seed;
		setRNG(name(), m_seed);
	}


	/// CPPONLY
	unsigned long generateRandomSeed();

	/** Maximum value of this RNG
	 */
	unsigned long max()
	{
		return gsl_rng_max(m_RNG);
	}


	string __repr__()
	{
		return toStr("<simuPOP::RNG ") + gsl_rng_name(m_RNG) + ">";
	}


	/** Return a random number in the range of <tt>[0, 2, ... max()-1]</tt>
	 */
	unsigned long int randGet()
	{
		return gsl_rng_get(m_RNG);
	}


	/** Return a random bit.</tt>
	 */
	bool randBit();

	/** return a random number in the range of <tt>[0, 1, 2, ... n-1]</tt>
	 */
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


	/** Generate a random number following a geometric distribution with
	 *  parameter \e p. Please check the documentation of \c gsl_ran_geometric
	 *  for details.
	 */
	int randGeometric(double p)
	{
		return gsl_ran_geometric(m_RNG, p);
	}


	/** Generate a random number following a uniform distribution between 0 and
	 *  1. Please check the documentation of \c gsl_ran_uniform for details.
	 */
	double randUniform01()
	{
		return gsl_rng_uniform(m_RNG);
	}


	/** Generate a random number following a normal distribution with mean
	 *  \e m and standard deviation \e v. Please check the documentation
	 *  of \c gsl_ran_gaussian for details.
	 */
	double randNormal(double m, double v)
	{
		return gsl_ran_gaussian(m_RNG, v) + m;
	}


	/** Generate a random number following a exponential distribution with
	 *  parameter \e v. Please check the documentation of
	 *  \c gsl_ran_exponential for details.
	 */
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


	/** Generate a random number following a binomial distribution with
	 *  parameters \e n and \e p. Please check the documentation of
	 *  \c gsl_ran_binomial for details.
	 */
	UINT randBinomial(UINT n, double p)
	{
		DBG_FAILIF(n <= 0, ValueError, "RandBinomial: n should be positive.");

		return gsl_ran_binomial(m_RNG, p, n);
	}


	/** Generate a random number following a multinomial distribution with
	 *  parameters \e N and \e p (a list of probabilities). Please check the
	 *  documentation of \c gsl_ran_multinomial for details.
	 */
	void randMultinomial(unsigned int N, const vectorf & p, vectoru::iterator n)
	{
		// if sum p_i != 1, it will be normalized.
		// the size of n is not checked!
		vector<unsigned int> val(p.size());
		gsl_ran_multinomial(m_RNG, p.size(), N, &p[0], &val[0]);
		for (size_t i = 0; i < p.size(); ++i)
			*(n++) = val[i];
	}


	/// CPPONLY
	vectoru randMultinomialVal(unsigned int N, const vectorf & p)
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


	/** Generate a random number following a Poisson distribution with
	 *  parameter \e p. Please check the documentation of \c gsl_ran_poisson
	 *  for details.
	 */
	UINT randPoisson(double p)
	{
		return gsl_ran_poisson(m_RNG, p);
	}


private:
	/// global random number generator
	gsl_rng * m_RNG;

	/// seed used
	unsigned long m_seed;
};

/// CPPONLY
void chisqTest(const vector<vectoru> & table, double & chisq, double & chisq_p);

/// CPPONLY
double armitageTrendTest(const vector<vectoru> & table, const vectorf & weight);

/// CPPONLY
double hweTest(const vectoru & cnt);

// weighted sampling using Walker's alias algorithm
class weightedSampler
{
public:
	// set up AliasMethod table
	weightedSampler(RNG & rng, const vectorf & weight = vectorf())
		: m_RNG(&rng), m_q(0), m_a(0), m_fixed(false), m_fixedValue(0)
	{
		set(weight);
	};

	~weightedSampler()
	{
	};

	void set(const vectorf & weight);

	// sample without replacement from 0,...,n-1,
	// with weight freq
	ULONG get()
	{
		if (m_fixed)
			return m_fixedValue;

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
	ULONG get(vectoru & res, ULONG shift = 0)
	{
		double rN;
		size_t K;

		if (m_fixed) {
			std::fill(res.begin(), res.end(), m_fixedValue);
			return 0;
		}
		for (vectoru::iterator it = res.begin(); it != res.end(); ++it) {
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

		if (m_fixed) {
			std::fill(beg, end, m_fixedValue + shift);
			return;
		}
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


	vectoru a()
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

	vectoru m_a;

	// handle special case
	bool m_fixed;
	ULONG m_fixedValue;
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
RNG & GetRNG();

/// set random number generator. If <tt>seed=0</tt> (default), a random seed will be given. If <tt>rng=""</tt>, seed will be set to the current random number generator.
void SetRNG(const string rng = string(), unsigned long seed = 0);

/// list the names of all available random number generators
vectorstr AvailableRNGs();

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
 *  \li \c maxNumSubPop: maximum number of subpopulations.
 *  \li \c maxIndex: maximum index size (limits population size * total number of marker).
 */
PyObject * ModuleInfo();

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

/// CPPONLY get a null stream that discard everything
ostream & cnull();

}
#endif
