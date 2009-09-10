/**
 *  $File: utility.cpp $
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

/**
   \file
   \brief implementation of some utility functions.
 */

#include "utility.h"
#include <cstdlib>
#include "time.h"

#include <bitset>
typedef std::bitset<DBG_CODE_LENGTH> DbgBitSet;

#include <sstream>
using std::stringstream;
using std::ostringstream;
using std::hex;

#include <fstream>
using std::fstream;
using std::ifstream;
using std::ofstream;

// for PySys_WriteStdout and python expressions
#include "swigpyrun.h"

// compile and eval enables compiling string to byte code
#include "compile.h"
#include "eval.h"

// for kbhit
#if  defined (_WIN32) || defined (__WIN32__)
#  include <conio.h>
#else
#  include <termios.h>
#endif

// for CryptGenRandom
#if defined (_WIN32) || defined (__WIN32__)
#  include <windows.h>
#endif

#include "boost/pending/lowest_bit.hpp"
using boost::lowest_bit;

#include "boost/regex.hpp"
using boost::regex;
using boost::regex_match;
using boost::cmatch;

//macros to provide a portable way to make macro-passed string to C++ string
#define MacroQuote_(x) # x
#define MacroQuote(x) MacroQuote_(x)

// these functions are defined in customizedTypes.c which is included
// in simuPOP_wrap.cpp
extern "C" PyObject * newcarrayobject(char * buf, char type, int size);

extern "C" PyObject * newcarrayiterobject(GenoIterator begin, GenoIterator end);

extern "C" bool is_carrayobject(PyObject *);

extern "C" int carray_length(PyObject * a);

extern "C" int carray_itemsize(PyObject * a);

extern "C" char carray_type(PyObject * a);

extern "C" char * carray_data(PyObject * a);

extern "C" PyObject * PyDefDict_New();

extern "C" bool is_defdict(PyTypeObject * type);

extern "C" int initCustomizedTypes(void);

extern "C" PyTypeObject Arraytype;

// for streambuf stuff
#include <streambuf>
using std::streambuf;

#include <numeric>
using std::accumulate;

#include <algorithm>
using std::min_element;
using std::max_element;
using std::find;
using std::sort;
using std::greater;

namespace simuPOP {

// Debug functions
// debug codes in a bitset.
DbgBitSet g_dbgCode;

char * g_debugCodes[] = {
	"DBG_ALL",
	"DBG_GENERAL",
	"DBG_UTILITY",
	"DBG_OPERATOR",
	"DBG_SIMULATOR",
	"DBG_INDIVIDUAL",
	"DBG_OUTPUTER",
	"DBG_MUTATOR",
	"DBG_TRANSMITTER",
	"DBG_INITIALIZER",
	"DBG_POPULATION",
	"DBG_STATOR",
	"DBG_TAGGER",
	"DBG_SELECTOR",
	"DBG_MATING",
	"DBG_MIGRATOR",
	"DBG_PROFILE",
	"DBG_BATCHTESTING",
	"DBG_VISUALIZER",
	"DBG_DEVEL",
	""
};


vectorstr DebugCodes()
{
	vectorstr ret;
	for (size_t i = 0; g_debugCodes[i][0]; ++i)
		ret.push_back(string(g_debugCodes[i]));
	return ret;
}


// set debug area, default to turn all code on
void TurnOnDebug(const string & codeString)
{
#ifndef OPTIMIZED
	if (codeString == "DBG_ALL")
		// set all
		g_dbgCode.set();
	else {
		vectorstr codes;
		size_t lastPos = 0;
		while (codeString.find(",", lastPos) != std::string::npos) {
			size_t pos = codeString.find(",", lastPos);
			codes.push_back(codeString.substr(lastPos, pos));
			lastPos = pos + 1;
		}
		// last piece
		codes.push_back(codeString.substr(lastPos));
		// handle each one
		for (size_t i = 0; i < codes.size(); ++i)
		{
			bool find = false;
			for (size_t c = 0; g_debugCodes[c][0]; ++c) {
				if (string(g_debugCodes[c]) == codes[i]) {
					g_dbgCode[c] = true;
					find = true;
					break;
				}
			}
			if (!find) {
				cerr << "Invalid debug code " << codes[i] << endl;
				exit(1);
			}
		}
	}
#endif
}


// turn off debug, default to turn all code off
void TurnOffDebug(const string & codeString)
{
#ifndef OPTIMIZED
	if (codeString == "DBG_ALL")
		g_dbgCode.reset();
	else {
		vectorstr codes;
		size_t lastPos = 0;
		while (codeString.find(",", lastPos) != std::string::npos) {
			size_t pos = codeString.find(",", lastPos);
			codes.push_back(codeString.substr(lastPos, pos));
			lastPos = pos + 1;
		}
		// last piece
		codes.push_back(codeString.substr(lastPos));
		// handle each one
		for (size_t i = 0; i < codes.size(); ++i)
		{
			bool find = false;
			for (size_t c = 0; g_debugCodes[c][0]; ++c) {
				if (string(g_debugCodes[c]) == codes[i]) {
					g_dbgCode[c] = false;
					find = true;
					break;
				}
			}
			if (!find) {
				cerr << "Invalid debug code " << codes[i] << endl;
				exit(1);
			}
		}
	}
#else
	cerr << "Debug info is ignored in optimized mode." << endl;
#endif
}


#ifndef OPTIMIZED
// test if one code is turned on
// in DEBUG section to make sure it will not be called
// in optimized mode
bool debug(DBG_CODE code)
{
	return g_dbgCode[code];
}


#endif


#ifdef Py_REF_DEBUG
long g_refTotal;
// give out at most these many warnings.
int g_refWarningCount;

// reference connt debug
void saveRefCount()
{
	g_refTotal = _Py_RefTotal;
}


void checkRefCount()
{
	if (_Py_RefTotal > g_refTotal and g_refWarningCount-- > 0)
		cerr << "Warning: Ref count increased from " << g_refTotal << " to " << _Py_RefTotal
		     << "\nThis may be a sign of memory leak, especially when refCount increase"
		     << "\nindefinitely in a loop. Please contact simuPOP deceloper and report"
		     << "\nthe problem.\n" << endl;
	g_refTotal = _Py_RefTotal;
}


#endif

// Some common functions/templates

#if  defined (_WIN32) || defined (__WIN32__)
int simuPOP_kbhit()
{
	return _kbhit();
}


int simuPOP_getch(void)
{
	return getch();
}


#else

int simuPOP_kbhit(void)
{
	// tcsetattr is executed very slowly. To avoid
	// significant slow down, this function will be
	// called every 1 sec. The response time should
	// still be good.
	static time_t last_time = 0;

	time_t cur_time = time(NULL);

	if (cur_time != last_time)
		last_time = cur_time;
	else
		return 0;

	struct termios term, oterm;
	int fd = 0;
	int c = 0;

	tcgetattr(fd, &oterm);
	memcpy(&term, &oterm, sizeof(term));
	term.c_lflag = term.c_lflag & (!ICANON);
	term.c_cc[VMIN] = 0;
	term.c_cc[VTIME] = 1;
	tcsetattr(fd, TCSANOW, &term);
	c = getchar();
	tcsetattr(fd, TCSANOW, &oterm);
	if (c != -1)
		ungetc(c, stdin);
	return (c != -1) ? 1 : 0;
}


int simuPOP_getch(void)
{
	struct termios oldt, newt;
	int ch;

	tcgetattr(STDIN_FILENO, &oldt);
	newt = oldt;
	newt.c_lflag &= ~(ICANON | ECHO);
	tcsetattr(STDIN_FILENO, TCSANOW, &newt);
	ch = getchar();
	tcsetattr(STDIN_FILENO, TCSANOW, &oldt);

	return ch;
}


#endif
}


namespace std {
// how to output a dictionary
ostream & operator<<(ostream & out, const strDict & dict)
{
	out << "{";
	if (!dict.empty() ) {
		strDict::const_iterator it = dict.begin();

		out << it->first << ":" << it->second;

		for (++it; it != dict.end(); ++it)
			out << ", " << it->first << ":" << it->second;
	}
	out << "}";
	return out;
}


// how to output a dictionary
ostream & operator<<(ostream & out, const intDict & dict)
{
	out << "{";
	if (!dict.empty() ) {
		intDict::const_iterator it = dict.begin();

		out << it->first << ":" << it->second;

		for (++it; it != dict.end(); ++it)
			out << ", " << it->first << ":" << it->second;
	}
	out << "}";
	return out;
}


// in msvc.
unsigned pow3(unsigned n)
{
	unsigned res = 1;

	for (unsigned i = 0; i < n; ++i)
		res *= 3;
	return res;
}


}


namespace simuPOP {

// additional types
stringList::stringList(PyObject * obj) : m_elems(), m_allAvail(false)
{
	if (obj == NULL || obj == Py_None)
		m_allAvail = true;
	else if (PyBool_Check(obj))
		// accept True/False
		m_allAvail = obj == Py_True;
	else if (PyString_Check(obj))
		addString(obj);
	else if (PySequence_Check(obj)) {
		// assign values
		UINT numStr = PySequence_Size(obj);
		for (size_t i = 0; i < numStr; ++i) {
			PyObject * item = PySequence_GetItem(obj, i);
			addString(item);
			Py_DECREF(item);
		}
	}
}


void stringList::addString(PyObject * str)
{
	PyObject * res = PyObject_Str(str);

	if (res == NULL)
		return;
	string value = string(PyString_AsString(res));
	m_elems.push_back(value);
	Py_DECREF(res);
}


void stringList::obtainFrom(const stringList & items, const char * allowedItems[],
                            const char * defaultItems[])
{
	if (items.empty()) {
		for (size_t i = 0; defaultItems[i][0]; ++i)
			m_elems.push_back(defaultItems[i]);
		return;
	}
	stringList allowed;
	for (size_t i = 0; allowedItems[i][0]; ++i)
		allowed.m_elems.push_back(allowedItems[i]);
	vectorstr::const_iterator it = items.elems().begin();
	vectorstr::const_iterator itEnd = items.elems().end();
	for (; it != itEnd; ++it)
		if (allowed.contains(*it))
			m_elems.push_back(*it);
	if (m_elems.empty())
		for (size_t i = 0; defaultItems[i][0]; ++i)
			m_elems.push_back(defaultItems[i]);
}


// additional types
stringMatrix::stringMatrix(PyObject * obj) : m_elems()
{
	if (obj == NULL)
		return;
	DBG_ASSERT(PySequence_Check(obj), ValueError,
		"A list or a nested list of strings is expected");

	UINT numItems = PySequence_Size(obj);
	for (size_t i = 0; i < numItems; ++i) {
		PyObject * item = PySequence_GetItem(obj, i);
		if (PyString_Check(item)) {
			DBG_FAILIF(m_elems.size() > 1, ValueError,
				"A mixture of string and list is not allowed.")
			if (m_elems.empty())
				m_elems.push_back(vectorstr());
			string value = string(PyString_AsString(item));
			m_elems[0].push_back(value);
		} else if (PySequence_Check(item)) {
			m_elems.push_back(vectorstr());
			int numStrs = PySequence_Size(item);
			for (int j = 0; j < numStrs; ++j) {
				PyObject * str = PySequence_GetItem(item, j);
				DBG_ASSERT(PyString_Check(str), ValueError,
					"A list or nested list of string is expected");
				string value = string(PyString_AsString(str));
				m_elems.back().push_back(value);
				Py_DECREF(str);
			}
		} else {
			DBG_FAILIF(true, ValueError, "Can not create a string matrix from input.");
		}
		Py_DECREF(item);
	}
}


uintList::uintList(PyObject * obj) : m_elems(), m_allAvail(false)
{
	if (obj == NULL || obj == Py_None)
		// accept NULL
		m_allAvail = true;
	else if (PyBool_Check(obj))
		// accept True/False
		m_allAvail = obj == Py_True;
	else if (PyNumber_Check(obj)) {
		// accept a number
		m_allAvail = false;
		m_elems.push_back(static_cast<UINT>(PyInt_AS_LONG(obj)));
	} else if (PySequence_Check(obj)) {
		m_elems.resize(PySequence_Size(obj));
		// assign values
		for (size_t i = 0, iEnd = m_elems.size(); i < iEnd; ++i) {
			PyObject * item = PySequence_GetItem(obj, i);
			DBG_ASSERT(PyNumber_Check(item), ValueError, "Invalid input for a list of loci.");
			m_elems[i] = static_cast<UINT>(PyInt_AS_LONG(item));
			Py_DECREF(item);
		}
	} else {
		DBG_FAILIF(true, ValueError, "Invalid input for a list of loci.");
	}
}


intList::intList(PyObject * obj) : m_elems(), m_allAvail(false)
{
	if (obj == NULL || obj == Py_None)
		// accept NULL
		m_allAvail = true;
	else if (PyBool_Check(obj))
		// accept True/False
		m_allAvail = obj == Py_True;
	else if (PyNumber_Check(obj)) {
		// accept a number
		m_allAvail = false;
		m_elems.push_back(PyInt_AS_LONG(obj));
	} else if (PySequence_Check(obj)) {
		m_elems.resize(PySequence_Size(obj));
		// assign values
		for (size_t i = 0, iEnd = m_elems.size(); i < iEnd; ++i) {
			PyObject * item = PySequence_GetItem(obj, i);
			DBG_ASSERT(PyNumber_Check(item), ValueError, "Invalid input for a list of rep.");
			m_elems[i] = PyInt_AS_LONG(item);
			Py_DECREF(item);
		}
	} else {
		DBG_FAILIF(true, ValueError, "Invalid input for a list of rep.");
	}
}


//
// shared variables
//

// utility functions: python => C++ conversion
void PyObj_As_Bool(PyObject * obj, bool & val)
{
	if (obj == NULL || obj == Py_False || obj == Py_None) {
		val = false;
		return;
	} else if (obj == Py_True)
		val = true;
	else
		val = PyInt_AS_LONG(obj) ? true : false;
}


void PyObj_As_Int(PyObject * obj, long int & val)
{
	if (obj == NULL) {
		val = 0;
		return;
	}
	// try to convert
	PyObject * res = PyNumber_Int(obj);
	if (res == NULL)
		throw ValueError("Can not convert object to an integer");

	val = static_cast<int>(PyInt_AS_LONG(res));
	Py_DECREF(res);
}


void PyObj_As_Double(PyObject * obj, double & val)
{
	if (obj == NULL) {
		val = 0.0;
		return;
	}
	// try to convert
	PyObject * res = PyNumber_Float(obj);
	if (res == NULL)
		throw ValueError("Can not convert object to a double number");

	val = PyFloat_AsDouble(res);
	Py_DECREF(res);
}


void PyObj_As_String(PyObject * obj, string & val)
{
	if (obj == NULL) {
		val = string();
		return;
	}

	PyObject * res = PyObject_Str(obj);
	if (res == NULL)
		throw ValueError("Can not convert to a string");

	val = string(PyString_AsString(res));
	Py_DECREF(res);
}


void PyObj_As_Array(PyObject * obj, vectorf & val)
{
	if (obj == NULL) {
		val = vectorf();
		return;
	}
	if (PySequence_Check(obj)) {
		val.resize(PySequence_Size(obj));

		// assign values
		for (size_t i = 0, iEnd = val.size(); i < iEnd; ++i) {
			PyObject * item = PySequence_GetItem(obj, i);
			PyObj_As_Double(item, val[i]);
			Py_DECREF(item);
		}
	} else {
		val.resize(1);
		PyObj_As_Double(obj, val[0]);
	}
}


void PyObj_As_IntArray(PyObject * obj, vectori & val)
{
	if (obj == NULL) {
		val = vectori();
		return;
	}
	if (PySequence_Check(obj)) {
		val.resize(PySequence_Size(obj));

		// assign values
		for (size_t i = 0, iEnd = val.size(); i < iEnd; ++i) {
			PyObject * item = PySequence_GetItem(obj, i);
			PyObj_As_Int(item, val[i]);
			Py_DECREF(item);
		}
	} else {
		val.resize(1);
		PyObj_As_Int(obj, val[0]);
	}
}


bool PyObj_Is_IntNumArray(PyObject * obj)
{
	return is_carrayobject(obj) &&
	       carray_type(obj) == 'i' ;
}


bool PyObj_Is_DoubleNumArray(PyObject * obj)
{
	return is_carrayobject(obj) &&
	       carray_type(obj) == 'd';
}


bool PyObj_Is_AlleleNumArray(PyObject * obj)
{
	return is_carrayobject(obj) && carray_type(obj) == 'a';
}


PyObject * Double_Vec_As_NumArray(vectorf::iterator begin, vectorf::iterator end)
{
	PyObject * res = newcarrayobject(reinterpret_cast<char *>(& * begin), 'd', end - begin);

	DBG_FAILIF(res == NULL, ValueError, "Can not convert vector to double num array");
	return res;
}


PyObject * Int_Vec_As_NumArray(vectori::iterator begin, vectori::iterator end)
{
	PyObject * res = newcarrayobject(reinterpret_cast<char *>(& * begin), 'i', end - begin);

	DBG_FAILIF(res == NULL, ValueError, "Can not convert vector to int num array");
	return res;
}


PyObject * Allele_Vec_As_NumArray(GenoIterator begin, GenoIterator end)
{
	PyObject * res = newcarrayiterobject(begin, end);

	DBG_FAILIF(res == NULL, ValueError, "Can not convert buf to Allele num array");
	return res;
}


int NumArray_Size(PyObject * obj)
{
	// return PyArray_Size(obj);
	return carray_length(obj);
}


char * NumArray_Data(PyObject * obj)
{
	// return reinterpret_cast<PyArrayObject*>(obj)->data;
	return carray_data(obj);
}


// copy constructor
SharedVariables::SharedVariables(const SharedVariables & rhs)
	: m_ownVars(rhs.m_ownVars)
{
	if (rhs.m_ownVars) {
		m_dict = PyDict_New();
		PyObject * key = 0, * value = 0;

		Py_ssize_t i = 0;
		while (PyDict_Next(rhs.m_dict, &i, &key, &value)) {
			//Py_INCREF(key);
			//Py_INCREF(value);
			PyDict_SetItem(m_dict, key, value);
		}
	} else
		m_dict = rhs.m_dict;
}


SharedVariables::~SharedVariables()
{
	if (m_ownVars) {
		PyDict_Clear(m_dict);
		Py_DECREF(m_dict);
	}
}


// setvars C++ ==> Python
PyObject * SharedVariables::setVar(const string & name, const PyObject * val)
{
	/* find the first piece */
	size_t i, s;

	for (i = 0; i < name.size() && name[i] != '[' && name[i] != '{'; ++i) ;

	if (i == 0)
		throw ValueError("Empty name? " + name);

	// we need to keep current parent, key/index, type (0 for array, 1 for dict)
	// keytype (if curType = 1) (0 for string, 1 for num)
	int curType = 1;
	PyObject * curParent = m_dict;                            // should be always valid
	PyObject * curKey = PyString_FromString(const_cast<char *>(name.substr(0, i).c_str()));
	int curIdx = 0;
	PyObject * curChild = NULL;

next:
	// get par[1] (dict), curChild can be null, or borrow ref
	if (curType == 1)
		curChild = PyDict_GetItem(curParent, curKey);
	// get par[1] (list)
	else
		curChild = PyList_GetItem(curParent, curIdx);

	// subPop[0]{'alleleNum'}[0]
	// curParent:  m_dict
	// curChild:   m_dirt['subPop']

	if (name[i] == '[') {                                             // we need to put an array in ...
		// m_dict['subPop'][0]

		// look for index
		s = ++i;
		for ( ; name[i] != ']'; i++)
			if (!isdigit(name[i]) )
				throw ValueError("Expecting numbers: " + name);

		// get index
		size_t idx = atoi(name.substr(s, i - s).c_str());
		// not exist
		//
		// create subPop, ...
		if (curChild == NULL || !PyList_Check(curChild) ) {
			// Py_XDECREF(curChild);
			// create a sequence with this name
			curChild = PyList_New(idx + 1);
			// keep this
			// Py_INCREF(curChild);

			if (curChild == NULL)
				throw  ValueError("Can not create list: " + name);

			for (size_t a = 0; a < idx; ++a) {
				Py_INCREF(Py_None);
				PyList_SetItem(curChild, a, Py_None);
			}
			// now set that one
			if (curType == 1) {                                       // dictionary?
				PyDict_SetItem(curParent, curKey, curChild);
				Py_XDECREF(curChild);
			} else
				// PyList_SetItem steals ref
				PyList_SetItem(curParent, curIdx, curChild);
		} else {                                                                  // if exist, if it is a sequence?
			size_t curSize = PyList_Size(curChild);
			// if the size if enough, get the item
			if (curSize <= idx) {
				while (curSize++ <= idx) {
					Py_INCREF(Py_None);
					PyList_Append(curChild, Py_None);
				}
			}
		}
		// ready for iteration
		curType = 0;
		curParent = curChild;
		Py_XDECREF(curKey);
		curKey = NULL;
		curIdx = idx;
		i++;

		goto next;
	} else if (name[i] == '{') {                                      // dictionary
		//	keytype can be numeric 0: subPop{0}
		//                 tuple   1: subPop{(0,1)}
		//                 string  2: subPop{0}{'alleleNum'}
		int keyType;

		// look for index,
		s = ++i;
		if (name[s] == '\'' || name[s] == '\"')
			keyType = 2;
		else if (name[s] == '(')
			keyType = 1;
		else
			keyType = 0;

		for ( ; name[i] != '}' && i < name.size(); ++i) ;

		DBG_ASSERT(name[i] == '}', ValueError, "Unmatched dictionary delimiter");

		PyObject * childKey;

		if (keyType == 0)
			childKey = PyInt_FromString(const_cast<char *>(name.substr(s, i - s).c_str()), NULL, 0);
		else if (keyType == 1) {
			vectori key;
			for (size_t j = s + 1, k = j; j < i; j = k + 1) {
				for (k = j + 1; k < i && name[k] != ',' && name[k] != ')'; ++k) ;
				key.push_back(atoi(name.substr(j, k - j).c_str()));
			}
			childKey = PyTuple_New(key.size());
			for (size_t j = 0; j < key.size(); ++j)
				PyTuple_SetItem(childKey, j, PyInt_FromLong(key[j]));
		} else
			childKey = PyString_FromString(const_cast<char *>(name.substr(s + 1, i - s - 2).c_str()));
		// not exist
		if (curChild == NULL || !PyDict_Check(curChild) ) {
			// create dictionary with given key, append
			curChild = PyDict_New();
			// keep this

			if (curChild == NULL)
				throw  ValueError("Can not create dictionary" + name);

			// now set that one
			if (curType == 1) {                                       // dictionary?
				PyDict_SetItem(curParent, curKey, curChild);
				Py_XDECREF(curChild);
			} else
				PyList_SetItem(curParent, curIdx, curChild);
		}

		// ready for iteration
		curType = 1;
		// Py_XDECREF(curParent);
		curParent = curChild;
		Py_XDECREF(curKey);
		curKey = childKey;
		++i;
		goto next;
	} else {                                                                          // end here
		// CASE:  'res' = 1
		//
		// regardless curChild == NULL
		// set value in
		if (curType == 1) {                                               // dictionary?
			PyDict_SetItem(curParent, curKey, const_cast<PyObject *>(val));
			Py_XDECREF(const_cast<PyObject *>(val));
		} else
			PyList_SetItem(curParent, curIdx, const_cast<PyObject *>(val));
		Py_XDECREF(curKey);
	}
	return const_cast<PyObject *>(val);
}


PyObject * SharedVariables::getVar(const string & name, bool nameError)
{
	DBG_ASSERT(m_dict != NULL, ValueError,
		"Shared variables are not associated with any Python variable. You populaiton might not be part of a simulator.");

	// go deep in to the string
	size_t i, s;
	for (i = 0; i < name.size() && name[i] != '[' && name[i] != '{'; ++i) ;

	if (i == 0)
		throw ValueError("Empty name? " + name);

	int curType = 1;
	PyObject * curParent = m_dict;                            // should be always valid
	PyObject * curKey = PyString_FromString(const_cast<char *>(name.substr(0, i).c_str()));
	int curIdx = 0;
	PyObject * curChild;

next:
	if (curType == 1)
		curChild = PyDict_GetItem(curParent, curKey);
	else
		curChild = PyList_GetItem(curParent, curIdx);

	if (curChild == NULL) {
		Py_XDECREF(curKey);
		// name does not exist
		if (nameError)
			throw ValueError("Shared variable name '" + name + "' does not exist.");
		else
			return NULL;
	}

	if (name[i] == '[') {                                             // we need to put an array in ...
		// get key
		// look for index
		s = ++i;
		for ( ; name[i] != ']'; i++)
			if (!isdigit(name[i]) )
				throw ValueError("Expecting numbers: " + name);

		// get index
		int idx = atoi(name.substr(s, i - s).c_str());

		curType = 0;
		curParent = curChild;
		Py_XDECREF(curKey);
		curKey = NULL;
		curIdx = idx;
		i++;
		goto next;
	} else if (name[i] == '{') {
		//	keytype can be numeric 0: subPop{0}
		//                 tuple   1: subPop{(0,1)}
		//                 string  2: subPop{0}{'alleleNum'}
		int keyType;

		// look for index,
		s = ++i;
		if (name[s] == '\'' || name[s] == '\"')
			keyType = 2;
		else if (name[s] == '(')
			keyType = 1;
		else
			keyType = 0;

		for ( ; name[i] != '}' && i < name.size(); ++i) ;

		assert(name[i] == '}');

		PyObject * childKey;
		if (keyType == 0)
			childKey = PyInt_FromString(const_cast<char *>(name.substr(s, i - s).c_str()), NULL, 0);
		else if (keyType == 1) {
			vectori key;
			for (size_t j = s + 1, k = j; j < i; j = k + 1) {
				for (k = j + 1; k < i && name[k] != ',' && name[k] != ')'; ++k) ;
				key.push_back(atoi(name.substr(j, k - j).c_str()));
			}
			childKey = PyTuple_New(key.size());
			for (size_t j = 0; j < key.size(); ++j)
				PyTuple_SetItem(childKey, j, PyInt_FromLong(key[j]));
		} else
			childKey = PyString_FromString(const_cast<char *>(name.substr(s + 1, i - s - 2).c_str()));

		// ready for iteration
		curType = 1;
		curParent = curChild;
		Py_XDECREF(curKey);
		curKey = childKey;
		++i;
		goto next;
	} else {
		Py_XDECREF(curKey);
		return curChild;
	}
}


void SharedVariables::removeVar(const string & name)
{
	DBG_ASSERT(m_dict != NULL, ValueError,
		"Shared variables are not associated with any Python variable. You populaiton might not be part of a simulator.");

	DBG_DO(DBG_UTILITY, cerr << "Removing variable " << name << endl);

	// go deep in to the string
	size_t i, s;
	for (i = 0; i < name.size() && name[i] != '[' && name[i] != '{'; ++i) ;

	if (i == 0)
		throw ValueError("Empty name? " + name);

	int curType = 1;
	PyObject * curParent = m_dict;                            // should be always valid
	PyObject * curKey = PyString_FromString(const_cast<char *>(name.substr(0, i).c_str()));
	int curIdx = 0;
	PyObject * curChild;

next:
	if (curType == 1)
		curChild = PyDict_GetItem(curParent, curKey);
	else {
		if (curIdx < PyList_Size(curParent))
			curChild = PyList_GetItem(curParent, curIdx);
		else
			return;
	}

	// maybe the item has been removed?
	if (curChild == NULL) {
		Py_XDECREF(curKey);
		return;
	}

	if (name[i] == '[') {                                             // we need to put an array in ...
		// get key
		// look for index
		s = ++i;
		for ( ; name[i] != ']'; i++)
			if (!isdigit(name[i]) )
				throw ValueError("Expecting numbers: " + name);

		// get index
		int idx = atoi(name.substr(s, i - s).c_str());

		curType = 0;
		curParent = curChild;
		Py_XDECREF(curKey);
		curKey = NULL;
		curIdx = idx;
		i++;
		goto next;
	} else if (name[i] == '{') {
		int keyType;

		// look for index,
		s = ++i;
		if (name[s] == '\'' || name[s] == '\"')
			keyType = 2;
		else if (name[s] == '(')
			keyType = 1;
		else
			keyType = 0;

		for ( ; name[i] != '}' && i < name.size(); ++i) ;

		assert(name[i] == '}');

		PyObject * childKey;

		if (keyType == 0)
			childKey = PyInt_FromString(const_cast<char *>(name.substr(s, i - s).c_str()), NULL, 0);
		else if (keyType == 1) {
			vectori key;
			for (size_t j = s + 1, k = j; j < i; j = k + 1) {
				for (k = j + 1; k < i && name[k] != ',' && name[k] != ')'; ++k) ;
				key.push_back(atoi(name.substr(j, k - j).c_str()));
			}
			childKey = PyTuple_New(key.size());
			for (size_t j = 0; j < key.size(); ++j)
				PyTuple_SetItem(childKey, j, PyInt_FromLong(key[j]));
		} else
			childKey = PyString_FromString(const_cast<char *>(name.substr(s + 1, i - s - 2).c_str()));

		// ready for iteration
		curType = 1;
		curParent = curChild;
		Py_XDECREF(curKey);
		curKey = childKey;
		++i;
		goto next;
	} else {                                                                          // now the last piece, we know the parents
		if (curType == 1)
			PyDict_DelItem(curParent, curKey);
		else {
			// to avoid affecting order of others, set to None
			Py_INCREF(Py_None);
			PyList_SetItem(curParent, curIdx, Py_None);
		}
		Py_XDECREF(curKey);
		return;
	}

}


PyObject * SharedVariables::setBoolVar(const string & name, const bool val)
{
	PyObject * obj = val ? Py_True : Py_False;

	Py_INCREF(obj);
	return setVar(name, obj);
}


PyObject * SharedVariables::setIntVar(const string & name, const int val)
{
	return setVar(name, PyInt_FromLong(val));
}


PyObject * SharedVariables::setDoubleVar(const string & name, const double val)
{
	return setVar(name, PyFloat_FromDouble(val));
}


PyObject * SharedVariables::setStringVar(const string & name, const string & val)
{
	return setVar(name, Py_BuildValue("s", name.c_str()));
}


PyObject * SharedVariables::setIntVectorVar(const string & name, const vectori & val)
{
	PyObject * obj = PyList_New(0);
	PyObject * item;

	for (vectori::const_iterator it = val.begin();
	     it < val.end(); ++it) {
		item = PyInt_FromLong(*it);
		PyList_Append(obj, item);
		Py_XDECREF(item);
	}
	return setVar(name, obj);
}


//CPPONLY
PyObject * SharedVariables::setDoubleVectorVar(const string & name, const vectorf & val)
{
	PyObject * obj = PyList_New(0);
	PyObject * item;

	for (vectorf::const_iterator it = val.begin();
	     it < val.end(); ++it) {
		item = PyFloat_FromDouble(*it);
		PyList_Append(obj, item);
		Py_XDECREF(item);
	}
	return setVar(name, obj);
}


PyObject * SharedVariables::setStrDictVar(const string & name, const strDict & val)
{
	PyObject * obj = PyDict_New();
	PyObject * v;

	for (strDict::const_iterator i = val.begin(); i != val.end(); ++i) {
		PyDict_SetItemString(obj, const_cast<char *>(i->first.c_str()),
			v = PyFloat_FromDouble(i->second));
		Py_XDECREF(v);
	}
	return setVar(name, obj);
}


PyObject * SharedVariables::setIntDictVar(const string & name, const intDict & val)
{
	PyObject * obj = PyDict_New();
	PyObject * u, * v;

	for (intDict::const_iterator i = val.begin(); i != val.end(); ++i) {
		PyDict_SetItem(obj,
			u = PyInt_FromLong(i->first),
			v = PyFloat_FromDouble(i->second));
		Py_XDECREF(u);
		Py_XDECREF(v);
	}
	return setVar(name, obj);
}


PyObject * SharedVariables::setIntDefDictVar(const string & name, const intDict & val)
{
	PyObject * obj = PyDefDict_New();
	PyObject * u, * v;

	for (intDict::const_iterator i = val.begin(); i != val.end(); ++i) {
		PyObject_SetItem(obj,
			u = PyInt_FromLong(i->first),
			v = PyFloat_FromDouble(i->second));
		Py_XDECREF(u);
		Py_XDECREF(v);
	}
	return setVar(name, obj);
}


PyObject * SharedVariables::setTupleDefDictVar(const string & name, const tupleDict & val)
{
	PyObject * obj = PyDefDict_New();
	PyObject * u, * v;

	tupleDict::const_iterator it = val.begin();
	tupleDict::const_iterator itEnd = val.end();

	for (; it != itEnd; ++it) {
		const vectori & key = it->first;
		u = PyTuple_New(key.size());
		for (size_t i = 0; i < key.size(); ++i)
			PyTuple_SetItem(u, i, PyInt_FromLong(key[i]));
		v = PyFloat_FromDouble(it->second);
		PyObject_SetItem(obj, u, v);
		Py_XDECREF(u);
		Py_XDECREF(v);
	}
	return setVar(name, obj);
}


void save_none(string & str)
{
	str += "n";
}


PyObject * load_none(const string & str, size_t & offset)
{
	offset++;
	Py_INCREF(Py_None);
	return Py_None;
}


void save_int(string & str, PyObject * args)
{
	long l = PyInt_AS_LONG((PyIntObject *)args);

	// type + string + ' '
	str += 'i' + toStr(l) + ' ';
}


PyObject * load_int(const string & str, size_t & offset)
{
	// search for blank
	int len = 0;

	while (str[offset + len + 1] != ' ') len++;
	PyObject * val = PyInt_FromString(
		const_cast<char *>(str.substr(offset + 1, len).c_str()),
		NULL, 0);
	offset += len + 2;
	return val;
}


void save_long(string & str, PyObject * args)
{
	long l = PyInt_AS_LONG(args);

	// type +  string + ' '
	str += 'i' + toStr(l) + ' ';
}


PyObject * load_long(const string & str, size_t & offset)
{
	int len = 0;

	while (str[offset + len + 1] != ' ') len++;
	PyObject * val = PyLong_FromString(
		const_cast<char *>(str.substr(offset + 1, len).c_str()),
		NULL, 0);
	offset += len + 2;
	return val;
}


void save_float(string & str, PyObject * args)
{
	double d = PyFloat_AsDouble(args);

	// type + string
	str += 'f' + toStr(d) + ' ';
}


PyObject * load_float(const string & str, size_t & offset)
{
	int len = 0;

	while (str[offset + len + 1] != ' ') len++;
	double d = atof(const_cast<char *>(str.substr(offset + 1, len).c_str()));
	offset += len + 2;
	return PyFloat_FromDouble(d);
}


void save_string(string & str, PyObject * args)
{
	char * s = PyString_AsString(args);

	str += 's' + string(s) + '\0';
}


PyObject * load_string(const string & str, size_t & offset)
{
	char * s = const_cast<char *>(str.c_str()) + offset + 1;
	size_t len = strlen(s);

	offset += len + 2;
	DBG_ASSERT(str[offset - 1] == '\0', SystemError,
		"Error while loading string from str");

	return PyString_FromString(s);
}


void saveObj(string & str, PyObject * args);

PyObject * loadObj(const string & vars, size_t & offset);

void save_dict(string & str, PyObject * args)
{
	PyObject * key = 0, * value = 0;

	str += 'd';                                                                       // dictionary
	Py_ssize_t i = 0;
	while (PyDict_Next(args, &i, &key, &value)) {
		saveObj(str, key);
		saveObj(str, value);
	}
	str += 'e';                                                                       // ending of a dictionary
}


PyObject * load_dict(const string & vars, size_t & offset)
{
	// skip 'd'
	offset++;
	PyObject * d = PyDict_New();
	while (vars[offset] != 'e') {
		// get key
		PyObject * key = loadObj(vars, offset);
		// get value
		PyObject * val = loadObj(vars, offset);
		PyDict_SetItem(d, key, val);
		// Is this needed???
		Py_DECREF(key);
		Py_DECREF(val);
	}
	offset++;                                                                         // skip 'e'
	return d;
}


void save_defdict(string & str, PyObject * args)
{
	PyObject * key = 0, * value = 0;

	str += 'D';                                                                       // dictionary
	Py_ssize_t i = 0;
	while (PyDict_Next(args, &i, &key, &value)) {
		saveObj(str, key);
		saveObj(str, value);
	}
	str += 'e';                                                                       // ending of a dictionary
}


PyObject * load_defdict(const string & vars, size_t & offset)
{
	// skip 'D'
	offset++;
	PyObject * d = PyDefDict_New();
	while (vars[offset] != 'e') {
		// get key
		PyObject * key = loadObj(vars, offset);
		// get value
		PyObject * val = loadObj(vars, offset);
		PyObject_SetItem(d, key, val);
		// Is this needed???
		Py_DECREF(key);
		Py_DECREF(val);
	}
	offset++;                                                                         // skip 'e'
	return d;
}


void save_list(string & str, PyObject * args)
{
	str += 'L';                                                                       // dictionary
	int len = PyList_Size(args);
	for (int i = 0; i < len; i++) {
		PyObject * elem = PyList_GET_ITEM((PyListObject *)args, i);
		saveObj(str, elem);
	}

	str += 'e';                                                                       // ending of a dictionary
}


PyObject * load_list(const string & vars, size_t & offset)
{
	// skip 'L'
	offset++;
	PyObject * d = PyList_New(0);
	while (vars[offset] != 'e') {
		PyObject * elem = loadObj(vars, offset);
		PyList_Append(d, elem);
		// There has been a weird bug related to this. Basically,
		// There is a loaded list [1154, b, c] and a number 1154 somewhere
		// However, the loaded list does not increase reference count
		// so the ref count of 1154 is one. Then, when the list element
		// is replaced, the other instance of 1154 is also changed...
		//
		Py_DECREF(elem);
	}
	offset++;                                                                         // skip 'e'
	return d;
}


void save_tuple(string & str, PyObject * args)
{
	str += 't';                                                                       // dictionary
	int len = PyTuple_Size(args);
	// save length
	str += toStr(len) + ' ';
	// save items
	for (int i = 0; i < len; i++) {
		PyObject * elem = PyTuple_GET_ITEM((PyTupleObject *)args, i);
		saveObj(str, elem);
	}
}


PyObject * load_tuple(const string & vars, size_t & offset)
{
	// skip 't'
	offset++;
	// search for blank
	int l = 0;
	while (vars[offset + l] != ' ') l++;
	int len = atoi(vars.substr(offset, l).c_str());
	offset += l + 1;
	PyObject * d = PyTuple_New(len);
	for (int i = 0; i < len; ++i) {
		PyObject * elem = loadObj(vars, offset);
		PyTuple_SetItem(d, i, elem);
		Py_DECREF(elem);
	}
	return d;
}


void saveObj(string & str, PyObject * args)
{
	PyTypeObject * type;

	if (args == Py_None) {
		save_none(str);
		return;
	}

	type = args->ob_type;
	if (type == &PyInt_Type)
		save_int(str, args);
	else if (type == &PyDict_Type)
		save_dict(str, args);
	else if (type == &PyString_Type)
		save_string(str, args);
	else if (type == &PyLong_Type)
		save_long(str, args);
	else if (type == &PyList_Type)
		save_list(str, args);
	else if (type == &PyTuple_Type)
		save_tuple(str, args);
	else if (type == &PyFloat_Type)
		save_float(str, args);
	else if (is_defdict(type))
		save_defdict(str, args);
	else {
		// some other unknown type
		DBG_DO(DBG_UTILITY, cerr << "Warning: object of type '" + toStr(type->tp_name) + "' cannot be saved. Use none.");
		save_none(str);
	}
}


PyObject * loadObj(const string & vars, size_t & offset)
{
	switch (vars[offset]) {
	case 'd':                                                                         //
		return load_dict(vars, offset);
	case 'D':
		return load_defdict(vars, offset);
	case 'i':
		return load_int(vars, offset);
	case 's':
		return load_string(vars, offset);
	case 'f':
		return load_float(vars, offset);
	case 'l':
		return load_long(vars, offset);
	case 'L':
		return load_list(vars, offset);
	case 't':
		return load_tuple(vars, offset);
	case 'n':
		return load_none(vars, offset);
	default:
	{
		DBG_DO(DBG_UTILITY, cerr << endl);
		throw ValueError("Unknown type code at offset " + toStr(offset));
	}
	}
}


string SharedVariables::asString() const
{
	// go through each variable and save
	string str;

	saveObj(str, m_dict);
	str += 'e';                                                                       // ending character
	return str;
}


void SharedVariables::fromString(const string & vars)
{
	size_t offset = 0;
	PyObject * obj = loadObj(vars, offset);

	DBG_ASSERT(obj == NULL || vars[offset] == 'e', SystemError,
		"Failed to load objects from string");

	// remove m_dict
	if (m_ownVars) {
		PyDict_Clear(m_dict);
		Py_XDECREF(m_dict);
	}
	m_ownVars = true;
	m_dict = obj;
}


// simuVars will hold replicate specific Shared Variables.

// global dictionary
// will be set by // initialize() function
// DO NOT OWN the dictionaries
SharedVariables g_main_vars, g_module_vars;
swig_type_info * g_swigPopType, * g_swigindividual;

SharedVariables & mainVars()
{
	return g_main_vars;
}


SharedVariables & moduleVars()
{
	return g_module_vars;
}


PyObject * pyPopObj(void * p)
{
	return SWIG_NewPointerObj(p, g_swigPopType, 0);
}


PyObject * pyIndObj(void * p)
{
	return SWIG_NewPointerObj(p, g_swigindividual, 0);
}


void * pyIndPointer(PyObject * obj)
{
	void * ptr = 0;

	SWIG_Python_ConvertPtr(obj, &ptr, g_swigindividual, 0);
	return ptr;
}


// Expression evaluation
// because of ref count, need to define copier

Expression::Expression(const Expression & rhs)
	: m_expr(rhs.m_expr), m_stmts(rhs.m_stmts), m_locals(rhs.m_locals)
{
	Py_XINCREF(m_expr);
	Py_XINCREF(m_stmts);
}


Expression::~Expression()
{
	// release compiled code
	Py_XDECREF(m_expr);
	Py_XDECREF(m_stmts);
}


void Expression::compileExpr(const string & expr)
{
	if (m_expr != NULL) {
		// discard previous
		Py_XDECREF(m_expr);
		m_expr = NULL;
	}

	if (expr.empty())
		return;

	m_expr = Py_CompileString(const_cast<char *>(expr.c_str()), "<embed>",
		Py_eval_input);

	if (m_expr == NULL)
		throw ValueError("Expression '" + expr + "' is not a valid python expression.");
}


void Expression::compileStmts(const string & stmts)
{
	if (m_stmts != NULL) {
		// discard previous statements
		Py_XDECREF(m_stmts);
		m_stmts = NULL;
	}

	if (stmts.empty())
		return;

	// clear leading blank etc for statements ? Not yet.
	m_stmts = Py_CompileString(const_cast<char *>(stmts.c_str()), "<embed>",
		Py_file_input);

	if (m_stmts == NULL)
		throw ValueError("statement '" + stmts + "' is not valid.");
}


// python expression
PyObject * Expression::evaluate()
{
	if (m_expr == NULL && m_stmts == NULL)
		return NULL;

	DBG_ASSERT(mainVars().dict() != NULL && m_locals != NULL,
		ValueError, "Can not evalulate. Dictionary is empty!");

	PyObject * res = NULL;
	if (m_stmts != NULL) {
		res = PyEval_EvalCode((PyCodeObject *)m_stmts,
			mainVars().dict(), m_locals);
		if (res == NULL) {
#ifndef OPTIMIZED
			if (debug(DBG_GENERAL)) {
				PyErr_Print();
				PyErr_Clear();
			}
#endif
			throw RuntimeError("Evalulation of statements failed");
		} else {
			Py_DECREF(res);
			res = NULL;
		}
	}

	if (m_expr != NULL) {
		res = PyEval_EvalCode((PyCodeObject *)m_expr,
			mainVars().dict(), m_locals);
		if (res == NULL) {
#ifndef OPTIMIZED
			if (debug(DBG_GENERAL)) {
				PyErr_Print();
				PyErr_Clear();
			}
#endif
			throw RuntimeError("Evalulation of expression failed");
		}
	}

	return res;
}


bool Expression::valueAsBool()
{
	PyObject * res = evaluate();

	if (res == NULL)
		return false;
	bool val;
	PyObj_As_Bool(res, val);
	Py_XDECREF(res);
	return val;
}


long int Expression::valueAsInt()
{
	PyObject * res = evaluate();

	if (res == NULL)
		return 0;
	long int val;
	PyObj_As_Int(res, val);
	Py_XDECREF(res);
	return val;
}


double Expression::valueAsDouble()
{
	PyObject * res = evaluate();

	if (res == NULL)
		return 0;
	double val;
	PyObj_As_Double(res, val);
	Py_XDECREF(res);
	return val;
}


string Expression::valueAsString()
{
	PyObject * res = evaluate();

	if (res == NULL)
		return "";
	string val;
	PyObj_As_String(res, val);
	Py_XDECREF(res);
	return val;
}


vectorf Expression::valueAsArray()
{
	PyObject * res = evaluate();

	if (res == NULL)
		return vectorf();
	vectorf val;
	PyObj_As_Array(res, val);
	Py_XDECREF(res);
	return val;
}


simpleStmt::simpleStmt(const string & stmt) : m_var(""),
	m_operation(NoOperation), m_value(0)
{
	const regex assign("\\s*([\\d\\w]+)\\s*=\\s*([-]*[1-9\\.]+)\\s*");
	const regex increment1("\\s*([\\d\\w]+)\\s*\\+=\\s*([-]*[1-9\\.]+)\\s*");
	const regex increment2("\\s*([\\d\\w]+)\\s*=\\s*\\1\\s*\\+\\s*([-]*[1-9\\.]+)\\s*");
	const regex increment3("\\s*([\\d\\w]+)\\s*=\\s*([-]*[1-9\\.]+)\\s*\\+\\s*\\1\\s*");
	const regex decrement1("\\s*([\\d\\w]+)\\s*-=\\s*([-]*[1-9\\.]+)\\s*");
	const regex decrement2("\\s*([\\d\\w]+)\\s*=\\s*\\1\\s*-\\s*([-]*[1-9\\.]+)\\s*");
	const regex multiplied1("\\s*([\\d\\w]+)\\s*\\*=\\s*([-]*[1-9\\.]+)\\s*");
	const regex multiplied2("\\s*([\\d\\w]+)\\s*\\=\\s*\\1\\s*\\*\\s*([-]*[1-9\\.]+)\\s*");
	const regex multiplied3("\\s*([\\d\\w]+)\\s*\\=\\s*\\s*([-]*[1-9\\.]+)\\s*\\*\\s*\\1\\s*");

	cmatch matches;

	DBG_DO(DBG_DEVEL, cerr << "Decipher statement " << stmt << endl);
	if (regex_match(stmt.c_str(), matches, assign))
		m_operation = Assignment;
	else if (regex_match(stmt.c_str(), matches, increment1) ||
	         regex_match(stmt.c_str(), matches, increment2) ||
	         regex_match(stmt.c_str(), matches, increment3))
		m_operation = Increment;
	else if (regex_match(stmt.c_str(), matches, decrement1) ||
	         regex_match(stmt.c_str(), matches, decrement2))
		m_operation = Decrement;
	else if (regex_match(stmt.c_str(), matches, multiplied1) ||
	         regex_match(stmt.c_str(), matches, multiplied2) ||
	         regex_match(stmt.c_str(), matches, multiplied3))
		m_operation = MultipliedBy;
	else
		return;

	m_var = string(matches[1].first, matches[1].second);
	try {
		m_value = boost::lexical_cast<double>(string(matches[2].first, matches[2].second));
	} catch (...) {
		m_operation = NoOperation;
		return;
	}
	DBG_DO(DBG_DEVEL, cerr << "Match statement with name " << m_var
		                   << " and value " << m_value << " with operation " << m_operation << endl);
}


// Stream element, can be of different types

StreamElem::StreamElem(const string & name, bool readable, bool realAppend, bool useString)
	: m_filename(name)
{

	DBG_DO(DBG_UTILITY, cerr << "creating " << name << " with parameter " <<
		readable << " " << realAppend << " " << useString << " " << endl);

	if (useString) {
		// existing file will be truncated...
		m_stream = new stringstream(std::ios::in | std::ios::out);
		m_type = SSTREAM;
		m_append = false;
	} else {                                                                          // no string, usual file
		m_append = realAppend;
		m_type = FSTREAM;
		if (readable) {
			if (realAppend)     // readable , append
				m_stream = new fstream(name.c_str(),  std::ios::in | std::ios::out | std::ios::ate);
			else                // readable, ! append
				// existing file will be truncated...
				m_stream = new fstream(name.c_str(),  std::ios::in | std::ios::trunc | std::ios::out);
		} else {
			if (realAppend) {     // ! readable, append
				m_stream = new ofstream(name.c_str(), std::ios::out | std::ios::app);
				m_type = OFSTREAM;
			} else                //  !readable, !append )
				// existing file will be truncated...
				m_stream = new fstream(name.c_str(),  std::ios::trunc | std::ios::out);
		}
	}

	if (m_stream == NULL || !*m_stream)
		throw ValueError("Can not open specified file:" + name);

	DBG_DO(DBG_UTILITY, cerr << "New file info: " << info() << endl);

}


// copy constructor, we need to clear rhs.m_stream to avoid closing file too early
// this is techniquely advanced (and dangerous)
StreamElem::StreamElem(const StreamElem & rhs)
{
	m_filename = rhs.m_filename;
	m_type = rhs.m_type;
	m_append = rhs.m_append;
	m_stream = rhs.m_stream;
	const_cast<StreamElem &>(rhs).m_stream = NULL;
}


// destructor
StreamElem::~StreamElem()
{
	// close file.
	if (m_stream != NULL) {

		DBG_DO(DBG_UTILITY, cerr << "Closing file " << m_filename << endl);

		if (m_type == OFSTREAM)
			static_cast<ofstream *>(m_stream)->close();
		else if (m_type == FSTREAM)
			static_cast<fstream *>(m_stream)->close();

		// do not do anything for fstream (no close function

		delete m_stream;
	}
}


// the file was write-only, re-open it as read-write
void StreamElem::makeReadable()
{

	DBG_DO(DBG_UTILITY, cerr << "File was opened write-only. Re-open it.  " << info() << endl);

	static_cast<ofstream *>(m_stream)->close();

	// have to re-create a stream since pointer type is different.
	delete m_stream;

	// try to keep file content
	m_stream = new fstream(m_filename.c_str(),  std::ios::in | std::ios::out | std::ios::ate);

	if (m_stream == NULL || !*m_stream)
		throw ValueError("Can not re-open specified file.");

	m_type = FSTREAM;
}


// change the append status
void StreamElem::makeAppend(bool append)
{

	DBG_DO(DBG_UTILITY, cerr << "File append status changes to " << append << endl);

	DBG_FAILIF(m_type == SSTREAM, ValueError, "String stream can not be maded appendable. ");

	m_append = append;

	// if append = true: no problem, keep writing just do not close file. otherwise
	// reopen file. Otherwise,
	if (!append) {

		DBG_DO(DBG_UTILITY, cerr << "Re-open the file " << endl);

		// re-open the file.

		if (m_type == FSTREAM) {
			static_cast<fstream *>(m_stream)->close();
			static_cast<fstream *>(m_stream)->open(m_filename.c_str(), std::ios::in | std::ios::out | std::ios::trunc);
		} else if (m_type == OFSTREAM) {
			static_cast< ofstream *>(m_stream)->close();
			static_cast< ofstream *>(m_stream)->open(m_filename.c_str(), std::ios::out | std::ios::trunc);
		}
	}
}


string StreamElem::info()
{
	ostringstream out;

	switch (m_type) {
	case OFSTREAM:
		out << m_filename << " : write only file stream. " << endl;

		DBG_DO(DBG_UTILITY, out << "(write pos: " << static_cast<ofstream *>(m_stream)->tellp() << ")");

		break;
	case FSTREAM:
		out << m_filename << " : read-write file stream ";

		DBG_DO(DBG_UTILITY, out << "(write pos: " << static_cast<fstream *>(m_stream)->tellp()
			                    << ", read pos: " << static_cast<fstream *>(m_stream)->tellg() << ")");

		break;
	case SSTREAM:
		out << m_filename << " : string stream. " << endl;

		DBG_DO(DBG_UTILITY, out << "(write pos: " << static_cast<stringstream *>(m_stream)->tellp()
			                    << ", read pos: " << static_cast<stringstream *>(m_stream)->tellg() << ")");

		break;
	}
	return out.str();
}


// OStream Manager

OstreamManager::OstreamManager() : m_ostreams()
{
}


OstreamManager::~OstreamManager()
{
	closeAll();
}


ostream * OstreamManager::getOstream(const string & name, bool readable,  bool realAppend, bool useString)
{

	ostreamMapIterator it = m_ostreams.find(name);

	if (it == m_ostreams.end() ) {                            // not found.

		DBG_DO(DBG_UTILITY, cerr << "Create new file " << name << endl);

		return m_ostreams.insert(ostreamMapValue(name,
				StreamElem(name, readable, realAppend, useString))).first->second.stream();
	} else {                                                                          // already exist

		DBG_DO(DBG_UTILITY, cerr << "Find existing ostream " << name << " of info " << it->second.info() << endl);

		// try to see if existing stream type matches what is requrested.
		if (useString && it->second.type() != StreamElem::SSTREAM)
			throw ValueError("file " + name + " has already opened as string file.");
		else if (!useString && it->second.type() == StreamElem::SSTREAM)
			throw ValueError("file " + name + " has already opened as normal file.");
		else if (readable && it->second.type() == StreamElem::OFSTREAM)
			it->second.makeReadable();
		else if (realAppend && !it->second.append() )
			it->second.makeAppend(true);
		else if (!realAppend && it->second.append() )
			it->second.makeAppend(false);

		return it->second.stream();
	}

	// will never reach here.
	return NULL;
}


bool OstreamManager::hasOstream(const string & filename)
{
	return m_ostreams.end() != m_ostreams.find(filename) ;
}


void OstreamManager::listAll()
{
	for (ostreamMapIterator it = m_ostreams.begin(), itEnd = m_ostreams.end(); it != itEnd;  ++it)
		cerr << it->first << " : " << it->second.info() << endl;
}


void OstreamManager::closeOstream(const string & filename)
{
	if (!hasOstream(filename))
		return;
	m_ostreams.erase(filename);
}


// close all files and clean the map
void OstreamManager::closeAll()
{
	m_ostreams.clear();
}


// global ostream  manager
OstreamManager g_ostreams;

// return ostream manager
OstreamManager & ostreamManager()
{
	return g_ostreams;
}


// Stream provider

// all flags will be cleared to 0
StreamProvider::StreamProvider(const string & output, const pyFunc & func)
	: m_filename(output), m_filenameExpr(), m_func(func), m_flags(0), m_filePtr(NULL)
{
	if (!m_filename.empty() && m_filename[0] == '!') {
		m_filenameExpr.setExpr(m_filename.substr(1));
		m_filename.clear();
	} else if (m_func.isValid()) {
		SETFLAG(m_flags, m_flagUseFunc);
		SETFLAG(m_flags, m_flagCloseAfterUse);
	} else
		analyzeOutputString(m_filename);
}


ostream & StreamProvider::getOstream(PyObject * dict, bool readable)
{
	DBG_FAILIF(readable && (ISSETFLAG(m_flags, m_flagNoOutput) || ISSETFLAG(m_flags, m_flagUseDefault)),
		SystemError, "A readable file is requested but this Opertor uses cerr or cnull.");

	if (ISSETFLAG(m_flags, m_flagNoOutput))
		return cnull();

	// if using cerr, return it.
	if (ISSETFLAG(m_flags, m_flagUseDefault))
		return cerr;

	if (ISSETFLAG(m_flags, m_flagUseFunc)) {
		m_filePtr = new ostringstream();
		return *m_filePtr;
	}

	// if use and close type
	string filename;
	if (!m_filenameExpr.empty()) {
		DBG_ASSERT(dict != NULL, ValueError,
			"Need to know local dictionary to evaluate filename expression");
		m_filenameExpr.setLocalDict(dict);
		m_filename = m_filenameExpr.valueAsString();

		DBG_DO(DBG_UTILITY, cerr << "Filename " << m_filename << endl);

		analyzeOutputString(m_filename);
		if (ISSETFLAG(m_flags, m_flagNoOutput) )
			return cnull();

		// if using cerr, return it.
		if (ISSETFLAG(m_flags, m_flagUseDefault) )
			return cerr;
	}
	filename = m_filename;

	if (ISSETFLAG(m_flags, m_flagAppend) ) {

		DBG_DO(DBG_UTILITY, cerr << "Get a persistent file: "
			                     << filename << endl);

		return *ostreamManager().getOstream(filename, readable,
			ISSETFLAG(m_flags, m_flagRealAppend), ISSETFLAG(m_flags, m_flagUseString));
	} else {                                                                          // not in append mode, but check if this file is alreay there

		DBG_DO(DBG_UTILITY, cerr << "File is not persistent : "
			                     << filename << endl);

		if (!ostreamManager().hasOstream(filename) ) {
			if (readable)
				SETFLAG(m_flags, m_flagReadable);
			else
				RESETFLAG(m_flags, m_flagReadable);

			if (readable)
				m_filePtr = new fstream(filename.c_str() );
			else
				m_filePtr = new ofstream(filename.c_str() );

			if (m_filePtr == NULL || !*m_filePtr)
				throw SystemError("Can not create file " + filename);

			return *m_filePtr;
		} else {

			DBG_DO(DBG_UTILITY, cerr << "file " + filename +
				" is already opened as appendable. Use that file instead." << endl);

			RESETFLAG(m_flags, m_flagCloseAfterUse);
			SETFLAG(m_flags, m_flagAppend);
			return *ostreamManager().getOstream(filename, readable, ISSETFLAG(m_flags, m_flagRealAppend),
				ISSETFLAG(m_flags, m_flagUseString));
		}
	}
}


// close ostream and delete ostream pointer.. if it is a ofstream.
void StreamProvider::closeOstream()
{
	if (ISSETFLAG(m_flags, m_flagCloseAfterUse)) {
		if (ISSETFLAG(m_flags, m_flagUseFunc)) {
			DBG_ASSERT(m_func.isValid(), SystemError,
				"Passed function object is invalid");
			string str = dynamic_cast<ostringstream *>(m_filePtr)->str();
			PyObject * arglist = Py_BuildValue("(s)", str.c_str());
			PyObject * pyResult = PyEval_CallObject(m_func.func(), arglist);
			if (pyResult == NULL) {
				PyErr_Print();
				PyErr_Clear();
				throw RuntimeError("Failed to send output to a function.");
			} else
				Py_DECREF(pyResult);
		} else if (ISSETFLAG(m_flags, m_flagReadable))
			dynamic_cast<fstream *>(m_filePtr)->close();
		else
			dynamic_cast<ofstream *>(m_filePtr)->close();
		delete m_filePtr;
	}
}


void StreamProvider::analyzeOutputString(const string & output)
{
	if (output.empty() ) {
		SETFLAG(m_flags, m_flagNoOutput);
		RESETFLAG(m_flags, m_flagCloseAfterUse);
		return;
	} else
		RESETFLAG(m_flags, m_flagNoOutput);

	int i;
	for (i = output.size() - 1; i >= 0; --i)
		if (output[i] == '>' || output[i] == '|')
			break;

	++i;
	string type;
	string format;
	if (i == 0) {                                                                             // no specification >
		type = ">";
		format = output;
	} else {
		// >>file i=2, type: 0, size 2, format: start 2, size 6-2
		type = output.substr(0, i);
		format = output.substr(i, output.size() - i);
	}

	RESETFLAG(m_flags, m_flagCloseAfterUse);
	if (type == ">") {
		RESETFLAG(m_flags, m_flagAppend);
		SETFLAG(m_flags, m_flagCloseAfterUse);
	} else if (type == ">>") {
		SETFLAG(m_flags, m_flagAppend);
		RESETFLAG(m_flags, m_flagRealAppend);
	} else if (type == ">>>") {
		SETFLAG(m_flags, m_flagAppend);
		SETFLAG(m_flags, m_flagRealAppend);
	} else if (type == "|") {
		SETFLAG(m_flags, m_flagAppend);
		SETFLAG(m_flags, m_flagUseString);
	} else
		throw ValueError("Ostream types can only be one of '', '>', '>>', '>>>', '|'");

	if (format.empty()) {
		SETFLAG(m_flags, m_flagUseDefault);
		RESETFLAG(m_flags, m_flagCloseAfterUse);
	} else
		RESETFLAG(m_flags, m_flagUseDefault);

	DBG_DO(DBG_UTILITY, cerr << "Analyzed string is " << output << endl
		                     << "Filename is " << format << endl);

	m_filename = format;
}


void CloseOutput(const string & output)
{
	if (output.empty())
		ostreamManager().closeAll();
	else {
		DBG_ASSERT(ostreamManager().hasOstream(output), RuntimeError,
			"Output " + output + " does not exist or has already been closed.");
		ostreamManager().closeOstream(output);
	}
}


// Random number generator
RNG::RNG(const char * rng, unsigned long seed) : m_RNG(NULL)
{
	setRNG(rng, seed);
}


RNG::~RNG()
{
	// free current RNG
	gsl_rng_free(m_RNG);
}


#if  defined (_WIN32) || defined (__WIN32__)
// the following code is adapted from python os.urandom
//
typedef BOOL (WINAPI * CRYPTACQUIRECONTEXTA)(HCRYPTPROV * phProv, \
                                             LPCSTR pszContainer, LPCSTR pszProvider, DWORD dwProvType, \
                                             DWORD dwFlags);
typedef BOOL (WINAPI * CRYPTGENRANDOM)(HCRYPTPROV hProv, DWORD dwLen, \
                                       BYTE * pbBuffer);

unsigned long RNG::generateRandomSeed()
{
	unsigned long seed;

	HINSTANCE hAdvAPI32 = NULL;
	CRYPTACQUIRECONTEXTA pCryptAcquireContext = NULL;

	/* Obtain handle to the DLL containing CryptoAPI
	   This should not fail	*/
	hAdvAPI32 = GetModuleHandle("advapi32.dll");
	if (hAdvAPI32 == NULL) {
		DBG_WARNING(true, "advapi32.dll can not be loaded");
		return static_cast<unsigned long>(time(NULL));
	}

	/* Obtain pointers to the CryptoAPI functions
	   This will fail on some early versions of Win95 */
	pCryptAcquireContext = (CRYPTACQUIRECONTEXTA)GetProcAddress(
		hAdvAPI32, "CryptAcquireContextA");
	if (pCryptAcquireContext == NULL) {
		DBG_WARNING(true, "Failed to get process address of CryptAcquireContextA");
		return static_cast<unsigned long>(time(NULL));
	}

	CRYPTGENRANDOM pCryptGenRandom = (CRYPTGENRANDOM)GetProcAddress(
		hAdvAPI32, "CryptGenRandom");
	if (pCryptGenRandom == NULL) {
		DBG_WARNING(true, "Failed to get process address of CryptGenRandom");
		return static_cast<unsigned long>(time(NULL));
	}

	/* Acquire context */
	HCRYPTPROV hCryptProv = 0;
	if (!pCryptAcquireContext(&hCryptProv, NULL, NULL,
			PROV_RSA_FULL, CRYPT_VERIFYCONTEXT)) {
		DBG_WARNING(true, "Can not acquire context of CryptAcquireContextA");
		return static_cast<unsigned long>(time(NULL));
	}

	/* Get random data */
	if (!pCryptGenRandom(hCryptProv, sizeof(seed), (unsigned char *)&seed)) {
		DBG_WARNING(true, "Failed to get random number");
		return static_cast<unsigned long>(time(NULL));
	}

	DBG_DO(DBG_UTILITY, cerr << "Get random seed " << hex << seed << endl);
	return seed;
}


#else
unsigned long RNG::generateRandomSeed()
{
	// now, I need to work hard to get a good seed, considering
	// the use of clusters, and several jobs may be started at the same
	// time
	unsigned long seed;
	FILE * devrandom;

	if ((devrandom = fopen("/dev/urandom", "r")) != NULL) {
		UINT sz = fread(&seed, sizeof(seed), 1, devrandom);
		(void)sz; // suppress a warning.
		DBG_FAILIF(sz != 1, RuntimeError,
			"Incorrect bits of random digits are read from /dev/urandom");
		fclose(devrandom);
	} else if ((devrandom = fopen("/dev/random", "r")) != NULL) {
		UINT sz = fread(&seed, sizeof(seed), 1, devrandom);
		(void)sz; // suppress a warning.
		DBG_FAILIF(sz != 1, RuntimeError,
			"Incorrect bits of random digits are read from /dev/urandom");
		fclose(devrandom);
	} else {
		// this is not the best method, but I am out of ideas
		// of portable ways to add some other noises
		seed = static_cast<unsigned long>(time(NULL));
	}
	return seed;
}


#endif

// choose an random number generator.
void RNG::setRNG(const char * rng, unsigned long seed)
{
	const char * rng_name = rng;

	// if RNG name is not given, try GSL_RNG_TYPE
	if (rng_name == NULL)
		rng_name = getenv("GSL_RNG_TYPE");

	if (rng_name != NULL && rng_name[0] != '\0') {
		// locate the RNG
		const gsl_rng_type ** t, ** t0 = gsl_rng_types_setup();

		gsl_rng_default = 0;

		/* check GSL_RNG_TYPE against the names of all the generators */

		for (t = t0; *t != 0; t++) {
			// require that a RNG can generate full range of integer from 0 to the max of unsigned long int
			if (strcmp(rng_name, (*t)->name) == 0) {
				// free current RNG
				if (m_RNG != NULL)
					gsl_rng_free(m_RNG);

				m_RNG = gsl_rng_alloc(*t);

				DBG_ASSERT(gsl_rng_max(m_RNG) >= MaxRandomNumber && gsl_rng_min(m_RNG) == 0,
					ValueError, "You chosen random number generator can not generate full range of int.");
				break;
			}
		}

		if (*t == 0)
			throw SystemError("GSL_RNG_TYPE=" + toStr(rng_name)
				+ " not recognized or can not generate full range (0-2^32-1) of integers.\n c.f. AvailableRNGs()");
	} else {
		// free current RNG
		if (m_RNG != NULL)
			gsl_rng_free(m_RNG);

		m_RNG = gsl_rng_alloc(gsl_rng_mt19937);
	}

	// generate seed
	if (seed == 0)
		m_seed = generateRandomSeed();
	else
		m_seed = seed;

	// set seed
	gsl_rng_set(m_RNG, m_seed);
}


bool RNG::randBit()
{
	static WORDTYPE randbyte = 0;
	static UINT index = 0;

	if (index == 16)
		index = 0;

	if (index == 0)
		randbyte = randInt(0xFFFF);

	return (randbyte & (1UL << index++)) != 0;
}


double pvalChiSq(double chisq, unsigned int df)
{
	return 1 - gsl_cdf_chisq_P(chisq, df);
}


void chisqTest(const vector<vectoru> & table, double & chisq, double & chisq_p)
{
	UINT m = table.size();
	UINT n = table[0].size();
	vectoru rowSum(m, 0);
	vectoru colSum(n, 0);
	double N = 0;

	DBG_DO(DBG_STATOR, cerr << "ChiSq test with table\n" <<
		table[0] << endl << table[1] << endl <<
		(table.size() <= 2 ? vectoru() : table[2]) << endl);
	//
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			rowSum[i] += table[i][j];
			colSum[j] += table[i][j];
			N += table[i][j];
		}
	}
	chisq = 0;
	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < n; ++j)
			if (rowSum[i] > 0 && colSum[j] > 0)
				chisq += pow(table[i][j] - rowSum[i] * colSum[j] / N, 2)
				         / (rowSum[i] * colSum[j] / N);
	chisq_p = 1 - gsl_cdf_chisq_P(chisq, (m - 1) * (n - 1));
}


double armitageTrendTest(const vector<vectoru> & table, const vectorf & s)
{
	DBG_FAILIF(table.size() != 2 || table[0].size() != 3, ValueError,
		"Current Cochran-Armitage test can only handle 2 by 3 tables.");

	UINT n = table[0].size();

	DBG_FAILIF(s.size() != n, ValueError,
		"Weight for Cochran-Armitage test should have length 3");
	//
	DBG_DO(DBG_STATOR, cerr << "Armitage trend test with table\n" <<
		table[0] << endl << table[1] << endl);
	// formula is copied from HelixTree
	// www.goldenhelix.com/SNP_Variation/Manual/svs7/general_statistics.html
	//
	double N = 0;
	vectorf p1i(3);
	for (size_t i = 0; i < n; ++i) {
		p1i[i] = table[1][i] * 1.0 / (table[0][i] + table[1][i]);
		for (size_t j = 0; j < 2; ++j)
			N += table[j][i];
	}
	// s_bar
	double s_bar = 0;
	for (size_t i = 0; i < n; ++i)
		s_bar += (table[0][i] + table[1][i]) * s[i];
	s_bar /= N;
	// p_case
	double p_case = (table[1][0] + table[1][1] + table[1][2]) / N;
	// b
	double b;
	double b1 = 0;
	double b2 = 0;
	for (size_t i = 0; i < n; ++i) {
		b1 += (table[0][i] + table[1][i]) * (p1i[i] - p_case) * (s[i] - s_bar);
		b2 += (table[0][i] + table[1][i]) * (s[i] - s_bar) * (s[i] - s_bar);
	}
	b = b1 / b2;
	//
	double z2 = 0;
	for (size_t i = 0; i < n; ++i)
		z2 += (table[0][i] + table[1][i]) * (s[i] - s_bar) * (s[i] - s_bar);
	z2 = z2 * b * b / (p_case * (1 - p_case));
	//
	return 1 - gsl_cdf_chisq_P(z2, 1);
}


double hweTest(const vectoru & cnt)
{
	// Calculates exact two-sided hardy-weinberg p-value. Parameters
	// are number of genotypes, number of rare alleles observed and
	// number of heterozygotes observed.
	//
	// (c) 2003 Jan Wigginton, Goncalo Abecasis
	int obsAA = cnt[2];                                             // in this algorithm, AA is rare.
	int obsAB = cnt[1];
	int obsBB = cnt[0];

	int diplotypes = obsAA + obsAB + obsBB;
	int rare = (obsAA * 2) + obsAB;
	int hets = obsAB;


	//make sure "rare" allele is really the rare allele
	if (rare > diplotypes)
		rare = 2 * diplotypes - rare;

	//make sure numbers aren't screwy
	if (hets > rare)
		throw ValueError("HW test: " + toStr(hets) + "heterozygotes but only " + toStr(rare) + "rare alleles.");

	vectorf tailProbs(rare + 1, 0.0);

	//start at midpoint
	//all the casting is to make sure we don't overflow ints if there are 10's of 1000's of inds
	int mid = (int)((double)rare * (double)(2 * diplotypes - rare) / (double)(2 * diplotypes));

	//check to ensure that midpoint and rare alleles have same parity
	if (((rare & 1) ^ (mid & 1)) != 0) {
		mid++;
	}
	int het = mid;
	int hom_r = (rare - mid) / 2;
	int hom_c = diplotypes - het - hom_r;

	//Calculate probability for each possible observed heterozygote
	//count up to a scaling constant, to avoid underflow and overflow
	tailProbs[mid] = 1.0;
	double sum = tailProbs[mid];
	for (het = mid; het > 1; het -= 2) {
		tailProbs[het - 2] = (tailProbs[het] * het * (het - 1.0)) / (4.0 * (hom_r + 1.0) * (hom_c + 1.0));
		sum += tailProbs[het - 2];
		//2 fewer hets for next iteration -> add one rare and one common homozygote
		hom_r++;
		hom_c++;
	}

	het = mid;
	hom_r = (rare - mid) / 2;
	hom_c = diplotypes - het - hom_r;
	for (het = mid; het <= rare - 2; het += 2) {
		tailProbs[het + 2] = (tailProbs[het] * 4.0 * hom_r * hom_c) / ((het + 2.0) * (het + 1.0));
		sum += tailProbs[het + 2];
		//2 more hets for next iteration -> subtract one rare and one common homozygote
		hom_r--;
		hom_c--;
	}

	for (size_t z = 0; z < tailProbs.size(); z++)
		tailProbs[z] /= sum;

	double top = tailProbs[hets];
	for (int i = hets + 1; i <= rare; i++)
		top += tailProbs[i];

	double otherSide = tailProbs[hets];
	for (int i = hets - 1; i >= 0; i--)
		otherSide += tailProbs[i];

	if (top > 0.5 && otherSide > 0.5) {
		return 1.0;
	} else {
		if (top < otherSide)
			return top * 2;
		else
			return otherSide * 2;
	}

	return 1.0;
}


void weightedSampler::set(const vectorf & weight)
{
	m_N = weight.size();

	if (m_N == 0)
		return;

	if (m_N == 1) {
		m_fixed = true;
		m_fixedValue = 0;
		return;
	}

	m_fixed = true;
	int prevIndex = -1;
	for (size_t i = 0; i < weight.size(); ++i) {
		if (weight[i] != 0) {
			if (prevIndex == -1) {
				m_fixedValue = i;
				prevIndex = i;
			} else { // two non-zero index, not fixed.
				m_fixed = false;
				break;
			}
		}
	}
	if (m_fixed)
		return;

	// sum of weight
	double w = accumulate(weight.begin(), weight.end(), 0.0);

	DBG_FAILIF(fcmp_le(w, 0), ValueError, "Sum of weight is <= 0.");

	w = m_N / w;

	// initialize p with N*p0,...N*p_k-1
	m_q.resize(m_N);
	for (size_t i = 0; i < m_N; ++i)
		m_q[i] = weight[i] * w;
	// initialize Y with values
	m_a.resize(m_N);
	for (size_t i = 0; i < m_N; ++i)
		m_a[i] = i;
	// use two sets H and L
	// for efficiency purpose, use a single vector.
	ULONG * HL = new ULONG[m_N];
	ULONG * L = HL;
	ULONG * H = HL + m_N - 1;                                         // point to the end.

	for (size_t i = 0; i < m_N; ++i) {
		if (m_q[i] > 1)
			*H-- = i;
		else
			*L++ = i;
	}

	//
	ULONG j, k;
	while (L != HL && H != HL + m_N - 1) {
		j = *(L - 1);
		k = *(H + 1);
		m_a[j] = k;
		m_q[k] += m_q[j] - 1;

		L--;                                                                            // remove j from L
		if (m_q[k] < 1.) {
			*L++ = k;                                                                   // add k to L
			++H;                                                                        // remove k from H
		}
	}
	delete[] HL;
}


// this is used for BernulliTrials and copyGenotype
WORDTYPE g_bitMask[WORDBIT];

BernulliTrials::BernulliTrials(RNG & rng)
	: m_RNG(&rng), m_N(0), m_prob(0), m_table(0), m_pointer(0),
	m_cur(npos)
{
}


BernulliTrials::BernulliTrials(RNG & rng, const vectorf & prob, ULONG trials)
	: m_RNG(&rng), m_N(trials), m_prob(prob), m_table(prob.size()), m_pointer(prob.size()),
	m_cur(npos)
{
	DBG_FAILIF(trials <= 0, ValueError, "trial number can not be zero.");
	DBG_FAILIF(prob.empty(), ValueError, "probability table can not be empty.");

	// initialize the table
	for (size_t i = 0; i < probSize(); ++i) {
		DBG_FAILIF(m_prob[i] < 0 || m_prob[i] > 1, ValueError,
			"Probability for a Bernulli trail should be between 0 and 1 (value "
			+ toStr(m_prob[i]) + " at index " + toStr(i) + ")");
		m_table[i].resize(trials);
		m_pointer[i] = BITPTR(m_table[i].begin());
	}
}


//
BernulliTrials::~BernulliTrials()
{
}


void BernulliTrials::setParameter(const vectorf & prob, ULONG trials)
{
	m_N = trials;
	m_prob = prob;
	m_table.resize(m_prob.size());
	m_pointer.resize(m_prob.size());
	m_cur = npos;                                                             // will trigger doTrial.

	DBG_FAILIF(trials <= 0, ValueError, "trial number can not be zero.");
	DBG_FAILIF(prob.empty(), ValueError, "probability table can not be empty.");

	for (size_t i = 0; i < probSize(); ++i) {
		DBG_FAILIF(m_prob[i] < 0 || m_prob[i] > 1, ValueError,
			"Probability for a Bernulli trail should be between 0 and 1 (value "
			+ toStr(m_prob[i]) + " at index " + toStr(i) + ")");
		m_table[i].resize(trials);
		m_pointer[i] = BITPTR(m_table[i].begin());
	}
}


// utility function.
void BernulliTrials::setAll(size_t idx, bool v)
{
	WORDTYPE * ptr = m_pointer[idx];

	DBG_ASSERT(BITOFF(m_table[idx].begin()) == 0, SystemError, "Start of a vector<bool> is not 0");
	DBG_ASSERT(BITPTR(m_table[idx].begin()) == m_pointer[idx],
		SystemError, "Pointers mismatch");

	size_t blk = m_N / WORDBIT;
	size_t rest = m_N - blk * WORDBIT;
	if (v) {
		// set all to 1
		for (size_t i = 0; i < blk; ++i)
			*ptr++ = ~WORDTYPE(0UL);
		if (rest > 0) {
			*ptr |= g_bitMask[rest];
			// upper to 0
			*ptr &= g_bitMask[rest];
		}
	} else {
		for (size_t i = 0; i < blk; ++i)
			*ptr++ = 0UL;
		if (rest > 0)
			*ptr = 0; //~g_bitMask[rest];
	}
}


#define setBit(ptr, i)    (*((ptr) + (i) / WORDBIT) |= 1UL << ((i) - ((i) / WORDBIT) * WORDBIT))
#define unsetBit(ptr, i)  (*((ptr) + (i) / WORDBIT) &= ~(1UL << ((i) - ((i) / WORDBIT) * WORDBIT)))
// use a != 0 to avoid compiler warning
#define getBit(ptr, i)    ((*((ptr) + (i) / WORDBIT) & (1UL << ((i) - ((i) / WORDBIT) * WORDBIT))) != 0)

void BernulliTrials::doTrial()
{
    DBG_ASSERT(m_N != 0, ValueError, "number of trials should be positive");

    DBG_DO(DBG_UTILITY, cerr << "n=" << m_N << " doTrial, cur trial: " << m_cur << endl);

    // for each column
    for (size_t cl = 0, clEnd = probSize(); cl < clEnd; ++cl) {
        WORDTYPE * ptr = m_pointer[cl];
        double prob = m_prob[cl];
        if (prob == 0.) {
            setAll(cl, false);
		} else if (prob == 0.5) {                                 // random 0,1 bit, this will be quicker
            // set to 0..
            setAll(cl, false);
            // treat a randInt as random bits and set them directly.
            // I.e., we will call 1/16 or 1/32 times of rng for this specifal case.
            // first several blocks
            // WORDTYPE * ptr = BITPTR(succ.begin());
            WORDTYPE tmp;
            size_t blk = m_N / WORDBIT;
            size_t rest = m_N - blk * WORDBIT;
            for (size_t i = 0; i < blk; ++i) {
                // even if the block size is large (I can not set it to int16_t)
                // I only take the last 16 bit of a rng
                // for the quality of random bits.
                *ptr = 0;
                for (size_t b = 0; b < WORDBIT / 16; ++b) {
                    // blocks[i] = static_cast<int16_t>(GetRNG().randGet());
                    tmp = GetRNG().randInt(0xFFFF);
                    *ptr |= (0xFFFF & tmp) << (b * 16);
				}
                ptr++;
			}
            // last block
            if (rest != 0) {
                size_t b = 0;
                for (b = 0; b < rest / 16; ++b) {
                    tmp = GetRNG().randInt(0xFFFF);
                    *ptr |= (0xFFFF & tmp) << (b * 16);
				}
                // last bits
                rest -= b * 16;
                if (rest != 0) {
                    tmp = GetRNG().randInt(0xFFFF);
                    *ptr |= (g_bitMask[rest] & tmp) << b * 16;
				}
			}
		}
        // algorithm i Sheldon Ross' book simulation (4ed), page 54
        else if (prob < 0.5) {
            // set all to 0, then set some to 1
            setAll(cl, false);
            // it may make sense to limit the use of this method to low p,
            UINT i = 0;
            while (true) {
                // i moves at least one.
                i += m_RNG->randGeometric(prob);
                if (i <= m_N)
					// succ[i-1] = true;
					setBit(ptr, i - 1);
                else
					break;
			}
		} else if (prob == 1.) {
            setAll(cl, true);
		} else {                                                                  // 1 > m_proc[cl] > 0.5
            // set all to 1, and then unset some.
            setAll(cl, true);
            // it may make sense to limit the use of this method to low p,
            UINT i = 0;
            prob = 1. - prob;
            while (true) {
                i += m_RNG->randGeometric(prob);
                if (i <= m_N)
					// succ[i-1] = false;
					unsetBit(ptr, i - 1);
                else
					break;
			}
		}
	}
    m_cur = 0;
}


UINT BernulliTrials::curTrial()
{
    DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
    return m_cur;
}


// get a trial corresponding to m_prob.
void BernulliTrials::trial()
{
    if (m_cur == npos || m_cur == m_N - 1)  // reach the last trial
		doTrial();
    else
		m_cur++;
    DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
}


bool BernulliTrials::trialSucc(size_t idx) const
{
    DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
    return getBit(m_pointer[idx], m_cur);
}


bool BernulliTrials::trialSucc(size_t idx, size_t cur) const
{
    return getBit(m_pointer[idx], cur);
}


size_t BernulliTrials::probFirstSucc() const
{
    DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
    size_t i = 0;
    const size_t sz = probSize();
    while (i < sz && !getBit(m_pointer[i], m_cur))
		++i;
    return i >= sz ? npos : i;
}


size_t BernulliTrials::probNextSucc(size_t pos) const
{
    DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
    const size_t sz = probSize();
    if (pos >= (sz - 1) || sz == 0)
		return npos;

    ++pos;
    while (pos < sz && !getBit(m_pointer[pos], m_cur))
		++pos;
    return pos >= sz ? npos : pos;
}


size_t BernulliTrials::trialFirstSucc(size_t idx) const
{
    size_t blk = m_N / WORDBIT;
    WORDTYPE * ptr = m_pointer[idx];

    size_t i = 0;

    while (i < blk && *ptr++ == 0)
		++i;

    if (i < blk) {                                                                          // not at the last blk
        return i * WORDBIT + lowest_bit(*(ptr - 1));
	} else {                                                                                // last block?
        size_t rest = m_N - blk * WORDBIT;
        size_t tmp = *ptr & g_bitMask[rest];
        if (tmp == 0)
			return npos;
        else
			return blk * WORDBIT + lowest_bit(tmp);
	}
}


size_t BernulliTrials::trialNextSucc(size_t idx, size_t pos) const
{
    const BitSet & bs = m_table[idx];

    if (pos >= (m_N - 1) || m_N == 0)
		return npos;

    ++pos;

    // first block
    BitSet::const_iterator it = bs.begin() + pos;
    WORDTYPE * ptr = BITPTR(it);
    int offset = BITOFF(it);
    size_t i = ptr - BITPTR(bs.begin());

    // mask out bits before pos
    WORDTYPE tmp = *ptr & ~g_bitMask[offset];

    size_t blk = m_N / WORDBIT;
    if (tmp != 0)
		return i * WORDBIT + lowest_bit(tmp);
    else if (blk == i)
		// if i is the last block, return no.
		return npos;

    // now, go from next block
    ++ptr;
    ++i;
    while (i < blk && *ptr++ == 0)
		++i;

    if (i < blk)                                                                            // not at the last blk
		return i * WORDBIT + lowest_bit(*(ptr - 1));
    else {                                                                                  // last block?
        size_t rest = m_N - blk * WORDBIT;
        // mask out bits after rest
        tmp = *ptr & g_bitMask[rest];
        if (tmp == 0)
			return npos;
        else
			return blk * WORDBIT + lowest_bit(tmp);
	}
}


void BernulliTrials::setTrialSucc(size_t idx, bool succ)
{
    DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
    if (succ)
		setBit(m_pointer[idx], m_cur);
    else
		unsetBit(m_pointer[idx], m_cur);
}


double BernulliTrials::trialSuccRate(UINT index) const
{
    // efficiency is not considered here
    size_t count = 0;

    for (size_t i = 0; i < trialSize(); ++i)
		if (getBit(m_pointer[index], i))
			count++;
    return count / static_cast<double>(m_table[index].size());
}


double BernulliTrials::probSuccRate() const
{
    DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
    UINT count = 0;
    for (size_t cl = 0, clEnd = probSize(); cl < clEnd; ++cl)
		count += getBit(m_pointer[cl], m_cur) ? 1 : 0;
    return count / static_cast<double>(probSize());
}


#undef setBit
#undef unsetBit
#undef getBit

// random number generator. a global variable.
// there might be multiple RNG later.
RNG g_RNG;

// return the global RNG
RNG & GetRNG()
{
    return g_RNG;
}


// set global rng
// this is temporary since GetRNG() might not exist in the future
void SetRNG(const string r, unsigned long seed)
{
    GetRNG().setRNG(r.c_str(), seed);
}


// list all available RNG.
vectorstr AvailableRNGs()
{
    vectorstr list;

    const gsl_rng_type ** t, ** t0;
    gsl_rng * rng;

    t0 = gsl_rng_types_setup();

    for (t = t0; *t != 0; t++) {
        rng = gsl_rng_alloc(*t);
        if (gsl_rng_min(rng) == 0 && gsl_rng_max(rng) >= MaxRandomNumber)
			list.push_back((*t)->name);
        gsl_rng_free(rng);
	}
    return list;
}


// Global debug and initialization related functions

void gsl_error_handler(const char * reason, const char *,
                       int, int gsl_errno)
{
    throw SystemError("GSL Error " + toStr(gsl_errno) + ":\t"
		+ reason);
}


// create a stream buf that print to python sys.stdout
// cerr will be redirected to this buf to really output
// to python console.
class PythonStdBuf : public streambuf
{
public:
	enum PythonBufType {
		StdOut = 1,
		StdErr = 2
	};
public:
    PythonStdBuf(PythonBufType type) : m_type(type)
    {
	}


protected:
    int overflow(int c)
    {
        // write out current buffer
        if (pbase() != pptr() ) {
            // the end of string might not be \0
            char_type * endPtr = pptr();
            char_type endChar = *endPtr;
            *endPtr = '\0';
            //          int len = pptr() - pbase();
            //          char * str = new char[len+1];
            //          strncpy(str, pbase(), len);
            //          str[len] = '\0';
			if (m_type == StdOut)
	            PySys_WriteStdout("%s", pbase());
			else
				PySys_WriteStderr("%s", pbase());
            // put original end character back, whatever it is.
            *endPtr = endChar;
            //          delete[] str;
            setp(pbase(), epptr());
		}
        // push character c in
        if (c != EOF) {
            // unbuffered, write out this character, do not put into buffer
            if (pbase() == epptr() ) {
				if (m_type == StdOut)
					PySys_WriteStdout("%c", c);
				else
					PySys_WriteStderr("%c", c);
			} else
				sputc(c);
		}
        return 0;
	}

protected:
	PythonBufType m_type;
};

// create a null stream buf that discard everything
/// CPPONLY
class NullStreamBuf : public streambuf
{
public:
    NullStreamBuf()
    {
	}


protected:
    int overflow(int)
    {
        return 0;
	}


};

/// CPPONLY
PythonStdBuf g_pythonStdoutBuf(PythonStdBuf::StdOut);

/// CPPONLY
PythonStdBuf g_pythonStderrBuf(PythonStdBuf::StdErr);

/// CPPONLY
NullStreamBuf g_nullStreamBuf;

// null stream
ostream g_cnull(&g_nullStreamBuf);

// return null stream
ostream & cnull()
{
    return g_cnull;
}


PyObject * ModuleInfo()
{
    // output a dictionary with many keys
    PyObject * dict = PyDict_New();

    // these macros will be passed from commandline, if not, use the default
#ifndef SIMUPOP_REV
#  define REVISION "9999"
#else
    // make passed macro to a real string
#  define REVISION MacroQuote(SIMUPOP_REV)
#endif

    // Revision
    const char * rev = REVISION;
    int num;
    // XX:XX or XX:XXM, or XX
    if (sscanf(rev, "%*d:%d", &num) != 1 && sscanf(rev, "%d", &num) != 1)
		num = 9999;
    PyDict_SetItem(dict, PyString_FromString("revision"), PyInt_FromLong(num));

    // Version
#ifndef SIMUPOP_VER
    const char * ver = "snapshot";
#else
    // convert name to a string
    const char * ver = MacroQuote(SIMUPOP_VER);
#endif
    PyDict_SetItem(dict, PyString_FromString("version"), PyString_FromString(ver));

    // optimized
#ifdef OPTIMIZED
    Py_INCREF(Py_True);
    PyDict_SetItem(dict, PyString_FromString("optimized"), Py_True);
#else
    Py_INCREF(Py_False);
    PyDict_SetItem(dict, PyString_FromString("optimized"), Py_False);
#endif

    // AlleleType
#ifdef LONGALLELE
    PyDict_SetItem(dict, PyString_FromString("alleleType"), PyString_FromString("long"));
#else
#  ifdef BINARYALLELE
    PyDict_SetItem(dict, PyString_FromString("alleleType"), PyString_FromString("binary"));
#  else
    PyDict_SetItem(dict, PyString_FromString("alleleType"), PyString_FromString("short"));
#  endif
#endif

#ifndef COMPILER
#  ifdef __GNUC__
#    define COMPILER "[GCC " __VERSION__ "]"
#  endif
#endif                                                                            /* !COMPILER */

#ifndef COMPILER
#  ifdef __cplusplus
#    define COMPILER "[C++]"
#  else
#    define COMPILER "[C]"
#  endif
#endif

#ifndef PLATFORM
#  define PLATFORM ""
#endif

    // compiler
    PyDict_SetItem(dict, PyString_FromString("compiler"), PyString_FromString(COMPILER));

    // date
    PyDict_SetItem(dict, PyString_FromString("date"), PyString_FromString(__DATE__));

    // version of python
    PyDict_SetItem(dict, PyString_FromString("python"), PyString_FromString(PY_VERSION));

    // platform
    PyDict_SetItem(dict, PyString_FromString("platform"), PyString_FromString(PLATFORM));

    // maxAllele
    PyDict_SetItem(dict, PyString_FromString("maxAllele"), PyInt_FromLong(ModuleMaxAllele));

    // limits
    PyDict_SetItem(dict, PyString_FromString("maxIndex"), PyLong_FromUnsignedLong(MaxIndexSize));

	// debug (code)
	PyObject * codes = PyList_New(0);
	for (size_t i = 0; i < DBG_CODE_LENGTH; ++i) {
		if (g_dbgCode[i]) {
			PyList_Append(codes, PyString_FromString(g_debugCodes[i]));
		}
	}
	PyDict_SetItem(dict, PyString_FromString("debug"), codes);

    //
    return dict;
}


#ifdef BINARYALLELE

// define a good way to copy long genotype sequence
void copyGenotype(GenoIterator fr, GenoIterator to, size_t n)
{
    DBG_ASSERT(BITOFF(fr) < WORDBIT, SystemError,
		"Your vector<bool> implementation is different...");

    WORDTYPE * fr_p = BITPTR(fr);
    WORDTYPE * to_p = BITPTR(to);
    unsigned int fr_off = BITOFF(fr);
    unsigned int to_off = BITOFF(to);

    // if offset is different, can not copy in block.
    if (n < WORDBIT) {
        for (size_t i = 0; i < n; ++i) {
            // set bit according to from bit
            if (*fr_p & (1UL << fr_off))
				*to_p |= (1UL << to_off);
            else
				*to_p &= ~(1UL << to_off);
            // next location
            if (fr_off++ == WORDBIT - 1) {
                fr_off = 0;
                ++fr_p;
			}
            if (to_off++ == WORDBIT - 1) {
                to_off = 0;
                ++to_p;
			}
		}
	} else if (fr_off == to_off) {
        // copy first block, fr_off + 1 bits
        // ABCDxxxx for off=4
        // what I am doing is
        // mask[4] = xxxx1111 <==== NOTICE mask[4] has four 1s
        //
        //  off = 4,
        //  from:  ABCDxxxx
        //  to:    EFGHxxxx
        // from & ~mask =  ABCD0000
        // to &    mask  = 0000xxxx
        // (from & mask) | (to & ~mask) = ABCDxxxx
        WORDTYPE mask = g_bitMask[fr_off];
        *to_p = (*fr_p & ~mask) | (*to_p & mask);

        size_t rest = n - (WORDBIT - fr_off);
        size_t blks = rest / WORDBIT;
        for (size_t i = 0; i < blks; ++i)
			*++to_p = *++fr_p;
        // the rest of the bits?
        rest -= blks * WORDBIT;
        if (rest != 0) {
            // rest = 3
            // from:   xxxxxABC
            // to:     xxxxxCDE
            // mask:   00000111
            //    &:   00000ABC
            //    & :  xxxxx000
            //    |:   xxxxxABC
            to_p++;
            fr_p++;
            // mask has top rest zeros.
            mask = g_bitMask[rest];
            *to_p = (*fr_p & mask) | (*to_p & ~mask);
		}
	} else if (fr_off < to_off) {
        WORDTYPE maskFrom = g_bitMask[fr_off];
        WORDTYPE maskTo = g_bitMask[to_off];
        WORDTYPE maskFrom1;
        size_t shift;

        shift = to_off - fr_off;
        // fr_off=5, to_off=7, shift=3
        // from:    ABCxxxxx,  maskFrom: 00011111
        // to:      Dxxxxxxx,  maskTo:   01111111
        //   (from & ~maskFrom) = ABC00000
        //   <<                   C0000000
        //   to & maskTo        = 0xxxxxxx
        //   |                  = Cxxxxxxx
        *to_p = ((*fr_p & ~maskFrom) << shift) |
                (*to_p & maskTo);
        // now. for other block bits
        // to other bits
        size_t rest = n - (WORDBIT - to_off);
        size_t blks = rest / WORDBIT;
        //
        //  already copied to_off+1 bits
        // from:   ABxxxxxx, maskFrom:    00111111
        // from1:  xxCDEFGH, maskFrom:    00111111
        // from & ~maskFrom            =  AB000000
        // >>                          =  000000AB
        // from1 & maskFrom            =  00CDEFGH
        // <<                          =  CDEFGHxx
        // |                           =  CDEFGHAB
        maskFrom = g_bitMask[WORDBIT - shift];
        for (size_t i = 0; i < blks; ++i) {
            to_p++;
            *to_p = ((*fr_p & ~maskFrom) >> (WORDBIT - shift)) |
                    ( (*(fr_p + 1) & maskFrom) << shift);
            fr_p++;
		}
        // the rest of the bits?
        rest -= blks * WORDBIT;
        if (rest != 0) {
            to_p++;
            if (rest <= shift) {
                // rest = 2, shift = 5
                // from:  ABCEDxxx, maskfrom: 00000111
                // to:    xxxxxxAB, maskto:   00000011
                // & ~                        ABCDE000
                // >>                         000ABCDE
                // & naskTo                   000000DE
                // & ~                        xxxxxx00
                // |                          xxxxxxDE
                maskFrom = g_bitMask[WORDBIT - shift];
                maskTo = g_bitMask[rest];
                *to_p = (((*(fr_p) & ~maskFrom) >> (WORDBIT - shift)) & maskTo)
                        | (*to_p & ~maskTo);
			} else {
                // rest = 5, shift = 2
                // from:     ABxxxxxx, maskFrom: 00111111
                // from:     xxxxxCDE, maskFrom1:00000111
                // to:       xxxCDEAB, maskTo:   00011111
                //   from & ~ mask  =  AB000000
                //     >>           =  000000AB .
                //   from 1 & mask1 =  00000CDE
                //     <<           =  000CDE00 .
                //  to * ~mask      =  xxx00000
                //  |               =  xxxCDEAB
                maskFrom = g_bitMask[WORDBIT - shift];
                maskFrom1 = g_bitMask[rest - shift];
                maskTo = g_bitMask[rest];

                *to_p = ((*(fr_p) & ~maskFrom) >> (WORDBIT - shift)) |
                        ((*(fr_p + 1) & maskFrom1) << shift) |
                        (*to_p & ~maskTo);
			}
		}
	} else {                                                                          // fr_off > to_off
        size_t shift = fr_off - to_off;
        WORDTYPE maskFrom = g_bitMask[fr_off];
        WORDTYPE maskFrom1 = g_bitMask[shift];
        WORDTYPE maskTo = g_bitMask[to_off];
        // from:   ABCxxxxx, maskFrom: 00011111
        // from1:  xxxxxxDE, maskFrom1:00000011
        // to:     DEABCxxx, maskTo:   00011111
        //
        *to_p = ((*fr_p & ~maskFrom) >> shift) |
                ( (*(fr_p + 1) & maskFrom1) << (WORDBIT - shift)) |
                (*to_p & maskTo);
        to_p++;
        fr_p++;
        //
        // to other bits
        size_t rest = n - (WORDBIT - to_off);
        size_t blks = rest / WORDBIT;
        //
        // already copied shift bits
        // from:   ABCDEFxx, maskFrom:   00000011
        // from1:  xxxxxxGH, maskFrom:   00000011
        // to:     GHABCDEF
        maskFrom = g_bitMask[shift];
        for (size_t i = 0; i < blks; ++i) {
            *to_p = ((*fr_p & ~maskFrom) >> shift) |
                    ( (*(fr_p + 1) & maskFrom) << (WORDBIT - shift));
            fr_p++;
            to_p++;
		}
        // the rest of the bits
        rest -= blks * WORDBIT;
        if (rest != 0) {
            if (rest < WORDBIT - shift) {
                // rest = 2, shift = 3,
                // from:     ABCDExxx, maskFrom: 00000111
                // to:       xxxxxxDE, maskTo:   00000011
                maskFrom = g_bitMask[shift];
                maskTo = g_bitMask[rest];
                *to_p = (((*fr_p & ~maskFrom) >> shift) & maskTo) |
                        (*to_p & ~maskTo) ;
			} else {
                // rest = 5, shift = 6
                // from:   ABxxxxxx, maskFrom: 00111111
                // from1:  xxxxxCDE, maskFrom1:00000111
                // rest:   xxxCDEAB, maskTo:   00011111
                maskFrom = g_bitMask[shift];
                maskFrom1 = g_bitMask[rest - (WORDBIT - shift)];
                maskTo = g_bitMask[rest];
                *to_p = ((*fr_p & ~maskFrom) >> shift) |
                        ((*(fr_p + 1) & maskFrom1) << (WORDBIT - shift)) |
                        (*to_p & ~maskTo);
			}
		}
	}
#  ifndef OPTIMIZED
    if (debug(DBG_UTILITY)) {
        if (vectora(fr, fr + n) != vectora(to, to + n)) {
            cerr << "Copy from " << vectora(fr, fr + n)
                 << " to " << vectora(to, to + n) << " failed " << endl;
            cerr << "Offsets are " << BITOFF(fr) << " and " << BITOFF(to) << endl;
		}
	}
#  endif
}


void clearGenotype(GenoIterator to, size_t n)
{
    WORDTYPE * to_p = BITPTR(to);
    unsigned int to_off = BITOFF(to);

    // This can be made more efficient.
    for (size_t i = 0; i < n; ++i) {
        // set bit according to from bit
        *to_p &= ~(1UL << to_off);
        if (to_off++ == WORDBIT - 1) {
            while (i + WORDBIT < n) {
                *++to_p = 0;
                i += WORDBIT;
			}
            ++to_p;
            to_off = 0;
		}
	}
}


#endif

/* This file is used to initialize simuPOP when being load into
   python. The swig interface file will has a init% % entry to
   include this file. */
bool initialize()
{
    // tie python stdout to cerr
    std::cout.rdbuf(&g_pythonStdoutBuf);
    std::cerr.rdbuf(&g_pythonStderrBuf);

#if __WORDSIZE == 32
    DBG_ASSERT(WORDBIT == 32, SystemError,
		"We are assuming 32 bit word size for this system but we are wrong."
		"Please report this problem to the simuPOP mailinglist.");
#endif
    // for example, if WORDBIT is 8 (most likely 32), we define
    // 0x0, 0x1, 0x3, 0x7, 0xF, 0x1F, 0x3F, 0x7F, 0xFF
    for (size_t i = 0; i < WORDBIT; ++i) {
        // g_bitMask[i] is the number of 1 count from right.
        for (size_t j = 0; j < i; ++j)
			g_bitMask[i] |= (1UL << j);
	}

    // give at most 100 ref count warnings.
#ifdef Py_REF_DEBUG
    g_refWarningCount = 100;
#endif

    // SIMUPOP_MODULE is passed as name, but we need it to be quoted.
    // Note that under gcc, I could pass the macro from command line
    // using \" \" but this trick does not work under VC.
    // the following process is safer.
#define SimuPOP_Module_Name "##SIMUPOP_MODULE##"

    // set global dictionary/variable
    PyObject * mm = PyImport_AddModule(SimuPOP_Module_Name);
    g_module_vars = SharedVariables(PyModule_GetDict(mm), false);

    // main dictionary
    mm = PyImport_AddModule("__main__");
    g_main_vars = SharedVariables(PyModule_GetDict(mm), false);

    // get population and individual type pointer
    g_swigPopType = SWIG_TypeQuery(PopSWIGType);
    g_swigindividual = SWIG_TypeQuery(IndSWIGType);
    //
    // g_swigOperator = SWIG_TypeQuery(OperatorSWIGType);
    if (g_swigPopType == NULL || g_swigindividual == NULL)
		throw SystemError("Can not get population and individual type pointer, your SWIG version may be run.");

    // load carray function and type
    if (initCustomizedTypes() < 0)
		throw SystemError("Failed to initialize carray and defdict types");

    // set gsl error handler
    gsl_set_error_handler(&gsl_error_handler);

#ifndef OPTIMIZED
#  ifdef BINARYALLELE
    // binary level genotype copy is compiler dependent and may
    // fail on some systems. Such a test will make sure the binary
    testCopyGenotype();
#  endif
#endif
    return true;
}


bool intList::match(UINT rep, const vector<bool> & activeRep)
{
    if (m_elems.empty())
		return m_allAvail;
    vectori::iterator it = m_elems.begin();
    vectori::iterator it_end = m_elems.end();
    for (; it != it_end; ++it) {
        // positive index is easy
        if (*it >= 0) {
            if (static_cast<UINT>(*it) == rep)
				return true;
            else
				continue;
		}
        // negative index
        DBG_ASSERT(!activeRep.empty() && activeRep[rep], SystemError,
			"Check is avtive should only be done for active replicates");
        // check the simple and most used case
        if (*it == -1 && activeRep.back() && rep + 1 == activeRep.size())
			return true;
        // find what exactly an negative index refer to
        int cnt = -*it;
        int curRep = activeRep.size() - 1;
        for (; curRep >= 0; --curRep) {
            if (activeRep[curRep])
				--cnt;
            if (cnt == 0)
				break;
		}
        if (cnt != 0)
			continue;
        if (curRep == static_cast<int>(rep))
			return true;
	}
    return false;
}


#ifndef OPTIMIZED

#  ifdef BINARYALLELE
void testCopyGenotype()
{
    vectora from(1000);
    vectora to(1000);

    for (size_t i = 0; i < 100; ++i) {
        for (size_t j = 0; j < 1000; ++j) {
            // use != 0 to reduce compiler warning
            from[j] = GetRNG().randInt(2) != 0;
            to[j] = 0;
		}
        size_t from_idx = GetRNG().randInt(300);
        size_t to_idx = GetRNG().randInt(300);
        if (from_idx > to_idx)
			continue;
        size_t length = GetRNG().randInt(500);
        copyGenotype(from.begin() + from_idx,
			to.begin() + to_idx, length);
        if (vectora(from.begin() + from_idx, from.begin() + from_idx + length) !=
            vectora(to.begin() + to_idx, to.begin() + to_idx + length)) {
            cerr << "Copying: " << vectora(from.begin() + from_idx, from.begin() + from_idx + length) << '\n'
                 << "Obtain:  " << vectora(to.begin() + to_idx, to.begin() + to_idx + length) << '\n'
                 << "Index From: " << from_idx << " to: " << to_idx << " length: " << length << endl;
            // the error message can not be shown
            throw SystemError("Allele copy test for your system fails.\n"
				              "Please email simuPOP mailing list with detailed os and compiler information");
		}
	}
}


#  endif
#endif

}
