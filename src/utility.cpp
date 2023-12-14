/**
 *  $File: utility.cpp $
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
   \brief implementation of some utility functions.
 */

#include "utility.h"
#include <cstdlib>
#include "time.h"

#include <bitset>

#include "gsl/gsl_machine.h"
#ifdef _OPENMP
#  include "omp.h"
#endif

#include "population.h"

#include <sstream>
using std::ostringstream;
using std::stringstream;
using std::ostringstream;
using std::hex;

#include <fstream>
using std::fstream;
using std::ifstream;
using std::ofstream;

#include "boost_pch.hpp"

// for data type lociList
#include "genoStru.h"

// for PySys_WriteStdout and python expressions
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#include "swigpyrun.h"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wmissing-field-initializers"

// The PyString_Check should be defined after swingpyrun.h because
// in swingpyrun.h, the PyString_Check is defined to PyBytes_Check
#if PY_VERSION_HEX >= 0x03000000

#  define PyString_Check PyUnicode_Check
#  define PyString_FromStringAndSize PyUnicode_FromStringAndSize
#  define PyString_FromString PyUnicode_FromString
#  define PyString_Concat PyUnicode_Concat
#  define PyString_ConcatAndDel PyUnicode_ConcatAndDel

#  define PyInt_Check(x) PyLong_Check(x)
#  define PyInt_AsLong(x) PyLong_AsLong(x)
#  define PyInt_AS_LONG(x) PyLong_AsLong(x)
#  define PyInt_FromLong(x) PyLong_FromLong(x)
#  define PyInt_FromSize_t(x) PyLong_FromSize_t(x)
#  define PyInt_FromString PyLong_FromString
#  define PyInt_Type PyLong_Type
#  define PyFloat_FromString(x, y)  PyFloat_FromString(x)
#  define PyString_Type PyUnicode_Type

#endif

// compile and eval enables compiling string to byte code
#include "compile.h"

#if PY_VERSION_HEX < 0x030b0000
#include "eval.h"
#endif

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

#include "boost/dynamic_bitset/detail/lowest_bit.hpp"
using boost::detail::lowest_bit;


using boost::regex;
using boost::regex_match;
using boost::cmatch;

//macros to provide a portable way to make macro-passed string to C++ string
#define MacroQuote_(x) # x
#define MacroQuote(x) MacroQuote_(x)

// these functions are defined in customizedTypes.c which is included
// in simuPOP_wrap.cpp

extern "C" PyObject * newcarrayobject(GenoIterator begin, GenoIterator end);

extern "C" PyObject * newcarrayobject_lineage(LineageIterator begin, LineageIterator end);

extern "C" PyObject * PyDefDict_New();

extern "C" bool is_defdict(PyTypeObject * type);

extern "C" int initCustomizedTypes(PyObject * m);

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

// global constant variables
#ifdef MUTANTALLELE
const Allele simuPOP::vectorm::zero_ = 0;
#endif
const unsigned long ModuleMaxAllele = std::numeric_limits<Allele>::max();
const unsigned long MaxRandomNumber = std::numeric_limits<int32_t>::max();
const unsigned char MaxTraitIndex = std::numeric_limits<TraitIndexType>::max();
const size_t InvalidValue = ~size_t(0);
const size_t MaxIndexSize = std::numeric_limits<size_t>::max();


namespace simuPOP {

// debug codes in a bitset.
std::bitset<DBG_CODE_LENGTH> g_dbgCode;


#ifndef OPTIMIZED
#  include <time.h>
clock_t g_clock;

void initClock()
{
	if (debug(DBG_PROFILE))
		g_clock = clock();
}


void elapsedTime(const string & name)
{
	if (debug(DBG_PROFILE)) {
		cerr << name << ": " << static_cast<double>(clock() - g_clock) / CLOCKS_PER_SEC << "\n";
		g_clock = clock();
	}
}


#endif


const char * g_debugCodes[] = {
	"DBG_ALL",                  // all debug information
	"DBG_GENERAL",              // general information, small warnings etc
	"DBG_UTILITY",              // utility classes and modules
	"DBG_POPULATION",           // debug information related to population
	"DBG_OPERATOR",             // debug information related to operators in general
	"DBG_SIMULATOR",            // debug information related to Simulators
	"DBG_INDIVIDUAL",           // debug information related to individuals
	"DBG_MUTATOR",              // debug information related to mutators
	"DBG_TRANSMITTER",          // debug information related to genotype transmitters
	"DBG_INITIALIZER",          // debug information related to initializer
	"DBG_STATOR",               // debug information related to statistics calculator
	"DBG_TAGGER",               // debug information realted to tagger
	"DBG_SELECTOR",             // debug information related to natural selection operators
	"DBG_MATING",               // debug information related to mating schemes
	"DBG_MIGRATOR",             // debug information related to migration
	"DBG_PROFILE",              // debug information used for performance profiling.
	"DBG_BATCHTESTING",         // debug information used for batch testing of scripts
	"DBG_INTEROPERABILITY",     // debug information used when interoperate with other applications
	"DBG_COMPATIBILITY",        // debug information (obsolete notice etc) about compatibility
	"DBG_DEVEL",                // debug information purely for development.
	"DBG_WARNING",              // warning information given for some particular usage of simuPOP.
	""
};


// set debug area, default to turn all code on
void turnOnDebug(const string & codeString)
{
#ifdef OPTIMIZED
	(void)codeString;  // avoid a warning about unused variable
#else
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
		for (size_t i = 0; i < codes.size(); ++i) {
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
void turnOffDebug(const string & codeString)
{
#ifdef OPTIMIZED
	(void)codeString;  // avoid a warning about unused variable
#else
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
		for (size_t i = 0; i < codes.size(); ++i) {
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


string g_warningMsg;

bool repeatedWarning(const string & message)
{
	if (message != g_warningMsg) {
		g_warningMsg = message;
		return false;
	}
	return true;
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
	if (_Py_RefTotal > g_refTotal && g_refWarningCount-- > 0)
		cerr	<< "Warning: Ref count increased from " << g_refTotal << " to " << _Py_RefTotal
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


/*
   BOOL WINAPI InterruptHandler(DWORD CEvent)
   {
    throw RuntimeError("KeyboardInterruption");
   }
 */

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


// thread number, global variable
UINT g_numThreads;

// random number generator. a global variable.
#ifdef _OPENMP
#  if THREADPRIVATE_SUPPORT == 0
vector<RNG *> g_RNGs;
#  else
RNG * g_RNG;
// each thread has its own g_RNG
#    pragma omp threadprivate(g_RNG)
#  endif
#else
RNG g_RNG;
#endif

void setOptions(const int numThreads, const char * name, unsigned long seed)
{
#ifdef _OPENMP
	// if numThreads is zero, all threads will be used.
	if (numThreads == 0) {
		g_numThreads = omp_get_max_threads();
	} else if (numThreads > 0) {
		omp_set_num_threads(numThreads);
		g_numThreads = numThreads;
	}
#  if THREADPRIVATE_SUPPORT == 0
	g_RNGs.resize(g_numThreads);
	if (seed == 0)
		seed = g_RNGs[0] == NULL ? RNG::generateRandomSeed() : g_RNGs[0]->seed();
	for (unsigned long i = 0; i < g_RNGs.size(); i++) {
		if (g_RNGs[i] == NULL) {
			g_RNGs[i] = new RNG(name, seed + i);
		} else {
			g_RNGs[i]->set(name, seed + i);
		}
	}
#  else
	if (seed == 0)
		seed = g_RNG == NULL ? RNG::generateRandomSeed() : g_RNG->seed();
#    pragma omp parallel
	{
		if (g_RNG == NULL) {
			g_RNG = new RNG(name, seed + omp_get_thread_num());
		} else {
			g_RNG->set(name, seed + omp_get_thread_num());
		}
	}
#  endif
#else
	(void)numThreads;  // avoid an unused parameter warning
	g_RNG.set(name, seed);
#endif
}


UINT numThreads()
{
#ifdef _OPENMP
	return g_numThreads;
#else
	return 1;
#endif
}


ATOMICLONG fetchAndIncrement(ATOMICLONG * val)
{
	if (g_numThreads == 1)
		return (*val)++;
	else
#ifdef _WIN64
		return InterlockedIncrement64(val) - 1;
#elif defined(_WIN32)
		return InterlockedIncrement(val) - 1;
#else
		    // for Intel C++, see page 164 of
		// http://softwarecommunity.intel.com/isn/downloads/softwareproducts/pdfs/347603.pdf
		//
		// for gcc, see
		// http://gcc.gnu.org/onlinedocs/gcc-4.1.0/gcc/Atomic-Builtins.html
		return __sync_fetch_and_add(val, 1);
#endif
}


// return the global RNG
RNG & getRNG()
{
#ifdef _OPENMP
#  if THREADPRIVATE_SUPPORT == 0
	return *g_RNGs[omp_get_thread_num()];
#  else
	return *g_RNG;
#  endif
#else
	return g_RNG;
#endif
}


}

namespace std {
// how to output a dictionary
ostream & operator<<(ostream & out, const strDict & dict)
{
	out << "{";
	if (!dict.empty()) {
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
	if (!dict.empty()) {
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
floatList::floatList(PyObject * obj) : m_elems()
{
	if (obj == NULL)
		return;

	if (PyNumber_Check(obj))
		m_elems.push_back(PyFloat_AsDouble(obj));
	else if (PySequence_Check(obj)) {
		size_t n = PySequence_Size(obj);
		for (size_t j = 0; j < n; ++j) {
			PyObject * val = PySequence_GetItem(obj, j);
			DBG_ASSERT(PyNumber_Check(val), ValueError,
				"A list of numbers is expected");
			m_elems.push_back(PyFloat_AsDouble(val));
			Py_DECREF(val);
		}
	} else {
		DBG_FAILIF(true, ValueError, "Can not create a float list from input.");
	}
}


stringList::stringList(PyObject * obj) : m_elems(), m_allAvail(false), m_trait(MaxTraitIndex)
{
	if (obj == NULL || obj == Py_None)
		m_allAvail = true;
	else if (PyBool_Check(obj))
		// accept True/False
		m_allAvail = obj == Py_True;
	else if (PyString_Check(obj)) {
		string value = PyObj_AsString(obj);
		m_elems.push_back(value);
	}
#if PY_VERSION_HEX >= 0x03000000
	else if (PyBytes_Check(obj)) {
		string value = PyBytes_AsString(obj);
		m_elems.push_back(value);
	}
#endif
	else if (PySequence_Check(obj)) {
		// assign values
		size_t numStr = PySequence_Size(obj);
		for (size_t i = 0; i < numStr; ++i) {
			PyObject * item = PySequence_GetItem(obj, i);
			if (PyString_Check(item)) {
				string value = PyObj_AsString(item);
				m_elems.push_back(value);
			}
#if PY_VERSION_HEX >= 0x03000000
			else if (PyBytes_Check(item)) {
				string value = PyBytes_AsString(item);
				m_elems.push_back(value);
			}
#endif
			else {
				cerr << "A string is expected" << endl;
				throw ValueError("A string is expected.");
			}
			Py_DECREF(item);
		}
	}
}


void stringList::addString(PyObject * str)
{
	PyObject * res = PyObject_Str(str);

	if (res == NULL)
		return;
	string value = PyObj_AsString(res);
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


const vectorstr & stringList::elems(const GenoStruTrait * trait) const
{
	if (trait == NULL || !allAvail())
		return m_elems;

	if (trait) {
		if (trait->genoStruIdx() == m_trait)
			return m_elems;
		m_elems.clear();
		const vectorstr fields = trait->infoFields();
		m_elems.insert(m_elems.end(), fields.begin(), fields.end());
		m_trait = trait->genoStruIdx();
	}
	return m_elems;
}


intMatrix::intMatrix(PyObject * obj) : m_elems()
{
	if (obj == NULL)
		return;
	// Exception will be converted to TypeError...
	if (!PySequence_Check(obj)) {
		cerr << "ERROR: A list or a nested list of integers is expected";
		DBG_ASSERT(false, ValueError,
			"A list or a nested list of intgers is expected");
	}

	size_t numItems = PySequence_Size(obj);
	bool oneDim = true;
	for (size_t i = 0; i < numItems; ++i) {
		PyObject * item = PySequence_GetItem(obj, i);
		oneDim = !PySequence_Check(item);
		Py_DECREF(item);
		if (!oneDim)
			break;
	}
	if (oneDim)
		m_elems.push_back(vectori());
	for (size_t i = 0; i < numItems; ++i) {
		PyObject * item = PySequence_GetItem(obj, i);
		if (PyNumber_Check(item)) {
			if (oneDim)
				m_elems[0].push_back(PyInt_AsLong(item));
			else
				m_elems.push_back(vectori(1, PyInt_AsLong(item)));
		} else if (PySequence_Check(item)) {
			m_elems.push_back(vectori());
			size_t n = PySequence_Size(item);
			for (size_t j = 0; j < n; ++j) {
				PyObject * val = PySequence_GetItem(item, j);
				if (!PyNumber_Check(val)) {
					cerr << "ERROR: A list or nested list of numbers is expected" << endl;
					DBG_ASSERT(false, ValueError,
						"ERROR: A list or nested list of numbers is expected");
				}
				m_elems.back().push_back(PyInt_AsLong(val));
				Py_DECREF(val);
			}
		} else {
			cerr << "ERROR: Can not create a int matrix from input." << endl;
			DBG_ASSERT(false, ValueError,
				"Can not create a int matrix from input.");
		}
		Py_DECREF(item);
	}
}


floatMatrix::floatMatrix(PyObject * obj) : m_elems()
{
	if (obj == NULL)
		return;
	if (!PySequence_Check(obj)) {
		cerr << "A list or a nested list of numbers is expected." << endl;
		DBG_ASSERT(false, ValueError,
			"A list or a nested list of numbers is expected");
	}

	size_t numItems = PySequence_Size(obj);
	for (size_t i = 0; i < numItems; ++i) {
		PyObject * item = PySequence_GetItem(obj, i);
		if (PyNumber_Check(item)) {
			if (m_elems.size() > 1) {
				cerr << "ERROR: A mixture of number and list is not allowed." << endl;
				DBG_ASSERT(false, ValueError,
					"A mixture of number and list is not allowed.")
			}
			if (m_elems.empty())
				m_elems.push_back(vectorf());
			m_elems[0].push_back(PyFloat_AsDouble(item));
		} else if (PySequence_Check(item)) {
			m_elems.push_back(vectorf());
			size_t n = PySequence_Size(item);
			for (size_t j = 0; j < n; ++j) {
				PyObject * val = PySequence_GetItem(item, j);
				if (!PyNumber_Check(val)) {
					cerr << "ERROR: A list or nested list of numbers is expected" << endl;
					DBG_ASSERT(false, ValueError,
						"A list or nested list of numbers is expected");
				}
				m_elems.back().push_back(PyFloat_AsDouble(val));
				Py_DECREF(val);
			}
		} else {
			cerr << "ERROR: Can not create a int matrix from input." << endl;
			DBG_FAILIF(true, ValueError, "Can not create a int matrix from input.");
		}
		Py_DECREF(item);
	}
}


stringMatrix::stringMatrix(PyObject * obj) : m_elems()
{
	if (obj == NULL)
		return;
	DBG_ASSERT(PySequence_Check(obj), ValueError,
		"A list or a nested list of strings is expected");

	size_t numItems = PySequence_Size(obj);
	for (size_t i = 0; i < numItems; ++i) {
		PyObject * item = PySequence_GetItem(obj, i);
		if (PyString_Check(item)) {
			if (m_elems.size() > 1) {
				cerr << "A mixture of string and list is not allowed." << endl;
				throw ValueError("A mixture of string and list is not allowed.");
			}
			if (m_elems.empty())
				m_elems.push_back(vectorstr());
			string value = PyObj_AsString(item);
			m_elems[0].push_back(value);
		}
#if PY_VERSION_HEX >= 0x03000000
		else if (PyBytes_Check(item)) {
			if (m_elems.size() > 1) {
				cerr << "A mixture of string and list is not allowed." << endl;
				throw ValueError("A mixture of string and list is not allowed.");
			}
			if (m_elems.empty())
				m_elems.push_back(vectorstr());
			string value = PyBytes_AsString(item);
			m_elems[0].push_back(value);
		}
#endif
		else if (PySequence_Check(item)) {
			m_elems.push_back(vectorstr());
			size_t numStrs = PySequence_Size(item);
			for (size_t j = 0; j < numStrs; ++j) {
				PyObject * str = PySequence_GetItem(item, j);
				if (PyString_Check(str)) {
					string value = PyObj_AsString(str);
					m_elems.back().push_back(value);
				}
#if PY_VERSION_HEX >= 0x03000000
				else if (PyBytes_Check(str)) {
					string value = PyBytes_AsString(str);
					m_elems.back().push_back(value);
				}
#endif
				else {
					cerr << "A list or nested list of string is expected" << endl;
					throw ValueError("A list or nested list of string is expected.");
				}
				Py_DECREF(str);
			}
		} else {
			DBG_FAILIF(true, ValueError, "Can not create a string matrix from input.");
		}
		Py_DECREF(item);
	}
}



/** A wrapper to a python function
 *  CPPONLY
 */
class CircularReferences
{
public:
	CircularReferences(): m_references()
	{
	}

	void register_ref(PyObject * obj)
	{
		// there is a slight chance that self has two functions...
		// in which case obj should have at least reference two.
		// I am ignoring this case here.
		m_references.insert(obj);
	}

	// this function will be called by the destructor to remove
	// the pointer from this global list. The destructor can
	// be triggered by cleanup(), or at the Python level if the user
	// decided to do something.
	void deregister_ref(PyObject * obj)
	{
		std::set<PyObject *>::iterator it = m_references.find(obj);
		if (it != m_references.end())
			m_references.erase(it);
	}

	// clean up circular refs. Because such objects might exist during
	// the creation of evolutionary scenario, we wait till the end of
	// each evolutionary cycle to clean up the mess
	void cleanup()
	{
		std::set<PyObject *>::iterator it = m_references.begin();
		std::set<PyObject *>::iterator it_end = m_references.end();
		std::vector<PyObject *> to_be_removed;
		for (; it != it_end; ++it)
			// a circular ref has a reference to the function self.func
			// so it has at least one reference. In another word, if we see
			// that self has only one reference, it means the object is not
			// referenced by anyone else and should be removed.
			if ((*it)->ob_refcnt == 1)
				to_be_removed.push_back(*it);

		// now really remove the objects. We do it not within the above loop
		// because the removal of objects will call deregister_ref and
		// change the container.
		std::vector<PyObject *>::iterator iv = to_be_removed.begin();
		std::vector<PyObject *>::iterator iv_end = to_be_removed.end();
		for (; iv != iv_end; ++iv) {
			Py_DECREF(*iv);
			deregister_ref(*iv);
		}
	}

private:
	std::set<PyObject *> m_references;
};

CircularReferences g_circular_refs;

void cleanupCircularRefs()
{
	g_circular_refs.cleanup();
}


pyFunc::pyFunc(PyObject * func) : m_func(func), m_numArgs(0), m_circular_self(NULL)
{
	if (!m_func.isValid())
		return;

	PyObject * obj = m_func.object();

	DBG_ASSERT(PyCallable_Check(obj), ValueError,
		"Passed parameter should be None or a Python function");

	if (PyObject_HasAttrString(obj, "__call__")) {
		if (PyObject_HasAttrString(obj, "__args__")) {
			// in this case, a WithArgs object must have been passed.
			PyObject * args = PyObject_GetAttrString(obj, "__args__");
			m_numArgs = PySequence_Size(args);
			for (size_t i = 0; i < m_numArgs; ++i) {
				PyObject * item = PySequence_GetItem(args, i);
				DBG_ASSERT(PyString_Check(item), ValueError,
					"Attribute args in a simuPOP WithArgs object should be a list of strings");
				m_args.push_back(PyObj_AsString(item));
				Py_DECREF(item);
			}
			Py_DECREF(args);
			// find its name.
			PyObject * func = PyObject_GetAttrString(obj, "__call__");
			DBG_ASSERT(PyCallable_Check(func), ValueError,
				"The func attribute of the passed object should be callable.");
			if (PyObject_HasAttrString(func, "__name__")) {
				PyObject * name = PyObject_GetAttrString(func, "__name__");
				m_name = PyObj_AsString(name);
				Py_DECREF(name);
			} else if (PyObject_HasAttrString(func, "__call__")) {
				PyObject * func1 = PyObject_GetAttrString(func, "__call__");
				if (PyObject_HasAttrString(func1, "__name__")) {
					PyObject * name = PyObject_GetAttrString(func1, "__name__");
					m_name = PyObj_AsString(name);
					Py_DECREF(name);
				}
				Py_DECREF(func1);
			}
			Py_DECREF(func);
			return;
#if PY_VERSION_HEX < 0x03000000
		} else if (!PyObject_HasAttrString(obj, "func_code")) {
#else
		} else if (!PyObject_HasAttrString(obj, "__code__")) {
#endif
			// if there is no arg so it is a member function of a class
			obj = PyObject_GetAttrString(obj, "__call__");
			/*
			   If a functional object is passed (an object with __call__), we keep the reference
			   to obj.__call__ (by GetAttrString) and decrease the reference of the object. This
			   essentially converts func=obj to func=obj.__call__. This allows obj to be
			   destructed after it is used.
			 */
			m_func = pyObject(obj);
			Py_DECREF(obj);
		}
	}
	/* The reference count of the passed function or functor is tricky here. simuPOP
	    allows the pass of functions and functors. There are several cases:

	    1. regular function: func. They are handled directly.
	    2. member function: A().func. This is handled directly.
	    3. functor: A() with __call__. In this case, we keep the reference
	        to obj.__call__ (by GetAttrString) and decrease the reference of the object. This
	        essentially converts func=obj to func=obj.__call__ and allows obj to be
	        destructed after it is used.
	    4. self.func:
	    5. self.__call__:
	    6. self

	    In the last two cases, the references are circular. When the function is passed
	    both the function and the classes will not be released, leading to memory
	    leak.

	    class MyOperator(PyOperator):
	        def __init__(self):
	            // case 4
	            PyOperator.__init__(self, func=self.my_func)
	            // case 5
	            PyOperator.__init__(self, func=self.__call__)
	            // case 6
	            PyOperator.__init__(self, func=self)

	        def __call__(self):
	            return True

	        def my_func(self):
	            return True

	    In both these cases, we can decrease the reference of self to
	    force the operator to be released.
	 */
	// is it bounded?
#if PY_VERSION_HEX < 0x03000000
#  define SELF_ATTR "im_self"
#  define CODE_ATTR "func_code"
#else
#  define SELF_ATTR "__self__"
#  define CODE_ATTR "__code__"
#endif
	int bounded = PyObject_HasAttrString(obj, SELF_ATTR);
	if (bounded) {
		PyObject * self = PyObject_GetAttrString(obj, SELF_ATTR);
		//
		// This check is not very rigrious because other classes and
		// also have these attributes, but it is difficult to really
		// check the super class (BaseOperator) because of the SWIG
		// interface
		if (PyObject_HasAttrString(self, "apply") &&
		    PyObject_HasAttrString(self, "describe")) {
			m_circular_self = self;
			g_circular_refs.register_ref(self);
		}
		Py_DECREF(self);
	}
	if (!PyObject_HasAttrString(obj, "__name__")) {
		cerr << "Cannot find name of the passed function. " << endl;
		throw ValueError("Cannot find name of the passed function.");
	}
	// find its name.
	PyObject * name = PyObject_GetAttrString(obj, "__name__");
	m_name = PyObj_AsString(name);
	Py_DECREF(name);

	// free python functions have a 'func_code' attribute
	// built-in functions might not have (e.g. random.random)
	if (!PyObject_HasAttrString(obj, CODE_ATTR))
		return;
	PyObject * code = PyObject_GetAttrString(obj, CODE_ATTR);

	if (!code) {
		cerr << "Invalid attribute func_code or __code__ for a function object" << endl;
		throw SystemError("Invalid attribute func_code or __code for a function object");
	}
	// probe number of parameters
	PyObject * co_argcount = PyObject_GetAttrString(code, "co_argcount");
	DBG_ASSERT(co_argcount, SystemError, "Invalid attribute co_argcount for a function object");
	// substract 1 if the method is bounded to remove the count for self.
	m_numArgs = PyInt_AsLong(co_argcount) - bounded;
	Py_DECREF(co_argcount);
	// probe parameter names
	PyObject * co_varnames = PyObject_GetAttrString(code, "co_varnames");
	DBG_ASSERT(co_varnames, SystemError, "Invalid attribute co_varnames for a function object");
	for (size_t i = 0; i < m_numArgs; ++i) {
		PyObject * item = PyTuple_GetItem(co_varnames, i + bounded);
		m_args.push_back(PyObj_AsString(item));
	}
	Py_DECREF(co_varnames);
	// accepting arbitrary number of parameters?
	/*
	   PyObject * co_flag = PyObject_GetAttrString(code, "co_flags");
	   DBG_ASSERT(co_flag, SystemError, "Invalid attribute co_flags for a function object");
	   int flags = static_cast<unsigned char>(PyInt_AsLong(co_flag));
	   Py_DECREF(co_flag);
	 */
	Py_DECREF(code);
}

void pyGenerator::set(PyObject * gen)
{
	Py_XDECREF(m_iterator);
	Py_XDECREF(m_generator);

	if (!gen) {
		m_iterator = NULL;
		m_generator = NULL;
		return;
	}

	m_generator = gen;

	// test if m_generator is a generator
	DBG_ASSERT(PyGen_Check(m_generator), ValueError,
		"Passed function is not a python generator");

	m_iterator = PyObject_GetIter(m_generator);

	// test if m_iterator is iteratable.
	DBG_ASSERT(m_iterator, RuntimeError, "Can not create an iterate from a generator");
}


PyObject * pyGenerator::next()
{
	PyObject * obj = PyIter_Next(m_iterator);

#ifndef OPTIMIZED
	if (PyErr_Occurred()) {
		PyErr_Print();
		PyErr_Clear();
	}
#endif
	DBG_ASSERT(obj, RuntimeError, "Iterator returns NULL item.");
	return obj;
}


uintList::uintList(PyObject * obj) : m_elems(), m_status(REGULAR)
{
	if (obj == NULL)
		// accept NULL
		m_status = UNSPECIFIED;
	else if (PyBool_Check(obj))
		// accept True/False
		m_status = obj == Py_True ? ALL_AVAIL : UNSPECIFIED;
	else if (PyNumber_Check(obj)) {
		// accept a number
		m_elems.push_back(static_cast<UINT>(PyInt_AsLong(obj)));
	} else if (PySequence_Check(obj)) {
		m_elems.resize(PySequence_Size(obj));
		// assign values
		for (size_t i = 0, iEnd = m_elems.size(); i < iEnd; ++i) {
			PyObject * item = PySequence_GetItem(obj, i);
			DBG_ASSERT(PyNumber_Check(item), ValueError, "Invalid input for a list of integers.");
			m_elems[i] = static_cast<UINT>(PyInt_AsLong(item));
			Py_DECREF(item);
		}
	} else {
		DBG_FAILIF(true, ValueError, "Invalid input for a list of integers.");
	}
}


lociList::lociList(PyObject * obj) : m_elems(), m_names(), m_func(NULL), m_func_gen(NOT_FOUND), m_status(REGULAR), m_trait(MaxTraitIndex), m_lociMap()
{
	if (obj == NULL)
		// accept NULL
		m_status = UNSPECIFIED;
	else if (PyBool_Check(obj))
		// accept True/False
		m_status = obj == Py_True ? ALL_AVAIL : UNSPECIFIED;
	else if (PyString_Check(obj)) {
		m_status = FROM_NAME;
		m_elems.resize(1);
		m_names.push_back(PyObj_AsString(obj));
#if PY_VERSION_HEX < 0x03000000
	} else if (PyUnicode_Check(obj)) {
		m_status = FROM_NAME;
		m_elems.resize(1);
		m_names.push_back(PyObj_AsString(obj));
#endif
		/* for some reason I do not understand, if we pass a
		   Python callable object (obj with _-call__, PyNumber_Check
		   will be ok so PyCallable_Check must be checked before
		   that case. */
	} else if (PyCallable_Check(obj)) {
		// if a function is provided
		m_status = FROM_FUNC;
		m_func = pyFunc(obj);
	} else if (PyNumber_Check(obj)) {
		m_status = REGULAR;
		// accept a number
		m_elems.push_back(static_cast<UINT>(PyInt_AsLong(obj)));
	} else if (PySequence_Check(obj)) {
		m_elems.resize(PySequence_Size(obj));
		// assign values
		for (size_t i = 0, iEnd = m_elems.size(); i < iEnd; ++i) {
			PyObject * item = PySequence_GetItem(obj, i);
			if (PyNumber_Check(item)) {
				if (m_status == FROM_NAME && m_elems.size() == 2) {
					// a special case for loci=(chr, pos)
					m_elems.resize(1);
					m_status = FROM_POSITION;
					m_positions.push_back(genomic_pos(m_names[0], PyFloat_AsDouble(item)));
					m_names.clear();
				} else {
					DBG_FAILIF(i != 0 && m_status != REGULAR, ValueError, "Cannot mix index and loci names.");
					m_status = REGULAR;
					m_elems[i] = static_cast<UINT>(PyInt_AsLong(item));
				}
			} else if (PyString_Check(item)) {
				DBG_FAILIF(i != 0 && m_status != FROM_NAME, ValueError, "Cannot mix index and loci names.");
				m_status = FROM_NAME;
				m_names.push_back(PyObj_AsString(item));
#if PY_VERSION_HEX < 0x03000000
			} else if (PyUnicode_Check(item)) {
				DBG_FAILIF(i != 0 && m_status != FROM_NAME, ValueError, "Cannot mix index and loci names.");
				m_status = FROM_NAME;
				m_names.push_back(PyObj_AsString(item));
#endif
			} else if (PySequence_Check(item)) {
				// a sequence (chr, pos) is acceptable
				DBG_FAILIF(PySequence_Size(item) != 2, ValueError, "Loci, if given as a nested list, should contain a list of (chr, pos) pairs.");
				PyObject * chr = PySequence_GetItem(item, 0);
				PyObject * pos = PySequence_GetItem(item, 1);
#if PY_VERSION_HEX < 0x03000000
				DBG_ASSERT((PyString_Check(chr) || PyUnicode_Check(chr)) && PyNumber_Check(pos), ValueError, "Loci, if given as a nested list, should contain a list of (chr, pos) pair");
#else
				DBG_ASSERT(PyString_Check(chr) && PyNumber_Check(pos), ValueError, "Loci, if given as a nested list, should contain a list of (chr, pos) pair");
#endif
				m_status = FROM_POSITION;
				m_positions.push_back(genomic_pos(PyObj_AsString(chr), PyFloat_AsDouble(pos)));
				Py_DECREF(chr);
				Py_DECREF(pos);
			} else {
				DBG_ASSERT(false, ValueError, "Invalid input for a list of loci (index or name should be used).");
			}
			Py_DECREF(item);
		}
	} else {
		DBG_FAILIF(true, ValueError, "Invalid input for a list of integers.");
	}
}


const vectoru & lociList::elems(const GenoStruTrait * trait) const
{
	if (trait) {
		if (trait->genoStruIdx() == m_trait)
			return m_elems;
		if (m_status == FROM_NAME)
			m_elems = trait->lociByNames(m_names);
		else if (m_status == FROM_POSITION)
			m_elems = trait->lociByPos(m_positions);
		else if (m_status == FROM_FUNC)
			throw ValueError("Calling a function for lociList from this operator is not allowed.");
		else if (m_status == ALL_AVAIL) {
			m_elems.resize(trait->totNumLoci());
			for (size_t i = 0; i < m_elems.size(); ++i)
				m_elems[i] = i;
		}
		m_trait = trait->genoStruIdx();
	}
	return m_elems;
}


const vectoru & lociList::elems(const Population * trait) const
{
	if (trait) {
		if (trait->genoStruIdx() == m_trait && m_status != FROM_FUNC)
			return m_elems;
		if (m_status == FROM_NAME)
			m_elems = trait->lociByNames(m_names);
		else if (m_status == FROM_POSITION)
			m_elems = trait->lociByPos(m_positions);
		else if (m_status == FROM_FUNC) {
			if (m_func_gen == trait->gen())
				return m_elems;

			PyObject * args = PyTuple_New(m_func.numArgs());
			DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

			for (size_t i = 0; i < m_func.numArgs(); ++i) {
				const string & arg = m_func.arg(i);
				if (arg == "pop")
					PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(const_cast<Population *>((trait)))));
				else {
					DBG_FAILIF(true, ValueError,
						"Only parameter pop are acceptable in function " + m_func.name());
				}
			}

			m_elems = m_func(PyObj_As_SizeTArray, args);
			Py_XDECREF(args);
			m_func_gen = trait->gen();
		} else if (m_status == ALL_AVAIL) {
			m_elems.resize(trait->totNumLoci());
			for (size_t i = 0; i < m_elems.size(); ++i)
				m_elems[i] = i;
		}
		m_trait = trait->genoStruIdx();
	}
	return m_elems;
}


size_t lociList::indexOf(size_t loc) const
{
	if (m_status == ALL_AVAIL)
		return loc;
	if (m_lociMap.empty())
		for (size_t i = 0; i < m_elems.size(); ++i)
			m_lociMap[m_elems[i]] = i;
	std::map<size_t, size_t>::iterator it = m_lociMap.find(loc);
	if (it == m_lociMap.end())
		return NOT_FOUND;
	else
		return it->second;
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
		m_elems.push_back(PyInt_AsLong(obj));
	} else if (PySequence_Check(obj)) {
		m_elems.resize(PySequence_Size(obj));
		// assign values
		for (size_t i = 0, iEnd = m_elems.size(); i < iEnd; ++i) {
			PyObject * item = PySequence_GetItem(obj, i);
			DBG_ASSERT(PyNumber_Check(item), ValueError, "Invalid input for a list of rep.");
			m_elems[i] = PyInt_AsLong(item);
			Py_DECREF(item);
		}
	} else {
		DBG_FAILIF(true, ValueError, "Invalid input for a list of rep.");
	}
}


floatListFunc::floatListFunc(PyObject * obj) :
	floatList(PyCallable_Check(obj) ? NULL : obj),
	m_func(PyCallable_Check(obj) ? obj : NULL)
{
}


//
// shared variables
//

// utility functions: python => C++ conversion
void PyObj_As_Bool(PyObject * obj, bool & val)
{
	val = PyObject_IsTrue(obj) == 1;
}


void PyObj_As_Int(PyObject * obj, long & val)
{
	if (obj == NULL) {
		val = 0;
		return;
	}
	// try to convert
	PyObject * res = PyNumber_Long(obj);
	if (res == NULL)
		throw ValueError("Can not convert object to an integer");

	val = PyInt_AsLong(res);
	Py_DECREF(res);
}


void PyObj_As_SizeT(PyObject * obj, size_t & val)
{
	if (obj == NULL) {
		val = 0;
		return;
	}
	// try to convert
	PyObject * res = PyNumber_Long(obj);
	if (res == NULL)
		throw ValueError("Can not convert object to an integer");

	val = PyInt_AsLong(res);
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

	val = PyObj_AsString(res);
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


void PyObj_As_SizeTArray(PyObject * obj, vectoru & val)
{
	if (obj == NULL) {
		val = vectoru();
		return;
	}
	if (PySequence_Check(obj)) {
		val.resize(PySequence_Size(obj));

		// assign values
		for (size_t i = 0, iEnd = val.size(); i < iEnd; ++i) {
			PyObject * item = PySequence_GetItem(obj, i);
			PyObj_As_SizeT(item, val[i]);
			Py_DECREF(item);
		}
	} else {
		val.resize(1);
		PyObj_As_SizeT(obj, val[0]);
	}
}


PyObject * Allele_Vec_As_NumArray(GenoIterator begin, GenoIterator end)
{
	PyObject * res = newcarrayobject(begin, end);

	DBG_FAILIF(res == NULL, ValueError, "Can not convert buf to Allele num array");
	return res;
}


PyObject * Lineage_Vec_As_NumArray(LineageIterator begin, LineageIterator end)
{
	PyObject * res = newcarrayobject_lineage(begin, end);

	DBG_FAILIF(res == NULL, ValueError, "Can not convert buf to Lineage num array");
	return res;
}


string PyObj_AsString(PyObject * str)
{
#if PY_VERSION_HEX >= 0x03000000
	char * cstr;
	char * newstr;
	Py_ssize_t len;
	str = PyUnicode_AsUTF8String(str);
	PyBytes_AsStringAndSize(str, &cstr, &len);
	newstr = (char *)malloc(len + 1);
	memcpy(newstr, cstr, len + 1);
	Py_XDECREF(str);
	string res(newstr);
	free(newstr);
	return res;
#else
	return string(PyString_AsString(str));
#endif
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
	// find the first piece
	size_t i, s;

	for (i = 0; i < name.size() && name[i] != '[' && name[i] != '{'; ++i) ;

	if (i == 0)
		throw ValueError("Empty name? " + name);

	// we need to keep current parent, key/index, type (0 for array, 1 for dict)
	// keytype (if curType = 1) (0 for string, 1 for num)
	size_t curType = 1;
	PyObject * curParent = m_dict;                            // should be always valid
	PyObject * curKey = PyString_FromString(const_cast<char *>(name.substr(0, i).c_str()));
	size_t curIdx = 0;
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
			if (!isdigit(name[i]))
				throw ValueError("Expecting numbers: " + name);

		// get index
		size_t idx = atoi(name.substr(s, i - s).c_str());
		// not exist
		//
		// create subPop, ...
		if (curChild == NULL || !PyList_Check(curChild)) {
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
		//                 float   3: Ne_LD[0.01]
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

		for ( ; name[i] != '}' && i < name.size(); ++i) {
			if (name[i] == '.' && keyType == 0)
				keyType = 3;
		}

		DBG_ASSERT(name[i] == '}', ValueError, "Unmatched dictionary delimiter");

		PyObject * childKey;

		if (keyType == 0)
			childKey = PyInt_FromString(const_cast<char *>(name.substr(s, i - s).c_str()), NULL, 0);
		else if (keyType == 3) {
			PyObject * key_str = PyString_FromString(const_cast<char *>(name.substr(s, i - s).c_str()));
			childKey = PyFloat_FromString(key_str, NULL);
		} else if (keyType == 1) {
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
		if (curChild == NULL || !PyDict_Check(curChild)) {
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


PyObject * SharedVariables::getVar(const string & name, bool nameError) const
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
			if (!isdigit(name[i]))
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
		//                 float   3:
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

		for ( ; name[i] != '}' && i < name.size(); ++i) {
			if (name[i] == '.' && keyType == '0')
				keyType = 3;
		}

		assert(name[i] == '}');

		PyObject * childKey;
		if (keyType == 0)
			childKey = PyInt_FromString(const_cast<char *>(name.substr(s, i - s).c_str()), NULL, 0);
		else if (keyType == 3) {
			PyObject * key_str = PyString_FromString(const_cast<char *>(name.substr(s, i - s).c_str()));
			childKey = PyFloat_FromString(key_str, NULL);
		}else if (keyType == 1) {
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
			if (!isdigit(name[i]))
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


PyObject * SharedVariables::setVar(const string & name, const bool val)
{
	PyObject * obj = val ? Py_True : Py_False;

	Py_INCREF(obj);
	return setVar(name, obj);
}


PyObject * SharedVariables::setVar(const string & name, const long val)
{
	return setVar(name, PyInt_FromLong(val));
}


PyObject * SharedVariables::setVar(const string & name, const size_t val)
{
	return setVar(name, PyInt_FromSize_t(val));
}


PyObject * SharedVariables::setVar(const string & name, const double val)
{
	return setVar(name, PyFloat_FromDouble(val));
}


PyObject * SharedVariables::setVar(const string & name, const string & val)
{
	return setVar(name, Py_BuildValue("s", val.c_str()));
}


PyObject * SharedVariables::setVar(const string & name, const vectori & val)
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


PyObject * SharedVariables::setVar(const string & name, const vectoru & val)
{
	PyObject * obj = PyList_New(0);
	PyObject * item;

	for (vectoru::const_iterator it = val.begin();
	     it < val.end(); ++it) {
		item = PyInt_FromSize_t(*it);
		PyList_Append(obj, item);
		Py_XDECREF(item);
	}
	return setVar(name, obj);
}


//CPPONLY
PyObject * SharedVariables::setVar(const string & name, const vectorf & val)
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


PyObject * SharedVariables::setVar(const string & name, const strDict & val)
{
	PyObject * obj = PyDict_New();
	PyObject * v;

	for (strDict::const_iterator i = val.begin(); i != val.end(); ++i) {
		PyDict_SetItemString(obj, const_cast<char *>(i->first.c_str()),
			// SetItem will increase ref count for v
			// for Py_None, this is good enough, for PyFloat, we need to
			// decrease its ref count
			v = i->second == MISSING_VALUE ? Py_None : PyFloat_FromDouble(i->second));
		if (v != Py_None)
			Py_XDECREF(v);
	}
	return setVar(name, obj);
}


PyObject * SharedVariables::setVar(const string & name, const intDict & val)
{
	PyObject * obj = PyDefDict_New();
	PyObject * u, * v;

	for (intDict::const_iterator i = val.begin(); i != val.end(); ++i) {
		PyDict_SetItem(obj,
			u = PyInt_FromLong(i->first),
			v = i->second == MISSING_VALUE ? Py_None : PyFloat_FromDouble(i->second));
		Py_XDECREF(u);
		if (v != Py_None)
			Py_XDECREF(v);
	}
	return setVar(name, obj);
}


void SharedVariables::getVarAsIntDict(const string & name, uintDict & res, bool nameError) const
{
	res.clear();
	PyObject * obj = getVar(name, nameError);

	PyObject * key, * value;
	Py_ssize_t pos = 0;

	while (PyDict_Next(obj, &pos, &key, &value)) {
		int k = PyInt_AS_LONG(key);
		double v = PyFloat_AsDouble(value);
		res[k] = v;
	}
}


void SharedVariables::getVectorVarAsIntDict(const string & name, uintDict & res, bool nameError) const
{
	res.clear();
	PyObject * obj = getVar(name, nameError);

	Py_ssize_t sz = PySequence_Size(obj);
	for (Py_ssize_t i = 0; i < sz; ++i) {
		PyObject * item = PySequence_GetItem(obj, i);
		size_t k = PyInt_AsLong(item);
		res[k] = 0;
		Py_XDECREF(item);
	}
}


PyObject * SharedVariables::setVar(const string & name, const uintDict & val)
{
	PyObject * obj = PyDefDict_New();
	PyObject * u, * v;

	for (uintDict::const_iterator i = val.begin(); i != val.end(); ++i) {
		PyDict_SetItem(obj,
			u = PyInt_FromSize_t(i->first),
			v = PyFloat_FromDouble(i->second));
		Py_XDECREF(u);
		Py_XDECREF(v);
	}
	return setVar(name, obj);
}


PyObject * SharedVariables::setVar(const string & name, const tupleDict & val)
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


void save_none(ostringstream & str)
{
	str << "n";
}


PyObject * load_none(const string & /* str */, size_t & offset)
{
	offset++;
	Py_INCREF(Py_None);
	return Py_None;
}


void save_int(ostringstream & str, PyObject * args)
{
	long l = PyInt_AsLong(args);

	// type + string + ' '
	str << 'i' << l << ' ';
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


void save_long(ostringstream & str, PyObject * args)
{
	long l = PyInt_AsLong(args);

	// type +  string + ' '
	str << 'l' << l << ' ';
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


void save_float(ostringstream & str, PyObject * args)
{
	double d = PyFloat_AsDouble(args);

	// type + string
	str << 'f' << boost::lexical_cast<std::string>(d) << ' ';
}


PyObject * load_float(const string & str, size_t & offset)
{
	int len = 0;

	while (str[offset + len + 1] != ' ') len++;
	double d = atof(const_cast<char *>(str.substr(offset + 1, len).c_str()));
	offset += len + 2;
	return PyFloat_FromDouble(d);
}


void save_string(ostringstream & str, PyObject * args)
{
	str << 's' << PyObj_AsString(args) << '\0';
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


void saveObj(ostringstream & str, PyObject * args);

PyObject * loadObj(const string & vars, size_t & offset);

void save_dict(ostringstream & str, PyObject * args)
{
	PyObject * key = 0, * value = 0;

	str << 'd';                                                                       // dictionary
	Py_ssize_t i = 0;
	while (PyDict_Next(args, &i, &key, &value)) {
		if (PyString_Check(key)) {
			string k = PyObj_AsString(key);
			// do not save things like __builtins__
			if (!k.empty() && k[0] == '_')
				continue;
		}
		saveObj(str, key);
		saveObj(str, value);
	}
	str << 'e';                                                                       // ending of a dictionary
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


void save_defdict(ostringstream & str, PyObject * args)
{
	PyObject * key = 0, * value = 0;

	str << 'D';                                                                       // dictionary
	Py_ssize_t i = 0;
	while (PyDict_Next(args, &i, &key, &value)) {
		saveObj(str, key);
		saveObj(str, value);
	}
	str << 'e';                                                                       // ending of a dictionary
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


void save_list(ostringstream & str, PyObject * args)
{
	str << 'L';                                                                       // dictionary
	size_t len = PyList_Size(args);
	for (size_t i = 0; i < len; i++) {
		PyObject * elem = PyList_GET_ITEM((PyListObject *)args, i);
		saveObj(str, elem);
	}

	str << 'e';                                                                       // ending of a dictionary
}


PyObject * load_list(const string & vars, size_t & offset)
{
	// skip 'L'
	offset++;
	PyObject * d = PyList_New(0);
	while (vars[offset] != 'e') {
		PyObject * elem = loadObj(vars, offset);
		PyList_Append(d, elem);
		Py_DECREF(elem);
	}
	offset++;                                                                         // skip 'e'
	return d;
}


void save_tuple(ostringstream & str, PyObject * args)
{
	str << 't';                                                                       // dictionary
	size_t len = PyTuple_Size(args);
	// save length
	str << len << ' ';
	// save items
	for (size_t i = 0; i < len; i++) {
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
		// this is not needed because PyTuple_SetItem steals a
		// referene from elem.
		//Py_DECREF(elem);
	}
	return d;
}


void saveObj(ostringstream & str, PyObject * args)
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
		DBG_DO(DBG_UTILITY, cerr << boost::format("Warning: object of type '%1%' cannot be saved. Use none.") % type->tp_name);
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
		throw ValueError((boost::format("Unknown type code at offset %1%") % offset).str());
	}
	}
}


string SharedVariables::to_pickle() const
{
#if PY_VERSION_HEX >= 0x03000000
	PyObject * pickle = PyImport_ImportModule("pickle");
#else
	PyObject * pickle = PyImport_ImportModule("cPickle");
#endif
	if (! pickle)
		throw RuntimeError("Failed to import module pickle to serialize population variables.");

	// here we use version 2 because this is the latest version that supported by
	// both python 2 and python 3, also because it is the one that handles simuPOP's
	// defdict type using its __reduce__ interface.
	char dumps[] = "dumps";
	char Oi[] = "(Oi)";
	PyObject * pres = PyObject_CallMethod(pickle, dumps, Oi, m_dict, 2);
	if (pres == NULL) {
		PyErr_Clear();
		/* If the dictionary is not pickleable, we have to go a longer way and test if
		   each item of the dictionary is pickleable. */
		PyObject * new_dict = PyDict_Copy(m_dict);
		// Some items in __builtins__ are not pickable so we will have to remove them
		// before picking.
		PyObject *key, *value;
		Py_ssize_t pos = 0;
		while (PyDict_Next(m_dict, &pos, &key, &value)) {
			// try to pickle each item
			// if the key is not pickleable, remove the key from the copied dictionary
			// to avoid a warning
			if (PyObject_CallMethod(pickle, dumps, Oi, key, 2) == NULL) {
				PyErr_Clear();
				PyDict_DelItem(new_dict, key);
			}
			if (PyObject_CallMethod(pickle, dumps, Oi, value, 2) == NULL) {
				PyErr_Clear();
				PyDict_DelItem(new_dict, key);
			}
		}
		// now try to pickale the whole new object, in any case, the original
		// dictionary is not touched.
		pres = PyObject_CallMethod(pickle, dumps, Oi, new_dict, 2);
		Py_DECREF(new_dict);
		if (pres == NULL) {
			PyErr_Print();
			PyErr_Clear();
			throw RuntimeError("Failed to call pickle.dumps to save population variables.");
		}
	}
	Py_ssize_t sz = 0;
	char * buf = NULL;
#if PY_VERSION_HEX >= 0x03000000
	PyBytes_AsStringAndSize(pres, &buf, &sz);
#else
	PyString_AsStringAndSize(pres, &buf, &sz);
#endif
	// need to copy the data out before releasing pres
	string res(buf, sz);
	Py_DECREF(pres);
	Py_DECREF(pickle);

	return res;
}


string SharedVariables::asString() const
{
	// go through each variable and save
	ostringstream str;

	saveObj(str, m_dict);
	str << 'e';                                                                       // ending character
	return str.str();
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

void SharedVariables::from_pickle(const string & vars)
{
#if PY_VERSION_HEX >= 0x03000000
	PyObject * pickle = PyImport_ImportModule("pickle");
#else
	PyObject * pickle = PyImport_ImportModule("cPickle");
#endif
	if (! pickle)
		throw RuntimeError("Failed to import module pickle to serialize population variables.");

	// use protocol 1 for compatibility
	string res;
#if PY_VERSION_HEX >= 0x03000000
	PyObject * args = PyBytes_FromStringAndSize(vars.c_str(), vars.size());
#else
	PyObject * args = PyString_FromStringAndSize(vars.c_str(), vars.size());
#endif
	// remove m_dict
	if (m_ownVars) {
		PyDict_Clear(m_dict);
		Py_XDECREF(m_dict);
	}
	m_ownVars = true;
	char loads[] = "loads";
	char O[] = "(O)";
	m_dict = PyObject_CallMethod(pickle, loads, O, args);
	if (m_dict == NULL) {
		PyErr_Print();
		PyErr_Clear();
		throw RuntimeError("Failed to call pickle.loads to load population variables.");
	}
	Py_DECREF(args);

	Py_DECREF(pickle);
}

// simuVars will hold replicate specific Shared Variables.

// global dictionary
// will be set by // initialize() function
// DO NOT OWN the dictionaries
SharedVariables g_main_vars, g_module_vars;

swig_type_info * g_swigPopType, * g_swigIndividual, * g_swigOpType;

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
	return SWIG_NewPointerObj(p, g_swigIndividual, 0);
}


void * pyIndPointer(PyObject * obj)
{
	void * ptr = 0;

	SWIG_Python_ConvertPtr(obj, &ptr, g_swigIndividual, SWIG_POINTER_DISOWN);
	return ptr;
}


void * pyOpPointer(PyObject * obj)
{
	void * ptr = 0;

	SWIG_Python_ConvertPtr(obj, &ptr, g_swigOpType, 0);
	return ptr;
}

void * pyPopPointer(PyObject * obj)
{
	void * ptr = 0;

	SWIG_Python_ConvertPtr(obj, &ptr, g_swigPopType, 0);
	return ptr;
}


string shorten(const string & val, size_t length)
{
	return (val.size() > length) ? val.substr(0, length) + "..." : val;
}


// Expression evaluation
// because of ref count, need to define copier

Expression::Expression(const Expression & rhs)
	: m_exprString(rhs.m_exprString), m_stmtsString(rhs.m_stmtsString),
	m_expr(rhs.m_expr), m_stmts(rhs.m_stmts), m_locals(rhs.m_locals)
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

	m_exprString = expr;

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

	m_stmtsString = stmts;

	if (stmts.empty())
		return;

	// clear leading blank etc for statements ? Not yet.
	m_stmts = Py_CompileString(const_cast<char *>(stmts.c_str()), "<embed>",
		Py_file_input);

	if (m_stmts == NULL)
		throw ValueError("statement '" + stmts + "' is not valid.");
}


// python expression
PyObject * Expression::evaluate() const
{
	if (m_expr == NULL && m_stmts == NULL)
		return NULL;

	DBG_ASSERT(mainVars().dict() != NULL && m_locals != NULL,
		ValueError, "Can not evalulate. Dictionary is empty!");

	PyObject * res = NULL;
	if (m_stmts != NULL) {
#if PY_VERSION_HEX >= 0x03020000
		res = PyEval_EvalCode((PyObject *)m_stmts, m_locals, m_locals);
#else
		res = PyEval_EvalCode((PyCodeObject *)m_stmts, m_locals, m_locals);
#endif
		if (res == NULL) {
#ifndef OPTIMIZED
			if (debug(DBG_GENERAL)) {
				PyErr_Print();
				PyErr_Clear();
			}
#endif
			throw RuntimeError("Evalulation of statements '" + m_stmtsString + "' failed");
		} else {
			Py_DECREF(res);
			res = NULL;
		}
	}

	if (m_expr != NULL) {
#if PY_VERSION_HEX >= 0x03020000
		res = PyEval_EvalCode((PyObject *)m_expr, m_locals, m_locals);
#else
		res = PyEval_EvalCode((PyCodeObject *)m_expr, m_locals, m_locals);
#endif
		if (res == NULL) {
#ifndef OPTIMIZED
			if (debug(DBG_GENERAL)) {
				PyErr_Print();
				PyErr_Clear();
			}
#endif
			throw RuntimeError("Evalulation of expression '" + m_exprString + "' failed");
		}
	}

	return res;
}


bool Expression::valueAsBool() const
{
	PyObject * res = evaluate();

	if (res == NULL)
		return false;
	bool val;
	PyObj_As_Bool(res, val);
	Py_XDECREF(res);
	return val;
}


long Expression::valueAsInt() const
{
	PyObject * res = evaluate();

	if (res == NULL)
		return 0;
	long val;
	PyObj_As_Int(res, val);
	Py_XDECREF(res);
	return val;
}


double Expression::valueAsDouble() const
{
	PyObject * res = evaluate();

	if (res == NULL)
		return 0;
	double val;
	PyObj_As_Double(res, val);
	Py_XDECREF(res);
	return val;
}


string Expression::valueAsString() const
{
	PyObject * res = evaluate();

	if (res == NULL)
		return "";
	string val;
	PyObj_As_String(res, val);
	Py_XDECREF(res);
	return val;
}


vectorf Expression::valueAsArray() const
{
	PyObject * res = evaluate();

	if (res == NULL)
		return vectorf();
	vectorf val;
	PyObj_As_Array(res, val);
	Py_XDECREF(res);
	return val;
}


simpleStmt::simpleStmt(const string & stmt, const string & indVar) : m_var(""),
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
	const regex multiplied3("\\s*([\\d\\w]+)\\s*\\=\\s*([-]*[1-9\\.]+)\\s*\\*\\s*\\1\\s*");
	const regex setSex("\\s*([\\d\\w]+)\\s*\\=\\s*([\\w]+).sex\\s*[(]\\s*[)]\\s*");
	const regex setAffected("\\s*([\\d\\w]+)\\s*\\=\\s*([\\w]+).affected\\s*[(]\\s*[)]\\s*");
	const regex setUnaffected("\\s*([\\d\\w]+)\\s*\\=\\s*not\\s+([\\w]+).affected\\s*[(]\\s*[)]\\s*");

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
	else if (regex_match(stmt.c_str(), matches, setSex)) {
		string name = string(matches[2].first, matches[2].second);
		if (name == indVar)
			m_operation = SetSex;
	} else if (regex_match(stmt.c_str(), matches, setAffected)) {
		string name = string(matches[2].first, matches[2].second);
		if (name == indVar)
			m_operation = SetAffected;
	} else if (regex_match(stmt.c_str(), matches, setUnaffected)) {
		string name = string(matches[2].first, matches[2].second);
		if (name == indVar)
			m_operation = SetUnaffected;
	} else
		return;

	m_var = string(matches[1].first, matches[1].second);
	if (m_operation != SetSex && m_operation != SetAffected && m_operation != SetUnaffected) {
		try {
			m_value = boost::lexical_cast<double>(string(matches[2].first, matches[2].second));
		} catch (...) {
			m_operation = NoOperation;
			return;
		}
	}
	DBG_DO(DBG_DEVEL, cerr	<< "Match statement with name " << m_var
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

	if (it == m_ostreams.end()) {                            // not found.

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
		else if (realAppend && !it->second.append())
			it->second.makeAppend(true);
		else if (!realAppend && it->second.append())
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
	for (ostreamMapIterator it = m_ostreams.begin(), itEnd = m_ostreams.end(); it != itEnd; ++it)
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
StreamProvider::StreamProvider(const string & output, const pyFunc & func, const string & mode)
	: m_filename(output), m_filenameExpr(), m_func(func), m_flags(0), m_mode(mode), m_filePtr(NULL)
{
	if (!m_filename.empty() && m_filename[0] == '!') {
		m_filenameExpr.setExpr(m_filename.substr(1));
		m_filename.clear();
	} else if (m_func.isValid()) {
		SETFLAG(m_flags, m_flagUseFunc);
		SETFLAG(m_flags, m_flagCloseAfterUse);
	} else
		analyzeOutputString(m_filename);
	if (m_mode != "" && m_mode != "b")
		throw ValueError("Only empty or b mode are supported.");
}


ostream & StreamProvider::getOstream(PyObject * dict, bool readable)
{
	DBG_FAILIF(readable && (ISSETFLAG(m_flags, m_flagNoOutput) || ISSETFLAG(m_flags, m_flagUseDefault)),
		SystemError, "A readable file is requested but this Opertor uses cerr or cnull.");

	if (ISSETFLAG(m_flags, m_flagNoOutput))
		return cnull();

	// if using cerr, return it.
	if (ISSETFLAG(m_flags, m_flagUseDefault))
		return std::cout;

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
		if (ISSETFLAG(m_flags, m_flagNoOutput))
			return cnull();

		// if using cerr, return it.
		if (ISSETFLAG(m_flags, m_flagUseDefault))
			return cerr;
	}
	filename = m_filename;

	if (ISSETFLAG(m_flags, m_flagAppend)) {

		DBG_DO(DBG_UTILITY, cerr	<< "Get a persistent file: "
			                        << filename << endl);

		return *ostreamManager().getOstream(filename, readable,
			ISSETFLAG(m_flags, m_flagRealAppend), ISSETFLAG(m_flags, m_flagUseString));
	} else {                                                                          // not in append mode, but check if this file is alreay there

		DBG_DO(DBG_UTILITY, cerr	<< "File is not persistent : "
			                        << filename << endl);

		if (!ostreamManager().hasOstream(filename)) {
			if (readable)
				SETFLAG(m_flags, m_flagReadable);
			else
				RESETFLAG(m_flags, m_flagReadable);

			if (readable)
				m_filePtr = new fstream(filename.c_str());
			else
				m_filePtr = new ofstream(filename.c_str());

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
			// in swingpyrun.h, the PyString_Check is defined to PyBytes_Check
#if PY_VERSION_HEX >= 0x03000000
			PyObject * arglist = NULL;
			PyObject * pyResult = NULL;
			if (m_mode == "b") {
				arglist = Py_BuildValue("(S)", PyBytes_FromString(str.c_str()));
				pyResult = PyEval_CallObject(m_func.func(), arglist);
			} else {
				arglist = Py_BuildValue("(s)", str.c_str());
				pyResult = PyEval_CallObject(m_func.func(), arglist);
			}
#else
			PyObject * arglist = Py_BuildValue("(s)", str.c_str());
			PyObject * pyResult = PyEval_CallObject(m_func.func(), arglist);
#endif
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
	if (output.empty()) {
		SETFLAG(m_flags, m_flagNoOutput);
		RESETFLAG(m_flags, m_flagCloseAfterUse);
		return;
	} else
		RESETFLAG(m_flags, m_flagNoOutput);

	ssize_t i;
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

	DBG_DO(DBG_UTILITY, cerr	<< "Analyzed string is " << output << endl
		                        << "Filename is " << format << endl);

	m_filename = format;
}


void closeOutput(const string & output)
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
	set(rng, seed);
}


RNG::RNG(const RNG & rhs) : m_RNG(NULL)
{
	// this will create a new instance of m_RNG.
	set(rhs.name(), rhs.seed());
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
		DBG_WARNIF(true, "advapi32.dll can not be loaded");
		return static_cast<unsigned long>(time(NULL));
	}

	/* Obtain pointers to the CryptoAPI functions
	   This will fail on some early versions of Win95 */
	pCryptAcquireContext = (CRYPTACQUIRECONTEXTA)GetProcAddress(
		hAdvAPI32, "CryptAcquireContextA");
	if (pCryptAcquireContext == NULL) {
		DBG_WARNIF(true, "Failed to get process address of CryptAcquireContextA");
		return static_cast<unsigned long>(time(NULL));
	}

	CRYPTGENRANDOM pCryptGenRandom = (CRYPTGENRANDOM)GetProcAddress(
		hAdvAPI32, "CryptGenRandom");
	if (pCryptGenRandom == NULL) {
		DBG_WARNIF(true, "Failed to get process address of CryptGenRandom");
		return static_cast<unsigned long>(time(NULL));
	}

	// Acquire context
	HCRYPTPROV hCryptProv = 0;
	if (!pCryptAcquireContext(&hCryptProv, NULL, NULL,
			PROV_RSA_FULL, CRYPT_VERIFYCONTEXT)) {
		DBG_WARNIF(true, "Can not acquire context of CryptAcquireContextA");
		return static_cast<unsigned long>(time(NULL));
	}

	// Get random data
	if (!pCryptGenRandom(hCryptProv, sizeof(seed), (unsigned char *)&seed)) {
		DBG_WARNIF(true, "Failed to get random number");
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
		size_t sz = fread(&seed, sizeof(seed), 1, devrandom);
		(void)sz; // suppress a warning.
		DBG_FAILIF(sz != 1, RuntimeError,
			"Incorrect bits of random digits are read from /dev/urandom");
		fclose(devrandom);
	} else if ((devrandom = fopen("/dev/random", "r")) != NULL) {
		size_t sz = fread(&seed, sizeof(seed), 1, devrandom);
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
void RNG::set(const char * rng, unsigned long seed)
{
	const char * rng_name = rng;

	// if RNG name is not given, try GSL_RNG_TYPE
	if (m_RNG == NULL && rng_name == NULL)
		rng_name = getenv("GSL_RNG_TYPE");

	// if a name is given ..... replace the existing RNG
	if (rng_name != NULL && rng_name[0] != '\0') {
		// locate the RNG
		const gsl_rng_type ** t, ** t0 = gsl_rng_types_setup();

		gsl_rng_default = 0;

		// check GSL_RNG_TYPE against the names of all the generators

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
			throw SystemError((boost::format("GSL_RNG_TYPE=%1% not recognized or can not generate full range (0-2^32-1) of integers.") % rng_name).str());
	} else if (m_RNG == NULL)
		// no name is given so we use a default one (mt19937)
		m_RNG = gsl_rng_alloc(gsl_rng_mt19937);

	// in the case that a name is not given, and m_RNG already exists, just set seed.

	// generate seed
	if (seed == 0)
		m_seed = generateRandomSeed();
	else
		m_seed = seed;

	// set seed
	gsl_rng_set(m_RNG, m_seed);
	m_bitByte = 0;
	m_bitIndex = 0;
}


bool RNG::randBit()
{
	if (m_bitIndex == 16)
		m_bitIndex = 0;

	if (m_bitIndex == 0)
		m_bitByte = randInt(0xFFFF);

	return (m_bitByte & (1UL << m_bitIndex++)) != 0;
}


ULONG RNG::search_poisson(UINT y, double * z, double p, double lambda)
{
	if (*z >= p) { // search to the left
		DBG_DO(DBG_UTILITY, cerr << "search to the left from " << y << endl);
		while (true) {
			DBG_FAILIF(y == 0, RuntimeError, "Zero should not be reached for a truncated Poisson distribution.");
			if ((*z = gsl_cdf_poisson_P(y - 1, lambda)) < p)
				return y;
			--y;
		}
	} else {        // search to the right
		DBG_DO(DBG_UTILITY, cerr << "search to the right from " << y << endl);
		while (true) {
			++y;
			if ((*z = gsl_cdf_poisson_P(y, lambda)) >= p)
				return y;
		}
	}
}


ULONG RNG::randTruncatedPoisson(double lambda)
{
	DBG_FAILIF(lambda == 0, ValueError, "Zero mu is not allowed for a truncated poisson distribution");

	if (lambda >= 2) {
		DBG_DO(DBG_UTILITY, cerr << "Using loop for truncated poisson with lambda = " << lambda << endl);
		// this is safe because the biggest zero-probability is Poisson(X=0,p=2)=0.13.
		// which means one out of ten times randPoisson will be called twice.
		while (true) {
			UINT numOff = randPoisson(lambda);
			if (numOff > 0)
				return numOff;
		}
	}
	// the following is adapted from R's qpois function

	// we have to use a special method that maps back from cdf.
	double p0 = gsl_ran_poisson_pdf(0, lambda);
	double p = gsl_rng_uniform(m_RNG) * (1 - p0) + p0;
	if (p + 1.01 * GSL_DBL_EPSILON >= 1.)
		return MaxRandomNumber;


	double mu = lambda;
	double sigma = sqrt(lambda);
	// gamma = sigma; PR#8058 should be kurtosis which is mu^-0.5
	double gamma = 1.0 / sigma;

	// y := approx.value (Cornish-Fisher expansion) :
	double z = gsl_cdf_ugaussian_Pinv(p);
	UINT y = static_cast<UINT>(floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5));

	if (y == 0)
		y = 1;

	z = gsl_cdf_poisson_P(y, lambda);

	DBG_DO(DBG_UTILITY, cerr	<< "Using inverse cdf with p0=" << p0 << " random quantile " << p
		                        << " and initial guess " << y << " with cdf " << z << endl);

	// fuzz to ensure left continuity; 1 - 1e-7 may lose too much :
	p *= 1 - 64 * GSL_DBL_EPSILON;

	// If the mean is not too large a simple search is OK
	return search_poisson(y, &z, p, lambda);
}


ULONG RNG::search_binomial(UINT y, double * z, double p, UINT n, double pr)
{
	if (*z >= p) {
		DBG_DO(DBG_UTILITY, cerr << "search to the left from " << y << endl);
		while (true) {
			DBG_FAILIF(y == 0, RuntimeError, "Zero should not be reached for a truncated binomial distribution");
			double newz;
			if ((newz = gsl_cdf_binomial_P(y - 1, pr, n)) < p)
				return y;
			--y;
			*z = newz;
		}
	} else {      // search to the right
		DBG_DO(DBG_UTILITY, cerr << "search to the right from " << y << endl);
		while (true) {
			++y;
			if (y >= n)
				return n;
			if ((*z = gsl_cdf_binomial_P(y, pr, n)) >= p)
				return y;
		}
	}
}


ULONG RNG::randTruncatedBinomial(UINT n, double pr)
{
	// only try once. whatever the probability is, only the successed one will be returned.
	if (n == 1)
		return 1;

	DBG_FAILIF(pr == 0. || n == 0, ValueError, "n=0 or p=0 is not allowed for truncated binomial distribution.");

	if (n * pr >= 2.) {
		// try our luck and see if we get any non-zero results ...
		int attempts = 0;
		while (++attempts < 2) {
			UINT numOff = randBinomial(n, pr);
			if (numOff > 0)
				return numOff;
		}
	}

	// the following is adapted from R's qbinom function

	double q = 1 - pr;
	if (q == 0.)
		return n;  // covers the full range of the distribution

	double p0 = gsl_ran_binomial_pdf(0, pr, n);
	double p = gsl_rng_uniform(m_RNG) * (1 - p0) + p0;
	if (p + 1.01 * GSL_DBL_EPSILON >= 1.)
		return n;

	double mu = n * pr;
	double sigma = sqrt(n * pr * q);
	double gamma = (q - pr) / sigma;

	// y := approx.value (Cornish-Fisher expansion) :
	double z = gsl_cdf_ugaussian_Pinv(p);
	UINT y = static_cast<UINT>(floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5));

	if (y == 0)
		y = 1;

	if (y > n) // way off
		y = n;

	z = gsl_cdf_binomial_P(y, pr, n);

	DBG_DO(DBG_UTILITY, cerr	<< "Using inverse cdf with p0=" << p0 << " random quantile "
		                        << p << " and initial guess " << y << " with cdf " << z << endl);

	// fuzz to ensure left continuity:
	p *= 1 - 64 * GSL_DBL_EPSILON;

	return search_binomial(y, &z, p, n, pr);
}


double pvalChiSq(double chisq, unsigned int df)
{
	return 1 - gsl_cdf_chisq_P(chisq, df);
}


void chisqTest(const vector<vectoru> & table, double & chisq, double & chisq_p)
{
	size_t m = table.size();
	size_t n = table[0].size();
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
	chisq_p = 1 - gsl_cdf_chisq_P(chisq, static_cast<double>((m - 1) * (n - 1)));
}


double armitageTrendTest(const vector<vectoru> & table, const vectorf & s)
{
	DBG_FAILIF(table.size() != 2 || table[0].size() != 3, ValueError,
		"Current Cochran-Armitage test can only handle 2 by 3 tables.");

	size_t n = table[0].size();

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
	size_t obsAA = cnt[2] > cnt[0] ? cnt[0] : cnt[2];                                             // in this algorithm, AA is rare.
	size_t obsAB = cnt[1];
	size_t obsBB = cnt[2] > cnt[0] ? cnt[2] : cnt[0];

	size_t diplotypes = obsAA + obsAB + obsBB;
	size_t rare = (obsAA * 2) + obsAB;
	size_t hets = obsAB;


	//make sure "rare" allele is really the rare allele
	if (rare > diplotypes)
		rare = 2 * diplotypes - rare;

	//make sure numbers aren't screwy
	if (hets > rare)
		throw ValueError((boost::format("HW test: %1% heterozygotes but only %2% rare alleles.") % hets % rare).str());

	vectorf tailProbs(rare + 1, 0.0);

	//start at midpoint
	//all the casting is to make sure we don't overflow ints if there are 10's of 1000's of inds
	size_t mid = (size_t)((double)rare * (double)(2 * diplotypes - rare) / (double)(2 * diplotypes));

	//check to ensure that midpoint and rare alleles have same parity
	if (((rare & 1) ^ (mid & 1)) != 0) {
		mid++;
	}
	size_t het = mid;
	size_t hom_r = (rare - mid) / 2;
	size_t hom_c = diplotypes - het - hom_r;

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
	for (het = mid; het + 2 <= rare; het += 2) {
		tailProbs[het + 2] = (tailProbs[het] * 4.0 * hom_r * hom_c) / ((het + 2.0) * (het + 1.0));
		sum += tailProbs[het + 2];
		//2 more hets for next iteration -> subtract one rare and one common homozygote
		hom_r--;
		hom_c--;
	}

	for (size_t z = 0; z < tailProbs.size(); z++)
		tailProbs[z] /= sum;

	double top = tailProbs[hets];
	for (size_t i = hets + 1; i <= rare; i++)
		top += tailProbs[i];

	double otherSide = tailProbs[hets];
	for (ssize_t i = hets - 1; i >= 0; i--)
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


string formatDescription(const string & text)
{
	vectorstr lines;
	// break from newline
	size_t pos = 0;
	size_t nextpos = 0;

	while ((nextpos = text.find('\n', pos)) != string::npos) {
		lines.push_back(text.substr(pos, nextpos - pos + 1));
		pos = nextpos + 1;
	}
	lines.push_back(text.substr(pos));
	// newtext
	string output;
	int indent = 0;
	for (size_t it = 0; it < lines.size(); ++it) {
		string line = lines[it];
		// remove leading blanks
		nextpos = line.find_first_not_of(' ');
		if (nextpos == string::npos)
			continue;
		line = line.substr(nextpos);
		string start = line.substr(0, 4);
		if (start == "<ul>") {
			++indent;
			continue;
		} else if (start == "</ul") {
			DBG_WARNIF(indent < 0, "Wrong description text caused by incorrect indent");
			if (indent != 0)
				--indent;
			continue;
		} else if (start == "<li>") {
			string indent_char("*#-.");
			string leading;
			if (indent >= 1 && static_cast<unsigned>(indent) <= indent_char.size())
				leading = string(1, indent_char[indent - 1]) + " ";
			else if (static_cast<unsigned>(indent) > indent_char.size())
				leading = ". ";
			line = leading + line.substr(4);
		} else if (start == "<ind") {
			line = "  " + line.substr(8);
		}
		// add blanks
		line = string(indent * 3, ' ') + line;
		// break lines
		if (line.size() > 78) {
			size_t lastblank = 0;
			pos = 0;
			for (size_t i = 0; i < line.size(); ++i) {
				if (line[i] == ' ')
					lastblank = i;
				if (i > 0 && i % 78 == 0) {
					output += line.substr(pos, lastblank - pos) + "\n" +
					          string(indent * 3 + (indent > 0 ? 2 : 0), ' ');
					pos = lastblank + 1;
				}
			}
			output += line.substr(pos);
		} else
			// last line (or the only line)
			output += line;
	}
	return output;
}


size_t WeightedSampler::draw()
{
	DBG_FAILIF(m_algorithm == 0, ValueError,
		"weighted sample is not initialized");
#ifdef _OPENMP
	ATOMICLONG index = 0;
#endif

	switch (m_algorithm) {
	case 1:
		// only return one
		return m_param;
	case 2:
		// all weights are the same
		return getRNG().randInt(static_cast<ULONG>(m_param));
	case 3: {
		double rN = getRNG().randUniform() * m_N;

		size_t K = static_cast<size_t>(rN);

		rN -= K;

		if (rN < m_q[K])
			return K;
		else
			return m_a[K];
	}
	case 4:
		// return according to proportion.
		if (m_index == m_sequence.size())
			m_index = 0;
#ifdef _OPENMP
		index = fetchAndIncrement(&m_index);
		return m_sequence[index];
#else

		return m_sequence[m_index++];
#endif
	default:
		throw RuntimeError("Invalid weighted sampler (empty weight?)");
	}
	// should never be reached
	return 0;
}


vectoru WeightedSampler::drawSamples(ULONG num)
{
	vectoru res(num);

	for (size_t i = 0; i < num; ++i)
		res[i] = draw();
	return res;
}


// this is used for Bernullitrials and copyGenotype
WORDTYPE g_bitMask[WORDBIT];

Bernullitrials::Bernullitrials(RNG & /* rng */)
	: m_N(0), m_prob(0), m_table(0), m_pointer(0),
	m_cur(npos)
{
}


Bernullitrials::Bernullitrials(RNG & /* rng */, const vectorf & prob, ULONG trials)
	: m_N(0), m_prob(prob), m_table(prob.size()), m_pointer(prob.size()),
	m_cur(npos)
{
	//DBG_FAILIF(trials <= 0, ValueError, "trial number can not be zero.");
	DBG_FAILIF(prob.empty(), ValueError, "probability table can not be empty.");

	if (trials == 0)
		if (*min_element(prob.begin(), prob.end()) < 0.0000001)
			m_N = 1024 * 4;
		else
			m_N = 1024;
	else
		m_N = trials;

	// initialize the table
	BitSet::iterator beg_it;
	for (size_t i = 0; i < probSize(); ++i) {
		DBG_FAILIF(m_prob[i] < 0 || m_prob[i] > 1, ValueError,
			(boost::format("Probability for a Bernulli trail should be between 0 and 1 (value %1% at index %2%)") % m_prob[i] % i).str());
		m_table[i].resize(m_N);
		beg_it = m_table[i].begin();
		m_pointer[i] = const_cast<WORDTYPE *>(BITPTR(beg_it));
	}
}


//
Bernullitrials::~Bernullitrials()
{
}


void Bernullitrials::setParameter(const vectorf & prob, size_t trials)
{
	if (trials == 0)
		if (*min_element(prob.begin(), prob.end()) < 0.0000001)
			m_N = 1024 * 4;
		else
			m_N = 1024;
	else
		m_N = trials;
	m_prob = prob;
	m_table.resize(m_prob.size());
	m_pointer.resize(m_prob.size());
	m_cur = npos;                                                             // will trigger doTrial.

	//DBG_FAILIF(trials <= 0, ValueError, "trial number can not be zero.");
	DBG_FAILIF(prob.empty(), ValueError, "probability table can not be empty.");

	BitSet::iterator beg_it;
	for (size_t i = 0; i < probSize(); ++i) {
		DBG_FAILIF(m_prob[i] < 0 || m_prob[i] > 1, ValueError,
			(boost::format("Probability for a Bernulli trail should be between 0 and 1 (value %1% at index %2%)") % m_prob[i] % i).str());
		m_table[i].resize(m_N);
		beg_it = m_table[i].begin();
		m_pointer[i] = const_cast<WORDTYPE *>(BITPTR(beg_it));
	}
}


// utility function.
void Bernullitrials::setAll(size_t idx, bool v)
{
	WORDTYPE * ptr = m_pointer[idx];

	BitSet::iterator beg_it = m_table[idx].begin();

	DBG_ASSERT(BITOFF(beg_it) == 0, SystemError, "Start of a vector<bool> is not 0");
	DBG_ASSERT(BITPTR(beg_it) == m_pointer[idx],
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
			*ptr = 0;  //~g_bitMask[rest];
	}
}


#define setBit(ptr, i)    (*((ptr) + (i) / WORDBIT) |= 1UL << ((i) - ((i) / WORDBIT) * WORDBIT))
#define unsetBit(ptr, i)  (*((ptr) + (i) / WORDBIT) &= ~(1UL << ((i) - ((i) / WORDBIT) * WORDBIT)))
// use a != 0 to avoid compiler warning
#define getBit(ptr, i)    ((*((ptr) + (i) / WORDBIT) & (1UL << ((i) - ((i) / WORDBIT) * WORDBIT))) != 0)

void Bernullitrials::doTrial()
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
					// blocks[i] = static_cast<int16_t>(getRNG().randGet());
					tmp = getRNG().randInt(0xFFFF);
					*ptr |= (0xFFFF & tmp) << (b * 16);
				}
				ptr++;
			}
			// last block
			if (rest != 0) {
				size_t b = 0;
				for (b = 0; b < rest / 16; ++b) {
					tmp = getRNG().randInt(0xFFFF);
					*ptr |= (0xFFFF & tmp) << (b * 16);
				}
				// last bits
				rest -= b * 16;
				if (rest != 0) {
					tmp = getRNG().randInt(0xFFFF);
					*ptr |= (g_bitMask[rest] & tmp) << b * 16;
				}
			}
		}
		// algorithm i Sheldon Ross' book simulation (4ed), page 54
		else if (prob < 0.5) {
			// set all to 0, then set some to 1
			setAll(cl, false);
			// it may make sense to limit the use of this method to low p,
			size_t i = 0;
			while (true) {
				// i moves at least one. (# trails until the first success)
				// 6,3 means (0 0 0 0 0 1) (0 0 1)
				ULONG step = getRNG().randGeometric(prob);
				if (step == 0)
					// gsl_ran_geometric sometimes return 0 when prob is really small.
					break;
				// i moves to 6 and 9
				i += step;
				if (i <= m_N)
					// set the 5th and 8th element to 1.
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
			size_t i = 0;
			prob = 1. - prob;
			while (true) {
				ULONG step = getRNG().randGeometric(prob);
				if (step == 0)
					// gsl_ran_geometric sometimes return 0 when prob is really small.
					break;
				i += step;
				if (i <= m_N)
					unsetBit(ptr, i - 1);
				else
					break;
			}
		}
	}
	m_cur = 0;
}


size_t Bernullitrials::curTrial()
{
	DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
	return m_cur;
}


// get a trial corresponding to m_prob.
void Bernullitrials::trial()
{
	if (m_cur == npos || m_cur == m_N - 1)  // reach the last trial
		doTrial();
	else
		m_cur++;
	DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
}


size_t Bernullitrials::trialFirstSucc(size_t idx) const
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


size_t Bernullitrials::trialNextSucc(size_t idx, size_t pos) const
{
	const BitSet & bs = m_table[idx];

	if (pos >= (m_N - 1) || m_N == 0)
		return npos;

	++pos;

	// first block
	BitSet::const_iterator it = bs.begin() + pos;
	WORDTYPE * ptr = const_cast<WORDTYPE *>(BITPTR(it));
	size_t offset = BITOFF(it);
	BitSet::const_iterator beg_it = bs.begin();
	size_t i = ptr - const_cast<WORDTYPE *>(BITPTR(beg_it));

	// mask out bits before pos
	WORDTYPE tmp = *ptr & ~g_bitMask[offset];

	size_t blk = m_N / WORDBIT;
	if (i == blk)
		// mask out bits after rest
		tmp &= g_bitMask[m_N - blk * WORDBIT];
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


void Bernullitrials::setTrialSucc(size_t idx, bool succ)
{
	DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
	if (succ)
		setBit(m_pointer[idx], m_cur);
	else
		unsetBit(m_pointer[idx], m_cur);
}


double Bernullitrials::trialSuccRate(UINT index) const
{
	// efficiency is not considered here
	size_t count = 0;

	for (size_t i = 0; i < trialSize(); ++i)
		if (getBit(m_pointer[index], i))
			count++;
	return double(count) / static_cast<double>(m_table[index].size());
}


double Bernullitrials::probSuccRate() const
{
	DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
	UINT count = 0;
	for (size_t cl = 0, clEnd = probSize(); cl < clEnd; ++cl)
		count += getBit(m_pointer[cl], m_cur) ? 1 : 0;
	return count / static_cast<double>(probSize());
}


// ###############################################

Bernullitrials_T::Bernullitrials_T(RNG & /* rng */)
	: m_N(1024), m_prob(0), m_table(0), m_pointer(0), m_cur(npos)
{
}


Bernullitrials_T::Bernullitrials_T(RNG & /* rng */, const vectorf & prob, size_t N)
	: m_N(N), m_prob(prob), m_table(N), m_pointer(N), m_cur(npos)
{
	//DBG_FAILIF(trials_T <= 0, ValueError, "trial number can not be zero.");
	DBG_FAILIF(prob.empty(), ValueError, "probability table can not be empty.");
}


//
Bernullitrials_T::~Bernullitrials_T()
{
}


void Bernullitrials_T::setParameter(const vectorf & prob, size_t N)
{
	//
	m_N = N == 0 ? 1024 : N;
	m_prob = prob;
	m_table.resize(m_N);
	m_pointer.resize(m_N);
	m_cur = npos;                                                             // will trigger doTrial.

	DBG_FAILIF(prob.empty(), ValueError, "probability table can not be empty.");

}


#define setBit(ptr, i)    (*((ptr) + (i) / WORDBIT) |= 1UL << ((i) - ((i) / WORDBIT) * WORDBIT))
#define unsetBit(ptr, i)  (*((ptr) + (i) / WORDBIT) &= ~(1UL << ((i) - ((i) / WORDBIT) * WORDBIT)))
// use a != 0 to avoid compiler warning
#define getBit(ptr, i)    ((*((ptr) + (i) / WORDBIT) & (1UL << ((i) - ((i) / WORDBIT) * WORDBIT))) != 0)

// utility function.
void Bernullitrials_T::setAll(size_t idx, bool v)
{
	if (v)
		for (size_t i = 0; i < m_N; ++i)
			setBit(m_pointer[i], idx);
	else
		for (size_t i = 0; i < m_N; ++i)
			unsetBit(m_pointer[i], idx);
}


void Bernullitrials_T::doTrial()
{
	// reset all values to 0
	BitSet::iterator beg_it;
	for (size_t i = 0; i < m_N; ++i) {
		m_table[i].clear();
		m_table[i].resize(m_prob.size(), 0);
		beg_it = m_table[i].begin();
		m_pointer[i] = const_cast<WORDTYPE *>(BITPTR(beg_it));
	}
	// for each column
	for (size_t cl = 0, clEnd = m_prob.size(); cl < clEnd; ++cl) {
		double prob = m_prob[cl];
		if (prob == 0.)
			continue;
		// algorithm i Sheldon Ross' book simulation (4ed), page 54
		else if (prob > 0 && prob < 0.2) {
			// it may make sense to limit the use of this method to low p,
			size_t i = 0;
			while (true) {
				// i moves at least one. (# trails until the first success)
				// 6,3 means (0 0 0 0 0 1) (0 0 1)
				ULONG step = getRNG().randGeometric(prob);
				if (step == 0)
					// gsl_ran_geometric sometimes return 0 when prob is really small.
					break;
				// i moves to 6 and 9
				i += step;
				if (i <= m_N)
					// set the 5th and 8th element to 1.
					setBit(m_pointer[i - 1], cl);
				else
					break;
			}
		} else if (prob == 1.) {
			setAll(cl, true);
		} else {                                                                  // 1 > m_proc[cl] > 0.5
			for (size_t i = 0; i < m_N; ++i)
				if (getRNG().randUniform() < prob)
					setBit(m_pointer[i], cl);
		}
	}
	m_cur = 0;
}


// get a trial corresponding to m_prob.
void Bernullitrials_T::trial()
{
	if (m_cur == npos || m_cur == m_N - 1)  // reach the last trial
		doTrial();
	else
		m_cur++;
	DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
}


size_t Bernullitrials_T::probFirstSucc() const
{
	size_t nProb = m_prob.size();
	size_t blk = nProb / WORDBIT;
	WORDTYPE * ptr = m_pointer[m_cur];

	size_t i = 0;

	while (i < blk && *ptr++ == 0)
		++i;

	if (i < blk)
		return i * WORDBIT + lowest_bit(*(ptr - 1));
	else {
		size_t rest = nProb - blk * WORDBIT;
		size_t tmp = *ptr & g_bitMask[rest];
		if (tmp == 0)
			return npos;
		else
			return blk * WORDBIT + lowest_bit(tmp);
	}
}


size_t Bernullitrials_T::probNextSucc(size_t pos) const
{
	const BitSet & bs = m_table[m_cur];
	size_t nProb = m_prob.size();

	if (pos >= nProb - 1 || nProb == 0)
		return npos;

	++pos;

	// first block
	BitSet::const_iterator it = bs.begin() + pos;
	WORDTYPE * ptr = const_cast<WORDTYPE *>(BITPTR(it));
	size_t offset = BITOFF(it);
	BitSet::const_iterator beg_it = bs.begin();
	size_t i = ptr - const_cast<WORDTYPE *>(BITPTR(beg_it));

	// mask out bits before pos
	WORDTYPE tmp = *ptr & ~g_bitMask[offset];

	size_t blk = nProb / WORDBIT;
	if (i == blk)
		// mask out bits after rest
		tmp &= g_bitMask[nProb - blk * WORDBIT];
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

	if (i < blk)
		return i * WORDBIT + lowest_bit(*(ptr - 1));
	else {                                                                                  // last block?
		size_t rest = nProb - blk * WORDBIT;
		// mask out bits after rest
		tmp = *ptr & g_bitMask[rest];
		if (tmp == 0)
			return npos;
		else
			return blk * WORDBIT + lowest_bit(tmp);
	}
}


void Bernullitrials_T::setTrialSucc(size_t idx, bool succ)
{
	DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
	if (succ)
		setBit(m_pointer[m_cur], idx);
	else
		unsetBit(m_pointer[m_cur], idx);
}


double Bernullitrials_T::trialSuccRate(UINT index) const
{
	UINT count = 0;

	for (size_t cl = 0, clEnd = m_N; cl < clEnd; ++cl)
		count += getBit(m_pointer[cl], index) ? 1 : 0;
	return count / static_cast<double>(m_N);
}


double Bernullitrials_T::probSuccRate() const
{
	// efficiency is not considered here
	size_t count = 0;

	for (size_t i = 0; i < m_prob.size(); ++i)
		if (getBit(m_pointer[m_cur], i))
			count++;
	return double(count) / static_cast<double>(m_prob.size());
}


#undef setBit
#undef unsetBit
#undef getBit


// Global debug and initialization related functions

void gsl_error_handler(const char * reason, const char *,
                       int, int gsl_errno)
{
	throw SystemError((boost::format("GSL Error %1%:\t %2%") % gsl_errno % reason).str());
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
		if (pbase() != pptr()) {
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
			if (pbase() == epptr()) {
				if (m_type == StdOut)
					PySys_WriteStdout("%c", c);
				else
					PySys_WriteStderr("%c", c);
			} else
				sputc(static_cast<char>(c));
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


PyObject * moduleInfo()
{
	// output a dictionary with many keys
	PyObject * dict = PyDict_New();
	PyObject * val = NULL;

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
	PyDict_SetItem(dict, PyString_FromString("revision"), val = PyInt_FromLong(num));
	Py_DECREF(val);

	// Version
#ifndef SIMUPOP_VER
	const char * ver = "snapshot";
#else
	// convert name to a string
	const char * ver = MacroQuote(SIMUPOP_VER);
#endif
	PyDict_SetItem(dict, PyString_FromString("version"), val = PyString_FromString(ver));
	Py_DECREF(val);

	// optimized
#ifdef OPTIMIZED
	PyDict_SetItem(dict, PyString_FromString("optimized"), Py_True);
#else
	PyDict_SetItem(dict, PyString_FromString("optimized"), Py_False);
#endif

	// AlleleType
#ifdef LONGALLELE
	PyDict_SetItem(dict, PyString_FromString("alleleType"), val = PyString_FromString("long"));
#elif defined BINARYALLELE
	PyDict_SetItem(dict, PyString_FromString("alleleType"), val = PyString_FromString("binary"));
#elif defined MUTANTALLELE
	PyDict_SetItem(dict, PyString_FromString("alleleType"), val = PyString_FromString("mutant"));
#elif defined LINEAGE
	PyDict_SetItem(dict, PyString_FromString("alleleType"), val = PyString_FromString("lineage"));
#else
	PyDict_SetItem(dict, PyString_FromString("alleleType"), val = PyString_FromString("short"));
#endif
	Py_DECREF(val);

#ifndef COMPILER
#  ifdef __GNUC__
#    define COMPILER "[GCC " __VERSION__ "]"
#  endif
#endif                                                                            // !COMPILER

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
	PyDict_SetItem(dict, PyString_FromString("compiler"), val = PyString_FromString(COMPILER));
	Py_DECREF(val);

	// date
	PyDict_SetItem(dict, PyString_FromString("date"), val = PyString_FromString(__DATE__));
	Py_DECREF(val);

	// version of python
	PyDict_SetItem(dict, PyString_FromString("python"), val = PyString_FromString(PY_VERSION));
	Py_DECREF(val);

	// platform
	PyDict_SetItem(dict, PyString_FromString("platform"), val = PyString_FromString(PLATFORM));
	Py_DECREF(val);

	// Number of Threads in openMP
#ifdef _OPENMP
	PyDict_SetItem(dict, PyString_FromString("threads"), val = PyLong_FromLong(numThreads()));
#else
	PyDict_SetItem(dict, PyString_FromString("threads"), val = PyLong_FromLong(0));
#endif
	Py_DECREF(val);

	// 32 or 64 bits
#ifdef _WIN64
	PyDict_SetItem(dict, PyString_FromString("wordsize"), val = PyLong_FromLong(64));
#else
#  ifdef _WIN32
	PyDict_SetItem(dict, PyString_FromString("wordsize"), val = PyLong_FromLong(32));
#  else
	PyDict_SetItem(dict, PyString_FromString("wordsize"), val = PyLong_FromLong(__WORDSIZE));
#  endif
#endif
	Py_DECREF(val);

#ifdef BINARYALLELE
	// bits for each allele
	PyDict_SetItem(dict, PyString_FromString("alleleBits"), val = PyLong_FromLong(1));
#else
	PyDict_SetItem(dict, PyString_FromString("alleleBits"), val = PyLong_FromLong(sizeof(Allele) * 8));
#endif
	Py_DECREF(val);

	// maxAllele
	PyDict_SetItem(dict, PyString_FromString("maxAllele"), val = PyLong_FromUnsignedLong(ModuleMaxAllele));
	Py_DECREF(val);

	// limits
	PyDict_SetItem(dict, PyString_FromString("maxIndex"), val = PyLong_FromUnsignedLong(static_cast<ULONG>(MaxIndexSize)));
	Py_DECREF(val);

	// debug (code)
	PyObject * codes = PyDict_New();
	for (size_t i = 0; g_debugCodes[i][0]; ++i) {
		if (g_dbgCode[i])
			PyDict_SetItemString(codes, g_debugCodes[i], Py_True);
		else
			PyDict_SetItemString(codes, g_debugCodes[i], Py_False);
	}
	PyDict_SetItem(dict, PyString_FromString("debug"), codes);
	Py_DECREF(codes);

	// availableRNGs
	PyObject * rngs = PyList_New(0);
	const gsl_rng_type ** t, ** t0;
	gsl_rng * rng;

	t0 = gsl_rng_types_setup();

	for (t = t0; *t != 0; t++) {
		rng = gsl_rng_alloc(*t);
		if (gsl_rng_min(rng) == 0 && gsl_rng_max(rng) >= MaxRandomNumber)
			PyList_Append(rngs, PyString_FromString((*t)->name));
		gsl_rng_free(rng);
	}
	PyDict_SetItem(dict, PyString_FromString("availableRNGs"), rngs);
	Py_DECREF(rngs);

	//
	return dict;
}


#ifdef BINARYALLELE

// define a good way to copy long genotype sequence
void copyGenotype(GenoIterator fr, GenoIterator to, size_t n)
{
	DBG_ASSERT(BITOFF(fr) < WORDBIT, SystemError,
		"Your vector<bool> implementation is different...");

	WORDTYPE * fr_p = const_cast<WORDTYPE *>(BITPTR(fr));
	WORDTYPE * to_p = const_cast<WORDTYPE *>(BITPTR(to));
	size_t fr_off = BITOFF(fr);
	size_t to_off = BITOFF(to);

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
			cerr	<< "Copy from " << vectora(fr, fr + n)
			        << " to " << vectora(to, to + n) << " failed " << endl;
			cerr << "Offsets are " << BITOFF(fr) << " and " << BITOFF(to) << endl;
		}
	}
#  endif
}


void clearGenotype(GenoIterator to, size_t n)
{
	WORDTYPE * to_p = const_cast<WORDTYPE *>(BITPTR(to));
	size_t to_off = BITOFF(to);

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

#ifdef MUTANTALLELE

void copyGenotype(const ConstGenoIterator begin, const ConstGenoIterator end,
                  GenoIterator it)
{
	it().copy_region(begin, end, it);
}


void clearGenotype(GenoIterator begin, GenoIterator end)
{
	begin().clear(begin.index(), end.index());
}


#endif

/* This file is used to initialize simuPOP when being load into
   python. The swig interface file will has a init% % entry to
   include this file. */
bool initialize(PyObject * module)
{
	setOptions(1);
	// tie python stdout to cerr
	std::cout.rdbuf(&g_pythonStdoutBuf);
	std::cerr.rdbuf(&g_pythonStderrBuf);

	/* Ctrl-C under windows still does not work... :-(
	   #if  defined (_WIN32) || defined (__WIN32__)
	    if (SetConsoleCtrlHandler((PHANDLER_ROUTINE)(InterruptHandler), TRUE) == FALSE)
	        cerr << "Warning: Unable to install keyboard interruption handler" << endl;
	   #endif
	 */
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

	// get population and Individual type pointer
	g_swigPopType = SWIG_TypeQuery(PopSWIGType);
	g_swigOpType = SWIG_TypeQuery(OpSWIGType);
	g_swigIndividual = SWIG_TypeQuery(IndSWIGType);
	//
	if (g_swigPopType == NULL || g_swigIndividual == NULL || g_swigOpType == NULL)
		throw SystemError("Can not get population, Individual, or Operator type pointer, your SWIG version may be wrong.");

	// load carray function and type
	if (initCustomizedTypes(module) < 0)
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


bool intList::match(ssize_t rep, const vector<bool> * activeRep) const
{
	if (m_elems.empty())
		return m_allAvail;
	vectori::const_iterator it = m_elems.begin();
	vectori::const_iterator it_end = m_elems.end();
	for (; it != it_end; ++it) {
		// positive index is easy
		if (*it >= 0) {
			if (*it == rep)
				return true;
			else
				continue;
		}
		// do not check active rep if it is not provided.
		if (activeRep == NULL)
			return true;
		// negative index
		DBG_ASSERT(!activeRep->empty() && (*activeRep)[rep], SystemError,
			"Active validation should only be done for active replicates");
		// check the simple and most used case
		if (*it == -1 && activeRep->back() && static_cast<size_t>(rep + 1) == activeRep->size())
			return true;
		// find what exactly an negative index refer to
		long cnt = -*it;
		ssize_t curRep = activeRep->size() - 1;
		for (; curRep >= 0; --curRep) {
			if ((*activeRep)[curRep])
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
			from[j] = getRNG().randInt(2) != 0;
			to[j] = 0;
		}
		size_t from_idx = getRNG().randInt(300);
		size_t to_idx = getRNG().randInt(300);
		if (from_idx > to_idx)
			continue;
		size_t length = getRNG().randInt(500);
		copyGenotype(from.begin() + from_idx,
			to.begin() + to_idx, length);
		if (vectora(from.begin() + from_idx, from.begin() + from_idx + length) !=
		    vectora(to.begin() + to_idx, to.begin() + to_idx + length)) {
			cerr	<< "Copying: " << vectora(from.begin() + from_idx, from.begin() + from_idx + length) << '\n'
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
