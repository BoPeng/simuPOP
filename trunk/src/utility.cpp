/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu                                                        *
 *                                                                         *
 *   $LastChangedDate$
 *   $Rev$
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
#if  defined(_WIN32) || defined(__WIN32__)
#include <conio.h>
#else
#include <termios.h>
#endif

#include "boost/pending/lowest_bit.hpp"
using boost::lowest_bit;

//macros to provide a portable way to make macro-passed string to C++ string
#define MacroQuote_(x) #x
#define MacroQuote(x) MacroQuote_(x)

// these functions are defined in arraymodule.c which is included
// in simuPOP_wrap.cpp
extern "C" PyObject* newcarrayobject(char* buf, char type, int size);
#ifdef SIMUMPI
#include "slave.h"
extern "C" PyObject* newcarrayiterobject(ULONG shift,
ULONG size, UINT piece_size, vectoru map);
#else
extern "C" PyObject* newcarrayiterobject(GenoIterator begin, GenoIterator end);
#endif
extern "C" bool   is_carrayobject(PyObject*);
extern "C" int    carray_length(PyObject*a);
extern "C" int    carray_itemsize(PyObject*a);
extern "C" char   carray_type(PyObject* a);
extern "C" char * carray_data(PyObject*a);
extern "C" void   initcarray(void);
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

namespace simuPOP
{

	//////////////////////////////////////////////////////////////
	/// Debug functions
	//////////////////////////////////////////////////////////////
	/// debug codes in a bitset.
	DbgBitSet g_dbgCode;

	/// debug code string. For output purpose only.
	/// debug code and DBG_CODE_LENGTH is defined in simuPOP_cfg.h
	string g_dbgString[DBG_CODE_LENGTH] =
	{
		"DBG_ALL",
		"DBG_GENERAL",
		"DBG_UTILITY",
		"DBG_OPERATOR",
		"DBG_SIMULATOR",
		"DBG_INDIVIDUAL",
		"DBG_OUTPUTER",
		"DBG_MUTATOR",
		"DBG_RECOMBINATOR",
		"DBG_INITIALIZER",
		"DBG_POPULATION",
		"DBG_STATOR",
		"DBG_TERMINATOR",
		"DBG_TAGGER",
		"DBG_VISUALIZER",
		"DBG_SELECTOR",
		"DBG_MATING",
		"DBG_MIGRATOR",
		"DBG_PROFILE",
		"DBG_MPI"
		"DBG_DEVEL"
	};

	/// set debug area, default to turn all code on
	void TurnOnDebug(DBG_CODE code)
	{
#ifndef OPTIMIZED
		if( code != DBG_ALL )
			g_dbgCode[static_cast<int>(code)] = true;
		else									  // set all
			g_dbgCode.set();
#endif
	}

	/// set debug area, default to turn all code on
	void TurnOnDebugWithName(string code)
	{
#ifndef OPTIMIZED
		for(int i=0; i < DBG_CODE_LENGTH; ++i)
		{
			if( code == g_dbgString[i])
			{
				TurnOnDebug(static_cast<DBG_CODE>(i));
				return;
			}
		}
		// this line should not have been reached
		cout << "Wrong DEBUG name " << code << endl
			<< "Please check it against the result of ListDebugCode()" << endl;
#endif
	}

	/// turn off debug, default to turn all code off
	void TurnOffDebug(DBG_CODE code)
	{
#ifndef OPTIMIZED
		if( code != DBG_ALL )
			g_dbgCode[static_cast<int>(code)] = false;
		else									  // reset all
			g_dbgCode.reset();

		if( debug( DBG_GENERAL) )
			cout << "Debug code " << g_dbgString[static_cast<int>(code)]
				<< " is turned off. cf. ListDebugCode(), TurnOnDebug()." << endl;
#else
		cout << "Debug info is ignored in optimized mode." << endl;
#endif
	}

#ifndef OPTIMIZED
	/// test if one code is turned on
	/// in DEBUG section to make sure it will not be called
	/// in optimized mode
	bool debug(DBG_CODE code)
	{
		return g_dbgCode[code];
	}
#endif

	void ListDebugCode()
	{
#ifndef OPTIMIZED
		cout << "Debug code \t On/Off" << endl;

		for(int i=0; i < DBG_CODE_LENGTH; ++i)
			cout << g_dbgString[i] << '\t' << debug(static_cast<DBG_CODE>(i)) << endl;

		cout << endl;
		if(debug(DBG_GENERAL))
			cout << "cf. TurnOnDebug(), TurnOffDebug(). " << endl;
#else
		cout << "Debug info is ignored in optimized mode." << endl;
#endif
	}

	string dbgString(DBG_CODE code)
	{
		return g_dbgString[code];
	}

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
		if(_Py_RefTotal > g_refTotal and g_refWarningCount-- > 0)
			cout << "Warning: Ref count increased from " << g_refTotal << " to " << _Py_RefTotal
				<< "\nThis may be a sign of memory leak, especially when refCount increase"
				<< "\nindefinitely in a loop. Please contact simuPOP deceloper and report"
				<< "\nthe problem.\n" << endl;
		g_refTotal = _Py_RefTotal;
	}
#endif

	//////////////////////////////////////////////////////////////
	/// Some common functions/templates
	//////////////////////////////////////////////////////////////

#if  defined(_WIN32) || defined(__WIN32__)
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
		return ((c != -1) ? 1 : 0);
	}

	int simuPOP_getch(void)
	{
		struct termios oldt, newt;
		int ch;

		tcgetattr( STDIN_FILENO, &oldt );
		newt = oldt;
		newt.c_lflag &= ~( ICANON | ECHO );
		tcsetattr( STDIN_FILENO, TCSANOW, &newt );
		ch = getchar();
		tcsetattr( STDIN_FILENO, TCSANOW, &oldt );

		return ch;
	}
#endif
}



namespace std
{
	/// how to output a dictionary
	ostream& operator<< (ostream& out, const strDict& dict)
	{
		out<<"{";
		if( ! dict.empty() )
		{
			strDict::const_iterator it=dict.begin();

			out <<  it->first << ":" << it->second;

			for(++it; it != dict.end(); ++it)
				out << ", " << it->first << ":" << it->second;
		}
		out<<"}";
		return out;
	}

	/// how to output a dictionary
	ostream& operator<< (ostream& out, const intDict& dict)
	{
		out<<"{";
		if( ! dict.empty() )
		{
			intDict::const_iterator it=dict.begin();

			out <<  it->first << ":" << it->second;

			for(++it; it != dict.end(); ++it)
				out << ", " << it->first << ":" << it->second;
		}
		out<<"}";
		return out;
	}

	/// CPPONLY: 3^n, can not use pow(3, n) because of overloading problem
	/// in msvc.
	unsigned pow3(unsigned n)
	{
		unsigned res = 1;
		for(unsigned i = 0; i < n; ++i)
			res *= 3;
		return res;
	}
}



namespace simuPOP
{

	// ///////////////////////////////////////////////////////////////
	//
	// shared variables
	//
	// ///////////////////////////////////////////////////////////////

	/// utility functions: python => C++ conversion
	/// CPPONLY
	void PyObj_As_Bool(PyObject *obj, bool& val)
	{
		if(obj==NULL || obj == Py_False || obj == Py_None)
		{
			val = false;
			return;
		}
		else if (obj == Py_True)
			val = true;
		else
			val = PyInt_AS_LONG(obj) ? true : false;
	}

	/// CPPONLY
	void PyObj_As_Int(PyObject *obj, int& val)
	{
		if(obj==NULL)
		{
			val = 0;
			return;
		}
		// try to convert
		PyObject* res = PyNumber_Int(obj);
		if( res==NULL)
			throw ValueError("Can not convert object to a integer");

		val = static_cast<int>(PyInt_AS_LONG(res));
		Py_DECREF(res);
	}

	/// CPPONLY
	void PyObj_As_Double(PyObject *obj, double& val)
	{
		if(obj==NULL)
		{
			val = 0.0;
			return;
		}
		// try to convert
		PyObject* res = PyNumber_Float(obj);
		if( res==NULL)
			throw ValueError("Can not convert object to a double number");

		val = PyFloat_AsDouble(res);
		Py_DECREF(res);
	}

	/// CPPONLY
	void PyObj_As_String(PyObject *obj, string& val)
	{
		if(obj==NULL)
		{
			val = "";
			return;
		}

		PyObject* res = PyObject_Str(obj);
		if( res==NULL)
			throw ValueError("Can not convert to a string");

		val = string(PyString_AsString(res));
		Py_DECREF(res);
	}

	/// CPPONLY
	void PyObj_As_StrDict(PyObject* obj, strDict& val)
	{
		if(obj==NULL)
		{
			val = strDict();
			return;
		}

		if( ! PyMapping_Check( obj ) )
			throw ValueError("Return value can not be converted to dictionary");

		// get key and vals
		PyObject* keys = PyMapping_Keys(obj);
		PyObject* vals = PyMapping_Values(obj);

		// number of items
		UINT sz = PyList_Size(keys);

		val.clear();

		// assign values
		for(size_t i=0; i<sz; ++i)
		{
			string k;
			double v;
			PyObj_As_String( PyList_GetItem(keys, i), k);
			PyObj_As_Double( PyList_GetItem(vals, i), v);

			val.insert( strDict::value_type(k, v));
		}
	}

	/// CPPONLY
	void PyObj_As_Array(PyObject* obj, vectorf& val)
	{
		if(obj==NULL)
		{
			val = vectorf();
			return;
		}
		DBG_ASSERT( PySequence_Check(obj), ValueError, "PyObj_As_Array: Expecting a sequence");

		val.resize( PySequence_Size( obj ));

		// assign values
		for(size_t i=0, iEnd=val.size(); i<iEnd; ++i)
		{
			PyObject* item = PySequence_GetItem(obj, i);
			PyObj_As_Double(item, val[i]);
			Py_DECREF(item);
		}
	}

	/// CPPONLY
	void PyObj_As_IntArray(PyObject* obj, vectori& val)
	{
		if(obj==NULL)
		{
			val = vectori();
			return;
		}
		DBG_ASSERT( PySequence_Check(obj), ValueError, "PyObj_As_IntArray: Expecting a sequence");

		val.resize( PySequence_Size( obj ));

		// assign values
		for(size_t i=0, iEnd=val.size(); i<iEnd; ++i)
		{
			PyObject * item = PySequence_GetItem(obj, i);
			PyObj_As_Int(item, val[i]);
			Py_DECREF(item);
		}
	}

	/// CPPONLY
	void PyObj_As_IntDict(PyObject* obj, intDict& val)
	{
		if(obj==NULL)
		{
			val = intDict();
			return;
		}

		if( ! PyMapping_Check( obj ) )
			throw ValueError("Return value can not be converted to dictionary");

		// get key and vals
		PyObject* keys = PyMapping_Keys(obj);
		PyObject* vals = PyMapping_Values(obj);

		// number of items
		UINT sz = PyList_Size(keys);

		val.clear();

		// assign values
		for(size_t i=0; i<sz; ++i)
		{
			int k;
			double v;
			PyObj_As_Int( PyList_GetItem(keys, i), k);
			PyObj_As_Double( PyList_GetItem(vals, i), v);

			val.insert( intDict::value_type(k, v));
		}
		Py_XDECREF(keys);
		Py_XDECREF(vals);
	}

	/// CPPONLY
	bool PyObj_Is_IntNumArray(PyObject * obj)
	{
		return is_carrayobject(obj) &&
			carray_type(obj) == 'i' ;
	}

	/// CPPONLY
	bool PyObj_Is_DoubleNumArray(PyObject * obj)
	{
		return is_carrayobject(obj) &&
			carray_type(obj) == 'd';
	}

	/// CPPONLY
	bool PyObj_Is_AlleleNumArray(PyObject * obj)
	{
		return is_carrayobject(obj) && carray_type(obj) == 'a';
	}

	/// CPPONLY
	PyObject* Int_Vec_As_NumArray(vectori::iterator begin, vectori::iterator end)
	{
		PyObject* res = newcarrayobject(reinterpret_cast<char*>(&*begin), 'i', end-begin);

		DBG_FAILIF(res==NULL, ValueError, "Can not convert buf to int num array");
		return res;
	}

	/// CPPONLY
	PyObject* Double_Vec_As_NumArray(vectorf::iterator begin, vectorf::iterator end)
	{
		PyObject* res = newcarrayobject(reinterpret_cast<char*>(&*begin), 'd', end-begin);

		DBG_FAILIF(res==NULL, ValueError, "Can not convert buf to double num array");
		return res;
	}

#ifdef SIMUMPI
	/// CPPONLY
	PyObject* Allele_Vec_As_NumArray(ULONG shift,
		ULONG size, UINT piece_size, vectoru map)
	{
		PyObject* res = newcarrayiterobject(shift, size, piece_size, map);
		DBG_FAILIF(res==NULL, ValueError, "Can not convert buf to Allele num array");
		return res;
	}

#else

	/// CPPONLY
	PyObject* Allele_Vec_As_NumArray(GenoIterator begin, GenoIterator end)
	{
		PyObject* res = newcarrayiterobject(begin, end);
		DBG_FAILIF(res==NULL, ValueError, "Can not convert buf to Allele num array");
		return res;
	}
#endif

	/// CPPONLY
	PyObject* Info_Vec_As_NumArray(InfoIterator begin, InfoIterator end)
	{
		PyObject* res = newcarrayobject(reinterpret_cast<char*>(&*begin), 'd', end-begin);
		DBG_FAILIF(res==NULL, ValueError, "Can not convert buf to info num array");
		return res;
	}

	/// CPPONLY
	int NumArray_Size(PyObject* obj)
	{
		// return PyArray_Size(obj);
		return carray_length(obj);
	}

	/// CPPONLY
	char* NumArray_Data(PyObject* obj)
	{
		// return reinterpret_cast<PyArrayObject*>(obj)->data;
		return carray_data(obj);
	}

	// //////////////////////////////////////////////////////

	// copy constructor
	/// CPPONLY
	SharedVariables::SharedVariables(const SharedVariables& rhs)
		: m_ownVars(rhs.m_ownVars)
	{
		if(rhs.m_ownVars)
		{
			m_dict = PyDict_New();
			PyObject *key = 0, *value = 0;

			Py_ssize_t i=0;
			while(PyDict_Next(rhs.m_dict, &i, &key, &value))
			{
				//Py_INCREF(key);
				//Py_INCREF(value);
				PyDict_SetItem(m_dict, key,value);
			}
		}
		else
			m_dict = rhs.m_dict;
	}

	/// CPPONLY
	SharedVariables::~SharedVariables()
	{
		if(m_ownVars)
		{
			PyDict_Clear(m_dict);
			Py_DECREF(m_dict);
		}
	}

	/// setvars C++ ==> Python
	/// CPPONLY
	PyObject* SharedVariables::setVar(const string& name, const PyObject* val)
	{
		/* find the first piece */
		size_t i, s;
		for( i=0; i<name.size() && name[i] != '[' && name[i] != '{'; ++i);

		if( i==0 )
			throw ValueError("Empty name? " + name);

		// we need to keep current parent, key/index, type (0 for array, 1 for dict)
		// keytype (if curType = 1) (0 for string, 1 for num)
		int curType = 1;
		PyObject* curParent = m_dict;			  // should be always valid
		PyObject* curKey = PyString_FromString(const_cast<char*>(name.substr(0, i).c_str()));
		int curIdx = 0;
		PyObject * curChild = NULL;

		next:
		// get par[1] (dict), curChild can be null, or borrow ref
		if( curType == 1 )
			curChild = PyDict_GetItem(curParent, curKey);
		// get par[1] (list)
		else
			curChild = PyList_GetItem(curParent, curIdx);

		// subPop[0]{'alleleNum'}[0]
		// curParent:  m_dict
		// curChild:   m_dirt['subPop']

		if( name[i] == '[' )					  // we need to put an array in ...
		{
			// m_dict['subPop'][0]

			// look for index
			s = ++i;
			for( ; name[i] != ']'; i++)
				if( !isdigit(name[i]) )
					throw ValueError("Expecting numbers: " + name);

			// get index
			size_t idx = atoi( name.substr(s,i-s).c_str());
			// not exist
			//
			// create subPop, ...
			if( curChild == NULL || ! PyList_Check(curChild) )
			{
				// Py_XDECREF(curChild);
				// create a sequence with this name
				curChild = PyList_New(idx+1);
				// keep this
				// Py_INCREF(curChild);

				if(curChild == NULL)
					throw  ValueError("Can not create list: " + name);

				for(size_t a = 0; a < idx; ++a )
				{
					Py_INCREF(Py_None);
					PyList_SetItem(curChild, a, Py_None);
				}
				// now set that one
				if( curType == 1 )				  // dictionary?
				{
					PyDict_SetItem(curParent, curKey, curChild);
					Py_XDECREF(curChild);
				}
				else
					// PyList_SetItem steals ref
					PyList_SetItem(curParent, curIdx, curChild);
			}
			else								  // if exist, if it is a sequence?
			{
				size_t curSize = PyList_Size(curChild);
				/// if the size if enough, get the item
				if( curSize <= idx )
				{
					while(curSize++ <= idx)
					{
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
		}
		else if( name[i] == '{' )				  // dictionary
		{
			//	curParent = subPop[0]
			//  curChild = subPop[0]{'alleleNum'}

			bool numKey;

			// look for index,
			s = ++i;
			if(name[s] == '\'' || name[s] == '\"')
				numKey = false;
			else
				numKey = true;

			for( ; name[i] != '}' && i < name.size(); ++i);

			assert(name[i] == '}');

			PyObject * childKey;

			if(numKey)
				childKey = PyInt_FromString( const_cast<char*>(name.substr(s,i-s).c_str()), NULL, 0);
			else
				childKey = PyString_FromString(const_cast<char*>(name.substr(s+1, i-s-2).c_str()));
												  // not exist
			if(curChild == NULL || !PyDict_Check(curChild) )
			{
				// create dictionary with given key, append
				curChild = PyDict_New();
				// keep this

				if( curChild == NULL)
					throw  ValueError("Can not create dictionary" + name);

				// now set that one
				if( curType == 1 )				  // dictionary?
				{
					PyDict_SetItem(curParent, curKey, curChild);
					Py_XDECREF(curChild);
				}
				else
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
		}
		else									  // end here
		{
			// CASE:  'res' = 1
			//
			// regardless curChild == NULL
			// set value in
			if( curType == 1 )					  // dictionary?
			{
				PyDict_SetItem(curParent, curKey, const_cast<PyObject*>(val));
				Py_XDECREF( const_cast<PyObject*>(val));
			}
			else
				PyList_SetItem(curParent, curIdx, const_cast<PyObject*>(val));
			Py_XDECREF(curKey);
		}
		return const_cast<PyObject*>(val);
	}

	/// CPPONLY
	PyObject* SharedVariables::getVar(const string& name, bool nameError)
	{
		DBG_ASSERT( m_dict != NULL, ValueError,
			"Shared variables are not associated with any Python variable. You populaiton might not be part of a simulator.");

		// go deep in to the string
		size_t i, s;
		for( i=0; i<name.size() && name[i] != '[' && name[i] != '{'; ++i);

		if( i==0 )
			throw ValueError("Empty name? " + name);

		int curType = 1;
		PyObject* curParent = m_dict;			  // should be always valid
		PyObject* curKey = PyString_FromString(const_cast<char*>(name.substr(0, i).c_str()));
		int curIdx = 0;
		PyObject * curChild;

		next:
		if( curType == 1 )
			curChild = PyDict_GetItem(curParent, curKey);
		else
			curChild = PyList_GetItem(curParent, curIdx);

		if(curChild == NULL)
		{
			Py_XDECREF(curKey);
			// name does not exist
			if( nameError)
				throw ValueError("Shared variable name '" + name + "' does not exist.");
			else
				return NULL;
		}

		if( name[i] == '[' )					  // we need to put an array in ...
		{
			// get key
			// look for index
			s = ++i;
			for( ; name[i] != ']'; i++)
				if( !isdigit(name[i]) )
					throw ValueError("Expecting numbers: " + name);

			// get index
			int idx = atoi( name.substr(s,i-s).c_str());

			curType = 0;
			curParent = curChild;
			Py_XDECREF(curKey);
			curKey = NULL;
			curIdx = idx;
			i++;
			goto next;
		}
		else if (name[i] == '{' )
		{
			bool numKey;

			// look for index,
			s = ++i;
			if(name[s] == '\'' || name[s] == '\"')
				numKey = false;
			else
				numKey = true;

			for( ; name[i] != '}' && i < name.size(); ++i);

			assert(name[i] == '}');

			PyObject * childKey;

			if(numKey)
				childKey = PyInt_FromString( const_cast<char*>(name.substr(s,i-s).c_str()), NULL, 0);
			else
				childKey = PyString_FromString(const_cast<char*>(name.substr(s+1, i-s-2).c_str()));

			// ready for iteration
			curType = 1;
			curParent = curChild;
			Py_XDECREF(curKey);
			curKey = childKey;
			++i;
			goto next;
		}
		else
		{
			Py_XDECREF(curKey);
			return curChild;
		}
	}

	/// CPPONLY
	void SharedVariables::removeVar(const string& name)
	{
		DBG_ASSERT( m_dict != NULL, ValueError,
			"Shared variables are not associated with any Python variable. You populaiton might not be part of a simulator.");

		DBG_DO(DBG_UTILITY, cout << "Removing variable " << name << endl);

		// go deep in to the string
		size_t i, s;
		for( i=0; i<name.size() && name[i] != '[' && name[i] != '{'; ++i);

		if( i==0 )
			throw ValueError("Empty name? " + name);

		int curType = 1;
		PyObject* curParent = m_dict;			  // should be always valid
		PyObject* curKey = PyString_FromString(const_cast<char*>(name.substr(0, i).c_str()));
		int curIdx = 0;
		PyObject * curChild;

		next:
		if( curType == 1 )
			curChild = PyDict_GetItem(curParent, curKey);
		else
			curChild = PyList_GetItem(curParent, curIdx);

		// maybe the item has been removed?
		if(curChild == NULL)
		{
			Py_XDECREF(curKey);
			return;
		}

		if( name[i] == '[' )					  // we need to put an array in ...
		{
			// get key
			// look for index
			s = ++i;
			for( ; name[i] != ']'; i++)
				if( !isdigit(name[i]) )
					throw ValueError("Expecting numbers: " + name);

			// get index
			int idx = atoi( name.substr(s,i-s).c_str());

			curType = 0;
			curParent = curChild;
			Py_XDECREF(curKey);
			curKey = NULL;
			curIdx = idx;
			i++;
			goto next;
		}
		else if (name[i] == '{' )
		{
			bool numKey;

			// look for index,
			s = ++i;
			if(name[s] == '\'' || name[s] == '\"')
				numKey = false;
			else
				numKey = true;

			for( ; name[i] != '}' && i < name.size(); ++i);

			assert(name[i] == '}');

			PyObject * childKey;

			if(numKey)
				childKey = PyInt_FromString( const_cast<char*>(name.substr(s,i-s).c_str()), NULL, 0);
			else
				childKey = PyString_FromString(const_cast<char*>(name.substr(s+1, i-s-2).c_str()));

			// ready for iteration
			curType = 1;
			curParent = curChild;
			Py_XDECREF(curKey);
			curKey = childKey;
			++i;
			goto next;
		}
		else									  // now the last piece, we know the parents
		{
			if( curType == 1 )
				PyDict_DelItem(curParent, curKey);
			else
			{
				// to avoid affecting order of others, set to None
				Py_INCREF(Py_None);
				PyList_SetItem(curParent, curIdx, Py_None);
			}
			Py_XDECREF(curKey);
			return;
		}

	}

	PyObject* SharedVariables::setBoolVar(const string& name, const bool val)
	{
		PyObject *obj = val ? Py_True : Py_False;
		Py_INCREF(obj);
		return setVar(name, obj);
	}

	PyObject* SharedVariables::setIntVar(const string& name, const int val)
	{
		return setVar(name, PyInt_FromLong(val));
	}

	PyObject* SharedVariables::setDoubleVar(const string& name, const double val)
	{
		return setVar(name, PyFloat_FromDouble(val));
	}

	PyObject* SharedVariables::setStringVar(const string& name, const string & val)
	{
		return setVar(name, Py_BuildValue("s", name.c_str()));
	}

	PyObject* SharedVariables::setIntVectorVar(const string& name, const vectori& val)
	{
		PyObject * obj = PyList_New(0);
		PyObject * item;
		for(vectori::const_iterator it=val.begin();
			it < val.end(); ++it)
		{
			item = PyInt_FromLong(*it);
			PyList_Append(obj, item);
			Py_XDECREF(item);
		}
		return setVar(name, obj);
	}

	///CPPONLY
	PyObject* SharedVariables::setDoubleVectorVar(const string& name, const vectorf& val)
	{
		PyObject * obj = PyList_New(0);
		PyObject * item;
		for(vectorf::const_iterator it=val.begin();
			it < val.end(); ++it)
		{
			item = PyFloat_FromDouble(*it);
			PyList_Append(obj, item);
			Py_XDECREF(item);
		}
		return setVar(name, obj);
	}

	PyObject* SharedVariables::setStrDictVar(const string& name, const strDict & val)
	{
		PyObject *obj = PyDict_New();
		PyObject * v;
		for (strDict::const_iterator i= val.begin(); i!= val.end(); ++i)
		{
			PyDict_SetItemString(obj, const_cast<char*>(i->first.c_str()),
				v=PyFloat_FromDouble(i->second));
			Py_XDECREF(v);
		}
		return setVar(name, obj);
	}

	PyObject* SharedVariables::setIntDictVar(const string& name, const intDict & val)
	{
		PyObject *obj = PyDict_New();
		PyObject * u, *v;
		for (intDict::const_iterator i= val.begin(); i!= val.end(); ++i)
		{
			PyDict_SetItem(obj,
				u=PyInt_FromLong( i->first ),
				v=PyFloat_FromDouble(i->second));
			Py_XDECREF(u);
			Py_XDECREF(v);
		}
		return setVar(name, obj);
	}

	/// CPPONLY
	void save_none(string& str)
	{
		str += "n";
	}

	/// CPPONLY
	PyObject* load_none(const string& str, size_t& offset)
	{
		offset++;
		Py_INCREF(Py_None);
		return Py_None;
	}

	/// CPPONLY
	void save_int(string& str, PyObject * args)
	{
		long l = PyInt_AS_LONG((PyIntObject *)args);
		// type + string + ' '
		str += 'i' + toStr(l) + ' ';
	}

	/// CPPONLY
	PyObject* load_int(const string& str, size_t& offset)
	{
		// search for blank
		int len = 0;
		while( str[offset+len+1] != ' ') len++;
		PyObject* val = PyInt_FromString(
			const_cast<char*>(str.substr(offset+1, len).c_str()),
			NULL, 0);
		offset += len + 2;
		return val;
	}

	/// CPPONLY
	void save_long(string& str, PyObject * args)
	{
		long l = PyInt_AS_LONG(args);
		// type +  string + ' '
		str += 'i' + toStr(l) + ' ';
	}

	/// CPPONLY
	PyObject* load_long(const string& str, size_t& offset)
	{
		int len = 0;
		while( str[offset+len+1] != ' ') len++;
		PyObject* val = PyLong_FromString(
			const_cast<char*>(str.substr(offset+1, len).c_str()),
			NULL, 0);
		offset += len + 2;
		return val;
	}

	/// CPPONLY
	void save_float(string& str, PyObject * args)
	{
		double d = PyFloat_AsDouble(args);
		// type + string
		str += 'f' + toStr(d) + ' ';
	}

	/// CPPONLY
	PyObject* load_float(const string& str, size_t& offset)
	{
		int len = 0;
		while( str[offset +len+1] != ' ') len++;
		double d = atof( const_cast<char*>(str.substr(offset+1, len).c_str()));
		offset += len + 2;
		return PyFloat_FromDouble(d);
	}

	/// CPPONLY
	void save_string(string& str, PyObject * args)
	{
		char * s = PyString_AsString(args);
		str += 's' + string(s)+'\0';
	}

	/// CPPONLY
	PyObject* load_string(const string& str, size_t& offset)
	{
		char * s = const_cast<char*>(str.c_str()) + offset + 1;
		size_t len = strlen(s);
		offset += len + 2;
		DBG_ASSERT( str[offset-1] == '\0', SystemError,
			"Error while loading string from str");

		return PyString_FromString(s);
	}

	/// CPPONLY
	void saveObj(string& str, PyObject* args);
	/// CPPONLY
	PyObject* loadObj(const string& vars, size_t& offset);

	/// CPPONLY
	void save_dict(string& str, PyObject* args)
	{
		PyObject *key = 0, *value = 0;

		str += 'd';								  // dictionary
		Py_ssize_t i = 0;
		while (PyDict_Next(args, &i, &key, &value))
		{
			saveObj(str, key);
			saveObj(str, value);
		}
		str += 'e';								  // ending of a dictionary
	}

	/// CPPONLY
	PyObject* load_dict(const string& vars, size_t& offset)
	{
		// skip 'd'
		offset++;
		PyObject * d = PyDict_New();
		while( vars[offset] != 'e' )
		{
			// get key
			PyObject* key = loadObj(vars, offset);
			// get value
			PyObject* val = loadObj(vars, offset);
			PyDict_SetItem(d, key, val);
		}
		offset++;								  // skip 'e'
		return d;
	}

	/// CPPONLY
	void save_list(string& str, PyObject* args)
	{
		str += 'L';								  // dictionary
		int len = PyList_Size(args);
		for (int i = 0; i < len; i++)
		{
			PyObject* elem = PyList_GET_ITEM((PyListObject *)args, i);
			saveObj(str, elem);
		}

		str += 'e';								  // ending of a dictionary
	}

	/// CPPONLY
	PyObject* load_list(const string& vars, size_t& offset)
	{
		// skip 'L'
		offset++;
		PyObject * d = PyList_New(0);
		while( vars[offset] != 'e' )
		{
			PyObject* elem = loadObj(vars, offset);
			PyList_Append(d, elem);
		}
		offset++;								  // skip 'e'
		return d;
	}

	/// CPPONLY
	void save_tuple(string& str, PyObject* args)
	{
		str += 't';								  // dictionary
		int len = PyTuple_Size(args);
		// save length
		str += toStr(len) + ' ';
		// save items
		for (int i = 0; i < len; i++)
		{
			PyObject* elem = PyTuple_GET_ITEM((PyTupleObject *)args, i);
			saveObj(str, elem);
		}
	}

	/// CPPONLY
	PyObject* load_tuple(const string& vars, size_t& offset)
	{
		// skip 't'
		offset++;
		// search for blank
		int l = 0;
		while( vars[offset+l] != ' ') l++;
		int len = atoi(vars.substr(offset, l).c_str());
		offset += l + 1;
		PyObject * d = PyTuple_New(len);
		for(int i=0; i<len; ++i)
		{
			PyObject* elem = loadObj(vars, offset);
			PyTuple_SET_ITEM(d, i, elem);
		}
		return d;
	}

	// can not save or load binary carray.
	/*
	void save_carray(string& str, PyObject* args)
	{
	  DBG_ASSERT(carray_type(args) != 'a', ValueError,
		"Binary carray can not be saved");

	  char* d = carray_data(args);
	  size_t len = carray_length(args);
	  str += 'a' + toStr(len) + ' ' + carray_type(args)
		+ string(d, d+len*carray_itemsize(args));
	}

	PyObject* load_carray(const string& vars, size_t& offset)
	{
	// skip 'a'
	size_t lenlen = 0;
	while( vars[offset+lenlen+1] != ' ' ) lenlen ++;
	size_t len = atoi( const_cast<char*>(vars.substr(offset+1, lenlen).c_str()));
	// get type
	offset += lenlen+2;
	char type = vars[offset];
	char * ptr = const_cast<char*>(vars.c_str()) + offset + 1;
	PyObject* arr = newcarrayobject(type, len, ptr, true);
	offset += len*carray_itemsize(arr) + 1;
	return arr;
	}
	*/

	/// CPPONLY
	void saveObj(string& str, PyObject* args)
	{
		PyTypeObject *type;

		if(args == Py_None)
		{
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
		else
		{
			// some other unknown type
			DBG_DO(DBG_UTILITY, cout << "Warning: object of type '" + toStr(type->tp_name) + "' cannot be saved. Use none.");
			save_none(str);
		}
	}

	/// CPPONLY
	PyObject* loadObj(const string& vars, size_t& offset)
	{
		switch( vars[offset])
		{
			case 'd':							  //
				return load_dict(vars, offset);
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
				DBG_DO(DBG_UTILITY, cout << endl);
				throw ValueError("Unknown type code at offset " + toStr(offset));
			}
		}
	}

	/// CPPONLY
	string SharedVariables::asString() const
	{
		// go through each variable and save
		string str;
		saveObj(str, m_dict);
		str += 'e';								  // ending character
		return str;
	}

	/// CPPONLY
	void SharedVariables::fromString(const string& vars)
	{
		size_t offset = 0;
		PyObject * obj;
		obj = loadObj(vars, offset);

		DBG_ASSERT( obj == NULL || vars[offset] == 'e', SystemError,
			"Failed to load objects from string");

		// remove m_dict
		if( m_ownVars )
		{
			PyDict_Clear(m_dict);
			Py_XDECREF(m_dict);
		}
		m_ownVars = true;
		m_dict = obj;
	}

	// //////////// END OF SharedVars ///////

	// //////////// BEGIN of SimuVars
	// simuVars will hold replicate specific Shared Variables.

	/// global dictionary
	/// will be set by // initialize() function
	/// DO NOT OWN the dictionaries
	SharedVariables g_main_vars, g_module_vars;
	swig_type_info* g_swigPopType, *g_swigindividual;

	SharedVariables& mainVars()
	{
		return g_main_vars;
	}

	SharedVariables& moduleVars()
	{
		return g_module_vars;
	}

	PyObject* pyPopObj(void*p)
	{
		return SWIG_NewPointerObj(p, g_swigPopType, 0);
	}

	PyObject* pyIndObj(void*p)
	{
		return SWIG_NewPointerObj(p, g_swigindividual, 0);
	}

	// ////////////////////////////////////////////////////////////
	// Expression evaluation
	// ////////////////////////////////////////////////////////////
	// because of ref count, need to define copier

	/// CPPONLY
	Expression::Expression(const Expression& rhs)
		: m_expr(rhs.m_expr), m_stmts(rhs.m_stmts), m_locals(rhs.m_locals)
	{
		if(m_expr != NULL)
			Py_INCREF(m_expr);
		if(m_stmts != NULL)
			Py_INCREF(m_stmts);
	}

	/// CPPONLY
	Expression::~Expression()
	{
		// release compiled code
		if( m_expr != NULL )
			Py_DECREF(m_expr);
		if( m_stmts != NULL )
			Py_DECREF(m_stmts);
	}

	/// CPPONLY
	void Expression::compileExpr(const string& expr)
	{
		if( m_expr != NULL)
		{
			// discard previous
			Py_XDECREF(m_expr);
			m_expr = NULL;
		}

		if(expr.empty())
			return;

		m_expr = Py_CompileString(const_cast<char*>(expr.c_str()), "<embed>",
			Py_eval_input);

		if(m_expr == NULL)
			throw ValueError("Expression '" + expr + "' is not a valid python expression.");
	}

	/// CPPONLY
	void Expression::compileStmts(const string& stmts)
	{
		if( m_stmts != NULL)
		{
			// discard previous statements
			Py_XDECREF(m_stmts);
			m_stmts = NULL;
		}

		if( stmts.empty())
			return;

		// clear leading blank etc for statements ? Not yet.
		m_stmts = Py_CompileString(const_cast<char*>(stmts.c_str()), "<embed>",
			Py_file_input);

		if( m_stmts == NULL)
			throw ValueError("statement '" + stmts + "' is not valid.");
	}

	/// python expression
	PyObject* Expression::evaluate()
	{
		if( m_expr == NULL && m_stmts == NULL )
			return NULL;

		DBG_ASSERT( mainVars().dict() != NULL && m_locals != NULL,
			ValueError, "Can not evalulate. Dictionary is empty!");

		PyObject * res=NULL;
		if( m_stmts != NULL)
		{
			res = PyEval_EvalCode((PyCodeObject*) m_stmts,
				mainVars().dict(), m_locals );
			if(res == NULL)
			{
				PyErr_Print();
				throw SystemError("Evalulation of statements failed");
			}
			else
			{
				Py_DECREF(res);
				res = NULL;
			}
		}

		if( m_expr != NULL)
		{
			res = PyEval_EvalCode((PyCodeObject*) m_expr,
				mainVars().dict(), m_locals );
			if(res  == NULL)
			{
				PyErr_Print();
				throw SystemError("Evalulation of expression failed");
			}
		}

		return res;
	}

	/// CPPONLY
	bool Expression::valueAsBool()
	{
		PyObject* res = evaluate();
		if( res==NULL )
			return false;
		bool val;
		PyObj_As_Bool(res, val);
		Py_XDECREF(res);
		return val;
	}

	/// CPPONLY
	int Expression::valueAsInt()
	{
		PyObject* res = evaluate();
		if( res==NULL )
			return 0;
		int val;
		PyObj_As_Int(res, val);
		Py_XDECREF(res);
		return val;
	}

	/// CPPONLY
	string Expression::valueAsString()
	{
		PyObject* res = evaluate();
		if( res==NULL )
			return "";
		string val;
		PyObj_As_String(res, val);
		Py_XDECREF(res);
		return val;
	}

	/// CPPONLY
	vectorf Expression::valueAsArray()
	{
		PyObject* res = evaluate();
		if( res==NULL )
			return vectorf();
		vectorf val;
		PyObj_As_Array(res, val);
		Py_XDECREF(res);
		return val;
	}

	/// CPPONLY
	strDict Expression::valueAsStrDict()
	{
		PyObject* res = evaluate();
		if( res==NULL )
			return strDict();
		strDict val;
		PyObj_As_StrDict(res, val);
		Py_XDECREF(res);
		return val;
	}

	/// CPPONLY
	intDict Expression::valueAsIntDict()
	{
		PyObject* res = evaluate();
		if( res==NULL )
			return intDict();
		intDict val;
		PyObj_As_IntDict(res, val);
		Py_XDECREF(res);
		return val;
	}

	//////////////////////////////////////////////////////////////
	/// Stream element, can be of different types
	//////////////////////////////////////////////////////////////

	StreamElem::StreamElem(const string& name, bool readable, bool realAppend, bool useString )
		:m_filename(name)
	{

		DBG_DO(DBG_UTILITY, cout << "creating " << name << " with parameter " <<
			readable << " " << realAppend << " " << useString << " " << endl );

		if( useString )
		{
			// existing file will be truncated...
			m_stream = new stringstream(std::ios::in | std::ios::out);
			m_type = SSTREAM;
			m_append = false;
		}
		else									  // no string, usual file
		{
			m_append = realAppend;
			if( readable )
			{
				m_type = FSTREAM;
				if (realAppend )				  // readable , append
					m_stream = new fstream(name.c_str(),  std::ios::in | std::ios::out | std::ios::ate );
				else							  // readable, ! append
					// existing file will be truncated...
					m_stream = new fstream(name.c_str(),  std::ios::in | std::ios::trunc | std::ios::out);
			}
			else
			{
				m_type = OFSTREAM;
				if( realAppend)					  // ! readable, append
					m_stream = new ofstream(name.c_str(), std::ios::out | std::ios::app );
				else							  //  !readable, !append )
					// existing file will be truncated...
					m_stream = new fstream(name.c_str(),  std::ios::trunc | std::ios::out);
			}
		}

		if( m_stream == NULL || !*m_stream )
			throw ValueError("Can not open specified file:" + name);

		DBG_DO(DBG_UTILITY, cout << "New file info: " << info() << endl );

	}

	/// copy constructor, we need to clear rhs.m_stream to avoid closing file too early
	/// this is techniquely advanced (and dangerous)
	StreamElem::StreamElem(const StreamElem& rhs)
	{
		m_filename = rhs.m_filename;
		m_type = rhs.m_type;
		m_append = rhs.m_append;
		m_stream = rhs.m_stream;
		const_cast<StreamElem&>(rhs).m_stream = NULL;
	}

	/// destructor
	StreamElem::~StreamElem()
	{
		// close file.
		if( m_stream != NULL)
		{

			DBG_DO(DBG_UTILITY, cout << "Closing file " << m_filename << endl);

			if( m_type == OFSTREAM)
				static_cast<ofstream*>(m_stream)->close();
			else if( m_type ==  FSTREAM)
				static_cast<fstream*>(m_stream)->close();

			// do not do anything for fstream (no close function

			delete m_stream;
		}
	}

	/// the file was write-only, re-open it as read-write
	void StreamElem::makeReadable()
	{

		DBG_DO(DBG_UTILITY, cout << "File was opened write-only. Re-open it.  " << info()  <<  endl);

		static_cast<ofstream*>(m_stream)->close();

		// have to re-create a stream since pointer type is different.
		delete m_stream;

		// try to keep file content
		m_stream = new fstream( m_filename.c_str(),  std::ios::in | std::ios::out | std::ios::ate );

		if( m_stream == NULL || !*m_stream )
			throw ValueError("Can not re-open specified file.");

		m_type = FSTREAM;
	}

	/// change the append status
	void StreamElem::makeAppend(bool append)
	{

		DBG_DO(DBG_UTILITY, cout << "File append status changes to " << append <<  endl);

		DBG_FAILIF( m_type == SSTREAM , ValueError, "String stream can not be maded appendable. ");

		m_append = append;

		// if append = true: no problem, keep writing just do not close file. otherwise
		// reopen file. Otherwise,
		if( !append )
		{

			DBG_DO(DBG_UTILITY, cout << "Re-open the file " <<  endl);

			// re-open the file.

			if( m_type == FSTREAM )
			{
				static_cast<fstream*>(m_stream)->close();
				static_cast<fstream*>(m_stream)->open( m_filename.c_str(), std::ios::in | std::ios::out | std::ios::trunc );
			}
			else if (m_type == OFSTREAM )
			{
				static_cast< ofstream*>(m_stream)->close();
				static_cast< ofstream*>(m_stream)->open( m_filename.c_str(), std::ios::out | std::ios::trunc);
			}
		}
	}

	string StreamElem::info()
	{
		ostringstream out;
		switch( m_type )
		{
			case OFSTREAM:
				out << m_filename << " : write only file stream. " << endl;

				DBG_DO(DBG_UTILITY, out << "(write pos: " << static_cast<ofstream*>(m_stream)->tellp() << ")");

				break;
			case FSTREAM:
				out << m_filename << " : read-write file stream ";

				DBG_DO(DBG_UTILITY, out << "(write pos: " << static_cast<fstream*>(m_stream)->tellp()
					<< ", read pos: " << static_cast<fstream*>(m_stream)->tellg() << ")");

				break;
			case SSTREAM:
				out << m_filename << " : string stream. " << endl;

				DBG_DO(DBG_UTILITY, out <<  "(write pos: " << static_cast<stringstream*>(m_stream)->tellp()
					<< ", read pos: " << static_cast<stringstream*>(m_stream)->tellg() << ")");

				break;
		}
		return out.str();
	}

	//////////////////////////////////////////////////////////////
	/// OStream Manager
	//////////////////////////////////////////////////////////////

	OstreamManager::OstreamManager():m_ostreams()
	{
	}

	OstreamManager::~OstreamManager()
	{
		closeAll();
	}

	ostream* OstreamManager::getOstream( const string& name, bool readable,  bool realAppend, bool useString)
	{

		ostreamMapIterator it = m_ostreams.find(name);

		if( it == m_ostreams.end() )			  // not found.
		{

			DBG_DO(DBG_UTILITY, cout << "Create new file " << name << endl);

			return m_ostreams.insert( ostreamMapValue(name,
				StreamElem( name, readable, realAppend, useString))).first->second.stream();
		}
		else									  // already exist
		{

			DBG_DO(DBG_UTILITY, cout << "Find existing ostream " << name << " of info " << it->second.info() << endl);

			// try to see if existing stream type matches what is requrested.
			if( useString && it->second.type() != StreamElem::SSTREAM )
				throw ValueError("file " + name + " has already opened as string file.");
			else if( !useString && it->second.type() == StreamElem::SSTREAM )
				throw ValueError("file " + name + " has already opened as normal file.");
			else if( readable && it->second.type() == StreamElem::OFSTREAM )
				it->second.makeReadable();
			else if( realAppend && ! it->second.append() )
				it->second.makeAppend(true);
			else if( ! realAppend && it->second.append() )
				it->second.makeAppend(false);

			return it->second.stream();
		}

		// will never reach here.
		return NULL;
	}

	bool OstreamManager::hasOstream(const string& filename )
	{
		return (m_ostreams.end() !=  m_ostreams.find(filename) );
	}

	void OstreamManager::listAll()
	{
		for(ostreamMapIterator it=m_ostreams.begin(), itEnd = m_ostreams.end(); it != itEnd;  ++it)
			cout << it->first << " : "  << it->second.info() << endl;
	}

	/// close all files and clean the map
	void OstreamManager::closeAll()
	{
		m_ostreams.clear();
	}

	/// global ostream  manager
	OstreamManager g_ostreams;

	/// return ostream manager
	OstreamManager& ostreamManager()
	{
		return g_ostreams;
	}

	//////////////////////////////////////////////////////////////
	/// Stream provider
	//////////////////////////////////////////////////////////////

	// all flags will be cleared to 0
	StreamProvider::StreamProvider(const string& output, const string& outputExpr)
		:m_filename(output), m_filenameExpr(outputExpr),
		m_flags(0), m_filePtr(NULL)
	{
		if(m_filenameExpr.empty())
			analyzeOutputString(output);
	}

	void StreamProvider::setOutput(const string& output, const string& outputExpr)
	{
		m_filename = output;
		m_filenameExpr = outputExpr;

		if(m_filenameExpr.empty())
			analyzeOutputString(output);
	}

	ostream& StreamProvider::getOstream( PyObject* dict, bool readable)
	{
		DBG_FAILIF( readable && ( ISSETFLAG(m_flags, m_flagNoOutput) || ISSETFLAG(m_flags, m_flagUseDefault)) ,
			SystemError, "A readable file is requested but this Opertor uses cout or cnull.");

		if( ISSETFLAG(m_flags,  m_flagNoOutput ) )
			return cnull();

		// if using cout, return it.
		if( ISSETFLAG(m_flags,  m_flagUseDefault) )
			return cout;

		// if use and close type
		string filename;
		if( !m_filenameExpr.empty())
		{
			DBG_ASSERT( dict != NULL, ValueError,
				"Need to know local dictionary to evaluate filename expression");
			m_filenameExpr.setLocalDict(dict);
			m_filename = m_filenameExpr.valueAsString();

			DBG_DO(DBG_UTILITY, cout << "Filename " << m_filename << endl);

			analyzeOutputString(m_filename);
			if( ISSETFLAG(m_flags, m_flagNoOutput ) )
				return cnull();

			// if using cout, return it.
			if( ISSETFLAG(m_flags, m_flagUseDefault) )
				return cout;
		}
		filename = m_filename;

		if( ISSETFLAG(m_flags, m_flagAppend ) )
		{

			DBG_DO(DBG_UTILITY, cout << "Get a persistent file: "
				<< filename << endl);

			return *ostreamManager().getOstream(filename, readable,
				ISSETFLAG(m_flags ,m_flagRealAppend), ISSETFLAG(m_flags, m_flagUseString));
		}
		else									  // not in append mode, but check if this file is alreay there
		{

			DBG_DO(DBG_UTILITY, cout << "File is not persistent : "
				<< filename << endl);

			if( ! ostreamManager().hasOstream( filename ) )
			{
				if(readable)
					SETFLAG(m_flags, m_flagReadable);
				else
					RESETFLAG(m_flags, m_flagReadable);

				if( readable )
					m_filePtr = new fstream( filename.c_str() );
				else
					m_filePtr = new ofstream( filename.c_str() );

				if( m_filePtr == NULL || !*m_filePtr )
					throw SystemError("Can not create file " + filename );

				return *m_filePtr;
			}
			else
			{

				DBG_DO(DBG_UTILITY, cout << "file " + filename +
					" is already opened as appendable. Use that file instead." << endl);

				RESETFLAG(m_flags, m_flagCloseAfterUse);
				SETFLAG(m_flags, m_flagAppend);
				return *ostreamManager().getOstream(filename, readable, ISSETFLAG(m_flags, m_flagRealAppend),
					ISSETFLAG(m_flags, m_flagUseString));
			}
		}
	}

	/// close ostream and delete ostream pointer.. if it is a ofstream.
	void StreamProvider::closeOstream()
	{
		if( ISSETFLAG(m_flags, m_flagCloseAfterUse ) )
		{
			if( ISSETFLAG(m_flags, m_flagReadable) )
				static_cast<fstream*>(m_filePtr)->close();
			else
				static_cast<ofstream*>(m_filePtr)->close();
			delete m_filePtr;
		}
	}

	void StreamProvider::analyzeOutputString(const string& output)
	{
		if( output.empty() )
		{
			SETFLAG(m_flags, m_flagNoOutput);
			RESETFLAG(m_flags, m_flagCloseAfterUse);
			return;
		}
		else
			RESETFLAG(m_flags, m_flagNoOutput);

		int i;
		for(i=output.size()-1; i>=0; --i)
			if( output[i] == '>' || output[i] == '|' )
				break;

		++i;
		string type;
		string format;
		if( i == 0 )
		{										  // no specification >
			type = ">";
			format=output;
		}
		else
		{
			// >>file i=2, type: 0, size 2, format: start 2, size 6-2
			type = output.substr(0, i);
			format = output.substr(i, output.size() - i);
		}

		RESETFLAG(m_flags, m_flagCloseAfterUse);
		if( type == ">" )
		{
			RESETFLAG(m_flags, m_flagAppend);
			SETFLAG(m_flags, m_flagCloseAfterUse);
		}
		else if( type == ">>" )
		{
			SETFLAG(m_flags, m_flagAppend);
			RESETFLAG(m_flags, m_flagRealAppend);
		}
		else if( type == ">>>")
		{
			SETFLAG(m_flags, m_flagAppend);
			SETFLAG(m_flags, m_flagRealAppend);
		}
		else if( type == "|")
		{
			SETFLAG(m_flags, m_flagAppend);
			SETFLAG(m_flags, m_flagUseString);
		}
		else
			throw ValueError("Ostream types can only be one of '', '>', '>>', '>>>', '|'");

		if(format == "")
		{
			SETFLAG(m_flags, m_flagUseDefault);
			RESETFLAG(m_flags, m_flagCloseAfterUse);
		}
		else
			RESETFLAG(m_flags, m_flagUseDefault);

		DBG_DO(DBG_UTILITY, cout << "Analyzed string is " << output << endl
			<< "Filename is " << format << endl);

		m_filename = format;
	}

	//////////////////////////////////////////////////////////////
	/// Random number generator
	//////////////////////////////////////////////////////////////
	RNG::RNG(const char* rng, unsigned long seed):m_RNG(NULL)
	{
		setRNG(rng, seed);
	}

	RNG::~RNG()
	{
		// free current RNG
		gsl_rng_free(m_RNG);
	}

	unsigned long RNG::generateRandomSeed()
	{
		// now, I need to work hard to get a good seed, considering
		// the use of clusters, and several jobs may be started at the same
		// time
		unsigned long seed;
		FILE *devrandom;
		if ((devrandom = fopen("/dev/urandom", "r")) != NULL)
		{
			fread(&seed, sizeof(seed), 1, devrandom);
			fclose(devrandom);
		}
		else if ((devrandom = fopen("/dev/random", "r")) != NULL)
		{
			fread(&seed, sizeof(seed), 1, devrandom);
			fclose(devrandom);
		}
		else
		{
			// this is not the best method, but I am out of ideas
			// of portable ways to add some other noises
			seed = static_cast<unsigned long>(time(NULL));
		}
		return seed;
	}

	/// choose an random number generator.
	/// This can be done by setting GSL_RNG_TYPE as well.
	void RNG::setRNG(const char * rng, unsigned long seed)
	{
		if (rng != NULL && rng[0] != '\0')
		{
			// locate the RNG
			const gsl_rng_type **t, **t0 = gsl_rng_types_setup ();

			gsl_rng_default = 0;

			/* check GSL_RNG_TYPE against the names of all the generators */

			for (t = t0; *t != 0; t++)
			{
				// require that a RNG can generate full range of integer from 0 to the max of unsigned long int
				if (strcmp (rng, (*t)->name) == 0)
				{
					// free current RNG
					if( m_RNG != NULL )
						gsl_rng_free(m_RNG);

					m_RNG = gsl_rng_alloc( *t );

					DBG_ASSERT(gsl_rng_max(m_RNG) >= MaxRandomNumber && gsl_rng_min(m_RNG) == 0,
						ValueError, "You chosen random number generator can not generate full range of int.");
					break;
				}
			}

			if ( *t == 0 )
				throw SystemError("GSL_RNG_TYPE=" + toStr(rng)
					+ " not recognized or can not generate full range (0-2^32-1) of integers.\n c.f. ListAllRNG()" );
		}
		else
		{
			// free current RNG
			if( m_RNG != NULL )
				gsl_rng_free(m_RNG);

			m_RNG = gsl_rng_alloc( gsl_rng_mt19937 );
		}

		// generate seed
		if( seed == 0 )
			m_seed = generateRandomSeed();
		else
			m_seed = seed;

		// set seed
		gsl_rng_set(m_RNG, m_seed);
	}

	///////////// Weighted sampler //////////////
	/// FIXME: consider adopting R's implementation.
	/// They may be quicker.
	void Weightedsampler::set(const vectorf& weight)
	{
		m_N = weight.size();

		if( m_N == 0)
			return;

		// if(m_fast)                                    // using the fast algorithm
		// {
		// sum of weight
		double w = accumulate(weight.begin(), weight.end(), 0.0);

		DBG_FAILIF( fcmp_le(w, 0), ValueError, "Sum of weight is <= 0.");

		w = m_N/w;

		// initialize p with N*p0,...N*p_k-1
		m_q.resize(m_N);
		for(size_t i=0; i<m_N; ++i)
			m_q[i] = weight[i]*w;
		// initialize Y with values
		m_a.resize(m_N);
		for(size_t i=0; i<m_N; ++i)
			m_a[i] = i;
		// use two sets H and L
		// for efficiency purpose, use a single vector.
		ULONG * HL = new ULONG[m_N];
		ULONG * L = HL;
		ULONG * H = HL+m_N-1;					  // point to the end.

		for(size_t i=0; i<m_N; ++i)
		{
			if( m_q[i] > 1)
				*H-- = i;
			else
				*L++ = i;
		}

		//
		ULONG j,k;
		while( L != HL && H != HL+m_N-1)
		{
			j = *(L-1);
			k = *(H+1);
			m_a[j] = k;
			m_q[k] += m_q[j] - 1;

			L--;								  // remove j from L
			if( m_q[k] < 1. )
			{
				*L++ = k;						  // add k to L
				++H;							  // remove k from H
			}
		}
		delete[] HL;
		//}
		// else                                          // using the bisection method
		// {
		// initialize p with N*p0,...N*p_k-1
		//   m_q.resize(m_N);
		//   for(size_t i=1; i<m_N; ++i)
		//     m_q[i] += m_q[i-1];
		// m_q[m_N-1] is the sum of all weights
		//   for(size_t i=0; i<m_N; ++i)
		//     m_q[i] /= m_q[m_N-1];
		// }
	}

	////////////// Bernulli trials ///////////

	// this is used for BernulliTrials and copyGenotype
	WORDTYPE g_bitMask[WORDBIT];

	BernulliTrials::BernulliTrials(RNG& rng)
		:m_RNG(&rng), m_N(0), m_prob(0), m_table(0), m_pointer(0),
		m_cur(npos)
	{
	}

	BernulliTrials::BernulliTrials(RNG& rng, const vectorf& prob, ULONG trials)
		:m_RNG(&rng), m_N(trials), m_prob(prob), m_table(prob.size()), m_pointer(prob.size()),
		m_cur(npos)
	{
		DBG_FAILIF(trials<=0 , ValueError, "trial number can not be zero.");
		DBG_FAILIF(prob.empty(), ValueError, "probability table can not be empty.");

		// initialize the table
		for(size_t i = 0; i < probSize(); ++i)
		{
			m_table[i].resize(trials);
			m_pointer[i] = BITPTR(m_table[i].begin());
		}
	}

	///
	BernulliTrials::~BernulliTrials()
	{
	}

	void BernulliTrials::setParameter(const vectorf & prob, ULONG trials)
	{
		m_N = trials;
		m_prob = prob;
		m_table.resize(m_prob.size());
		m_pointer.resize(m_prob.size());
		m_cur = npos;							  // will trigger doTrial.

		DBG_FAILIF(trials<=0, ValueError, "trial number can not be zero.");
		DBG_FAILIF(prob.empty(), ValueError, "probability table can not be empty.");

		for(size_t i = 0; i < probSize(); ++i)
		{
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
		if(v)
		{
			// set all to 1
			for(size_t i = 0; i < blk; ++i)
				*ptr++ = ~WORDTYPE(0UL);
			if(rest > 0)
			{
				*ptr |= g_bitMask[rest];
				// upper to 0
				*ptr &= g_bitMask[rest];
			}
		}
		else
		{
			for(size_t i = 0; i < blk; ++i)
				*ptr++ = 0UL;
			if(rest > 0)
				*ptr = 0;						  //~g_bitMask[rest];
		}
	}

#define setBit(ptr, i)    ( *((ptr)+(i)/WORDBIT) |= 1UL << ((i) - ((i)/WORDBIT)*WORDBIT))
#define unsetBit(ptr, i)  ( *((ptr)+(i)/WORDBIT) &= ~ (1UL << ((i) - ((i)/WORDBIT)*WORDBIT)))
	// use a != 0 to avoid compiler warning
#define getBit(ptr, i)    (( *((ptr)+(i)/WORDBIT) & (1UL << ((i) - ((i)/WORDBIT)*WORDBIT))) != 0)

	void BernulliTrials::doTrial()
	{
		DBG_ASSERT(m_N != 0, ValueError, "number of trials should be positive");

		DBG_DO(DBG_UTILITY, cout << "n=" << m_N << " doTrial, cur trial: " << m_cur << endl );

		// for each column
		for(size_t cl = 0, clEnd = probSize(); cl < clEnd; ++cl)
		{
			WORDTYPE * ptr = m_pointer[cl];
			double prob = m_prob[cl];
			if(prob == 0.)
			{
				setAll(cl, false);
			}
			else if( prob == 0.5)				  // random 0,1 bit, this will be quicker
			{
				// set to 0..
				setAll(cl, false);
				// treat a randInt as random bits and set them directly.
				// I.e., we will call 1/16 or 1/32 times of rng for this specifal case.
				// first several blocks
				// WORDTYPE * ptr = BITPTR(succ.begin());
				WORDTYPE tmp;
				size_t blk = m_N / WORDBIT;
				size_t rest = m_N - blk * WORDBIT;
				for(size_t i=0; i < blk; ++i)
				{
					// even if the block size is large (I can not set it to int16_t)
					// I only take the last 16 bit of a rng
					// for the quality of random bits.
					*ptr = 0;
					for(size_t b=0; b < WORDBIT/16; ++b)
					{
						// blocks[i] = static_cast<int16_t>(rng().randGet());
						tmp = rng().randInt(0xFFFF);
						*ptr |= (0xFFFF & tmp) << (b*16);
					}
					ptr ++;
				}
				// last block
				if (rest != 0)
				{
					size_t b = 0;
					for(b=0; b < rest / 16; ++b)
					{
						tmp = rng().randInt(0xFFFF);
						*ptr |= (0xFFFF & tmp) << (b*16);
					}
					// last bits
					rest -= b*16;
					if (rest != 0)
					{
						tmp = rng().randInt(0xFFFF);
						*ptr |= (g_bitMask[rest] & tmp) << b*16;
					}
				}
			}
			// algorithm i Sheldon Ross' book simulation (4ed), page 54
			else if( prob < 0.5)
			{
				// set all to 0, then set some to 1
				setAll(cl, false);
				// it may make sense to limit the use of this method to low p,
				UINT i = 0;
				while( true )
				{
					// i moves at least one.
					i += m_RNG->randGeometric(prob);
					if ( i <= m_N )
						// succ[i-1] = true;
						setBit(ptr, i-1);
					else
						break;
				}
			}
			else if(prob == 1.)
			{
				setAll(cl, true);
			}
			else								  // 1 > m_proc[cl] > 0.5
			{
				// set all to 1, and then unset some.
				setAll(cl, true);
				// it may make sense to limit the use of this method to low p,
				UINT i = 0;
				prob = 1. - prob;
				while( true )
				{
					i += m_RNG->randGeometric(prob);
					if ( i <= m_N )
						// succ[i-1] = false;
						unsetBit(ptr, i-1);
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

	/// get a trial corresponding to m_prob.
	void BernulliTrials::trial()
	{
		if(m_cur == npos || m_cur == m_N - 1)	  // reach the last trial
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
		while(i < sz && !getBit(m_pointer[i], m_cur))
			++i;
		return i >= sz ? npos : i;
	}

	size_t BernulliTrials::probNextSucc(size_t pos) const
	{
		DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
		const size_t sz = probSize();
		if (pos >= (sz-1) || sz == 0)
			return npos;

		++pos;
		while(pos < sz && !getBit(m_pointer[pos], m_cur))
			++pos;
		return pos >= sz ? npos : pos;
	}

	size_t BernulliTrials::trialFirstSucc(size_t idx) const
	{
		size_t blk = m_N / WORDBIT;
		WORDTYPE * ptr = m_pointer[idx];

		size_t i = 0;
		while(i < blk && *ptr++ == 0)
			++i;

		if (i < blk)							  // not at the last blk
		{
			return i * WORDBIT + lowest_bit(*(ptr-1));
		}
		else									  // last block?
		{
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
		WORDTYPE tmp = *ptr & ~ g_bitMask[offset];

		size_t blk =  m_N / WORDBIT;
		if (tmp != 0)
			return i * WORDBIT + lowest_bit(tmp);
		else if (blk == i)
			// if i is the last block, return no.
			return npos;

		// now, go from next block
		++ptr;
		++i;
		while(i < blk && *ptr++ == 0)
			++i;

		if (i < blk)							  // not at the last blk
			return i * WORDBIT + lowest_bit(*(ptr-1));
		else									  // last block?
		{
			size_t rest = m_N - blk * WORDBIT;
			// mask out bits after rest
			size_t tmp = *ptr & g_bitMask[rest];
			if (tmp == 0)
				return npos;
			else
				return blk * WORDBIT + lowest_bit(tmp);
		}
	}

	void BernulliTrials::setTrialSucc(size_t idx, bool succ)
	{
		DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
		if(succ)
			setBit(m_pointer[idx], m_cur);
		else
			unsetBit(m_pointer[idx], m_cur);
	}

	double BernulliTrials::trialSuccRate(UINT index) const
	{
		// efficiency is not considered here
		size_t count = 0;
		for(size_t i=0; i<trialSize(); ++i)
			if(getBit(m_pointer[index], i))
				count++;
		return count/static_cast<double>(m_table[index].size());
	}

	double BernulliTrials::probSuccRate() const
	{
		DBG_ASSERT(m_cur < m_N, ValueError, "Wrong trial index");
		UINT count = 0;
		for(size_t cl = 0, clEnd = probSize(); cl < clEnd; ++cl)
			count += getBit(m_pointer[cl], m_cur) ? 1 : 0;
		return count/static_cast<double>(probSize());
	}

#undef setBit
#undef unsetBit
#undef getBit

	/// random number generator. a global variable.
	/// there might be multiple RNG later.
	RNG g_RNG;

	/// return the global RNG
	RNG& rng()
	{
		return g_RNG;
	}

	/// set global rng
	/// this is temporary since rng() might not exist in the future
	void SetRNG(const string r, unsigned long seed)
	{
		rng().setRNG(r.c_str(), seed);
	}

	void setRNG(const string r, unsigned long seed)
	{
		DBG_WARNING(true, "This function has been renamed to SetRNG() and will be removed at the next major release");
		rng().setRNG(r.c_str(), seed);
	}

	/// list all available RNG.
	vectorstr ListAllRNG()
	{
		vectorstr list;

		const gsl_rng_type **t, **t0;
		gsl_rng * rng;

		t0 = gsl_rng_types_setup();

		for(t=t0; *t !=0; t++)
		{
			rng = gsl_rng_alloc( *t );
			if (gsl_rng_min(rng) == 0 && gsl_rng_max(rng) >= MaxRandomNumber)
				list.push_back((*t)->name);
			gsl_rng_free(rng);
		}

		return list;
	}

	vectorstr listAllRNG()
	{
		DBG_WARNING(true, "This function has been renamed to SetRNG() and will be removed at the next major release");
		return ListAllRNG();
	}

	//////////////////////////////////////////////////////////////
	///  Global debug and initialization related functions
	//////////////////////////////////////////////////////////////

	void gsl_error_handler (const char * reason, const char * ,
		int , int gsl_errno)
	{
		throw SystemError("GSL Error " + toStr( gsl_errno ) + ":\t"
			+ reason);
	}

	/// create a stream buf that print to python sys.stdout
	/// cout will be redirected to this buf to really output
	/// to python console.
	class PythonCoutBuf: public streambuf
	{
		public:
			PythonCoutBuf(){}

		protected:
			int overflow(int c)
			{
				/* In the maste slave mode, slave can also output
				#ifdef SIMUMPI
								// only head node can output to python output
								if(mpiRank() != 0)
									return 0;
				#endif
				*/
				// write out current buffer
				if( pbase() != pptr() )
				{
					// the end of string might not be \0
					char_type * endPtr = pptr();
					char_type endChar = * endPtr;
					*endPtr = '\0';
					//          int len = pptr() - pbase();
					//          char * str = new char[len+1];
					//          strncpy(str, pbase(), len);
					//          str[len] = '\0';
					PySys_WriteStdout("%s",pbase());
					// put original end character back, whatever it is.
					*endPtr = endChar;
					//          delete[] str;
					setp(pbase(), epptr());

				}
				// push character c in
				if( c != EOF )
				{
					// unbuffered, write out this character, do not put into buffer
					if ( pbase() == epptr() )
						PySys_WriteStdout("%c", c);
					else
						sputc(c);
				}
				return 0;
			}
	};

	/// create a null stream buf that discard everything
	class NullStreamBuf: public streambuf
	{
		public:
			NullStreamBuf(){}

		protected:
			int overflow(int)
			{
				return 0;
			}
	};

	/// create an object
	PythonCoutBuf g_pythonCoutBuf;

	/// null stream buf
	NullStreamBuf g_nullStreamBuf;

	/// null stream
	ostream g_cnull( &g_nullStreamBuf );

	/// return null stream
	ostream& cnull()
	{
		return g_cnull;
	}

	/// set default output ("" means standard output)
	void setLogOutput(const string filename)
	{
		static ofstream * outputFile=NULL;
		// close old file is necessary
		if( outputFile != NULL)
		{
			outputFile -> close();
			delete outputFile;
			outputFile = NULL;
		}
		// use stanard output
		if( filename == "" )
			cout.rdbuf( &g_pythonCoutBuf );
		else
		{
			// use a file.
			outputFile = new ofstream(filename.c_str());
			if( outputFile == NULL || !*outputFile )
				throw ValueError("Can not open file " + filename + " to store standard output");
			cout.rdbuf( outputFile->rdbuf());
		}
	}

	/* !COMPILER */

	// record COMPILER, PY_VERSION and __DATE__ , these info will
	// be displayed when simuPOP is loaded.

#ifndef COMPILER
#ifdef __GNUC__
#define COMPILER "[GCC " __VERSION__ "]"
#endif
#endif										  /* !COMPILER */

#ifndef COMPILER
#ifdef __cplusplus
#define COMPILER "[C++]"
#else
#define COMPILER "[C]"
#endif
#endif

#ifndef PLATFORM
#define PLATFORM ""
#endif

	// these macros will be passed from commandline, if not, use the default
#ifndef SIMUPOP_REV
#define REVISION "9999"
#else
	// make passed macro to a real string
#define REVISION MacroQuote(SIMUPOP_REV)
#endif

	int simuRev()
	{
		char * rev = REVISION;
		// can have form xx:xxM etc, or simply a number
		// we certainly like a single number but it is often the case that
		// svn local copy is not up to date.
		int num;
		// first try XX:XX or XX:XXM
		if( sscanf(rev, "%*d:%d", &num) == 1)	  //success
			return num;
		else if( sscanf(rev, "%d", &num) == 1)	  // XX or XXM
			return num;
		else
		{
			cout << "Can not extract revision information from " << REVISION << endl;
			return 0;
		}
	}

	string simuVer()
	{
#ifndef SIMUPOP_VER
		return "snapshot";
#else
		// convert name to a string
		return MacroQuote(SIMUPOP_VER);
#endif
	}

	bool optimized()
	{
#ifdef OPTIMIZED
		return(true);
#else
		return(false);
#endif
	}

	bool mpi()
	{
#ifdef SIMUMPI
		return(true);
#else
		return(false);
#endif
	}

#ifdef SIMUMPI
	// global MPI communicator
	// mpi will be finalized if the communicator is destructed.
	/** this is rediculus */
	int g_mpiArgc = 0;
	char * g_mpiArgv = "";
	char ** g_mpiArgvv = & g_mpiArgv;
	mpi::environment g_mpiEnv(g_mpiArgc, g_mpiArgvv);
	comm g_mpiComm;
	ULONG g_uniqueID = 1;

	comm::~comm()
	{
		if (mpiRank() == 0)
		{
			int action = SLAVE_TERMINATE;
			for(size_t node=1; node<mpiSize(); ++node)
				send(node, 0, action);
		}
	}

	const comm mpiComm()
	{
		return g_mpiComm;
	}

	ULONG uniqueID()
	{
		return g_uniqueID++;
	}
#endif

	UINT mpiRank()
	{
#ifdef SIMUMPI
		return g_mpiComm.rank();
#else
		return(0);
#endif
	}

	UINT mpiSize()
	{
#ifdef SIMUMPI
		return g_mpiComm.size();
#else
		return(0);
#endif
	}

	void mpiBarrier()
	{
#ifdef SIMUMPI
		g_mpiComm.barrier();
#endif
	}

	bool supportXML()
	{
#ifdef __NO_XML_SUPPORT__
		return false;
#else
		return true;
#endif
	}

	string alleleType()
	{
#ifdef LONGALLELE
		return "long";
#else
#ifdef BINARYALLELE
		return "binary";
#else
		return "short";
#endif
#endif
	}

	string compileCompiler()
	{
		return COMPILER;
	}

	string compileDate()
	{
		return __DATE__;
	}

	string compilePyVersion()
	{
		return PY_VERSION;
	}

	string compilePlatForm()
	{
		return PLATFORM;
	}

	// GZIP	\037\213	http://www.ietf.org/rfc/rfc1952.txt
	bool isGzipped(const string & filename)
	{
		// paranoia check
		if (filename.empty())
			return false;

		ifstream ifs(filename.c_str());
		if (!ifs)
			// Couldn't open file...
			return false;

		string str;
		getline(ifs, str);
		return str.substr(0,2) == "\037\213";
	}

	const string fileExtension(const string & filename)
	{
		const string::size_type last_slash = filename.rfind('/');
		string::size_type last_dot;
		if (filename.size() > 3 && filename.substr(filename.size()-3,3) == ".gz")
			last_dot = filename.rfind('.', filename.size()-4);
		else
			last_dot = filename.rfind('.');
		if (last_dot != string::npos &&
			(last_slash == string::npos || last_dot > last_slash))
			return filename.substr(last_dot + 1,
				filename.size() - (last_dot + 1));
		else
			return string();
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
		if ( n < WORDBIT )
		{
			for(size_t i=0; i<n; ++i)
			{
				// set bit according to from bit
				if (*fr_p & (1UL << fr_off))
					*to_p |= (1UL << to_off);
				else
					*to_p &= ~(1UL << to_off);
				// next location
				if (fr_off++ == WORDBIT - 1)
				{
					fr_off = 0;
					++fr_p;
				}
				if (to_off++ == WORDBIT - 1)
				{
					to_off = 0;
					++to_p;
				}
			}
		}
		else if (fr_off == to_off)
		{
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
			for(size_t i=0; i<blks; ++i)
				*++to_p = *++fr_p;
			// the rest of the bits?
			rest -= blks * WORDBIT;
			if(rest != 0)
			{
				// rest = 3
				// from:   xxxxxABC
				// to:     xxxxxCDE
				// mask:   00000111
				//    &:   00000ABC
				//    & :  xxxxx000
				//    |:   xxxxxABC
				to_p ++;
				fr_p ++;
				// mask has top rest zeros.
				mask = g_bitMask[rest];
				*to_p = (*fr_p & mask) | (*to_p & ~mask);
			}
		}
		else if(fr_off < to_off)
		{
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
			maskFrom  = g_bitMask[WORDBIT - shift];
			for(size_t i=0; i < blks; ++i)
			{
				to_p++;
				*to_p = ((*fr_p & ~maskFrom) >> (WORDBIT - shift)) |
					( (*(fr_p+1) & maskFrom) << shift);
				fr_p++;
			}
			// the rest of the bits?
			rest -= blks * WORDBIT;
			if(rest != 0)
			{
				to_p ++;
				if (rest <= shift)
				{
					// rest = 2, shift = 5
					// from:  ABCEDxxx, maskfrom: 00000111
					// to:    xxxxxxAB, maskto:   00000011
					// & ~                        ABCDE000
					// >>                         000ABCDE
					// & naskTo                   000000DE
					// & ~                        xxxxxx00
					// |                          xxxxxxDE
					maskFrom = g_bitMask[WORDBIT - shift];
					maskTo   = g_bitMask[rest];
					*to_p =  (((*(fr_p) & ~maskFrom) >> (WORDBIT - shift)) & maskTo)
						| (*to_p & ~maskTo);
				}
				else
				{
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
					maskFrom1= g_bitMask[rest - shift];
					maskTo   = g_bitMask[rest];

					*to_p =  ((*(fr_p) & ~maskFrom) >> (WORDBIT - shift)) |
						((*(fr_p+1) & maskFrom1) << shift) |
						(*to_p & ~maskTo);
				}
			}
		}
		else									  // fr_off > to_off
		{
			size_t shift = fr_off - to_off;
			WORDTYPE maskFrom = g_bitMask[fr_off];
			WORDTYPE maskFrom1= g_bitMask[shift];
			WORDTYPE maskTo   = g_bitMask[to_off];
			// from:   ABCxxxxx, maskFrom: 00011111
			// from1:  xxxxxxDE, maskFrom1:00000011
			// to:     DEABCxxx, maskTo:   00011111
			//
			*to_p = ((*fr_p & ~maskFrom) >> shift) |
				( (*(fr_p+1) & maskFrom1) << (WORDBIT - shift)) |
				(*to_p & maskTo);
			to_p ++;
			fr_p ++;
			//
			// to other bits
			size_t rest = n - (WORDBIT - to_off);
			size_t blks = rest / WORDBIT;
			//
			// already copied shift bits
			// from:   ABCDEFxx, maskFrom:   00000011
			// from1:  xxxxxxGH, maskFrom:   00000011
			// to:     GHABCDEF
			maskFrom  = g_bitMask[shift];
			for(size_t i=0; i<blks; ++i)
			{
				*to_p = ((*fr_p & ~maskFrom) >> shift) |
					( (*(fr_p+1) & maskFrom) << (WORDBIT - shift));
				fr_p ++;
				to_p ++;
			}
			// the rest of the bits
			rest -= blks * WORDBIT;
			if(rest != 0)
			{
				if (rest < WORDBIT - shift)
				{
					// rest = 2, shift = 3,
					// from:     ABCDExxx, maskFrom: 00000111
					// to:       xxxxxxDE, maskTo:   00000011
					maskFrom = g_bitMask[shift];
					maskTo   = g_bitMask[rest];
					*to_p = (((*fr_p & ~maskFrom) >> shift) & maskTo) |
						(*to_p & ~maskTo) ;
				}
				else
				{
					// rest = 5, shift = 6
					// from:   ABxxxxxx, maskFrom: 00111111
					// from1:  xxxxxCDE, maskFrom1:00000111
					// rest:   xxxCDEAB, maskTo:   00011111
					maskFrom  = g_bitMask[shift];
					maskFrom1 = g_bitMask[rest - (WORDBIT - shift)];
					maskTo    = g_bitMask[rest];
					*to_p = ((*fr_p & ~maskFrom) >> shift) |
						((*(fr_p+1) & maskFrom1) << (WORDBIT - shift)) |
						(*to_p & ~maskTo);
				}
			}
		}
#ifndef OPTIMIZED
		if(debug(DBG_UTILITY))
		{
			if(vectora(fr, fr+n) != vectora(to, to + n))
			{
				cout << "Copy from " << vectora(fr, fr + n)
					<< " to " << vectora(to, to + n) << " failed " << endl;
				cout << "Offsets are " << BITOFF(fr) << " and " << BITOFF(to) << endl;
			}
		}
#endif
	}
#endif

	/** This file is used to initialize simuPOP when being load into
	   python. The swig interface file will has a init% % entry to
	   include this file. */
	bool initialize()
	{
		// tie python stdout to cout
		setLogOutput();

		// for example, if WORDBIT is 8 (most likely 32), we define
		// 0x0, 0x1, 0x3, 0x7, 0xF, 0x1F, 0x3F, 0x7F, 0xFF
		for(size_t i=0; i<WORDBIT; ++i)
		{
			// g_bitMask[i] is the number of 1 count from right.
			for(size_t j=0; j < i; ++j)
				g_bitMask[i] |= (1UL << j);
		}

#ifndef OPTIMIZED
		// turn on some debug info
		TurnOnDebug(DBG_GENERAL);
		// give at most 100 ref count warnings.
#endif
#ifdef Py_REF_DEBUG
		g_refWarningCount = 100;
#endif

		// SIMUPOP_MODULE is passed as name, but we need it to be quoted.
		// Note that under gcc, I could pass the macro from command line
		// using \" \" but this trick does not work under VC.
		// the following process is safer.
#define SimuPOP_Module_Name "##SIMUPOP_MODULE##"

		// set global dictionary/variable
		PyObject* mm = PyImport_AddModule(SimuPOP_Module_Name);
		g_module_vars = SharedVariables(PyModule_GetDict(mm), false);

		// main dictionary
		mm = PyImport_AddModule("__main__");
		g_main_vars = SharedVariables(PyModule_GetDict(mm), false);

		// get population and individual type pointer
		g_swigPopType = SWIG_TypeQuery(PopSWIGType);
		g_swigindividual = SWIG_TypeQuery(IndSWIGType);
		//
		// g_swigOperator = SWIG_TypeQuery(OperatorSWIGType);
		if(g_swigPopType == NULL || g_swigindividual == NULL)
			throw SystemError("Can not get population and individual type pointer, your SWIG version may be run.");

		/// load carray function and type
		initcarray();

		// set gsl error handler
		gsl_set_error_handler(& gsl_error_handler);

#ifndef OPTIMIZED
#ifdef BINARYALLELE
		// binary level genotype copy is compiler dependent and may
		// fail on some systems. Such a test will make sure the binary
		// library work fine.
		testCopyGenotype();
#endif
#endif
		return true;
	}

#ifndef OPTIMIZED
	bool testGappedIterator()
	{
		vectorinfo a(20);
		size_t i;
		for(i=0; i<20; ++i)
			a[i] = i;

		// access to the right elements?
		GappedInfoIterator a5(a.begin(), 5);
		for(i=0; i<4; ++i)
			if(*a5++ != i*5)
				return false;

		// not start from the first?
		GappedInfoIterator a4(a.begin()+1, 4);
		for(i=0; i<4; ++i)
			if(*a4++ != 1+i*4)
				return false;

		// can I get a vector from it?
		vectorinfo b( GappedInfoIterator(a.begin()+1, 4),
			GappedInfoIterator(a.end()+1, 4));
		vectorinfo::iterator it = b.begin();
		for(i=0; it != b.end(); ++it, ++i)
			if(*it != 1+i*4)
				return false;

		return true;
	}

#ifdef BINARYALLELE
	void testCopyGenotype()
	{
		vectora from(1000);
		vectora to(1000);
		for(size_t i = 0; i<100; ++i)
		{
			for(size_t j=0; j<1000; ++j)
			{
				// use != 0 to reduce compiler warning
				from[j] = rng().randInt(2) != 0;
				to[j] = 0;
			}
			size_t from_idx = rng().randInt(300);
			size_t to_idx = rng().randInt(300);
			if (from_idx > to_idx)
				continue;
			size_t length = rng().randInt(500);
			copyGenotype(from.begin() + from_idx,
				to.begin() + to_idx, length);
			if( vectora(from.begin() + from_idx, from.begin() + from_idx + length) !=
				vectora(to.begin() + to_idx, to.begin() + to_idx + length))
			{
				cout << "Copying: " << vectora(from.begin() + from_idx, from.begin() + from_idx + length) << '\n'
					<<  "Obtain:  " << vectora(to.begin() + to_idx, to.begin() + to_idx + length) << '\n'
					<<  "Index From: " << from_idx << " to: " << to_idx << " length: " << length << endl;
				// the error message can not be shown
				throw SystemError("Allele copy test for your system fails.\n"
					"Please email simuPOP mailing list with detailed os and compiler information");
			}
		}
	}
#endif
#endif

}
