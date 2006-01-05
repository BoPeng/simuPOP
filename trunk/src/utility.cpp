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

#include <cstdlib>
#include "time.h"

#include "utility.h"

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
#include "Python.h"
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

// these functions are defined in arraymodule.c which is included
// in simuPOP_wrap.cpp
extern "C" PyObject* newcarrayobjectfrommem(char type, int size, char * ptr, bool copyOver);
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
  /// debug code and DBG_CODE_LENGTH is defined in simupop_cfg.h
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
    "DBG_DEVEL"
  };

  /// set debug area, default to turn all code on
  void TurnOnDebug(DBG_CODE code)
  {
#ifndef OPTIMIZED
    if( code != DBG_ALL )
      g_dbgCode[static_cast<int>(code)] = true;
    else                                          // set all
      g_dbgCode.set();
#endif

  }

  /// turn off debug, default to turn all code off
  void TurnOffDebug(DBG_CODE code)
  {
#ifndef OPTIMIZED
    if( code != DBG_ALL )
      g_dbgCode[static_cast<int>(code)] = false;
    else                                          // reset all
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
    if(obj==NULL)
    {
      val = false;
      return;
    }
    if (obj == Py_True)
      val = true;
    else if (obj == Py_False)
      val = false;
    else
      val = PyInt_AS_LONG(obj) ? true : false;
  }

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

  void PyObj_As_Array(PyObject* obj, vectorf& val)
  {
    if(obj==NULL)
    {
      val = vectorf();
      return;
    }
    if( !PySequence_Check(obj))
      throw ValueError("Expecting a sequence");

    val.resize( PySequence_Size( obj ));

    // assign values
    for(size_t i=0, iEnd=val.size(); i<iEnd; ++i)
      PyObj_As_Double(PySequence_GetItem(obj, i), val[i]);
  }

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
#ifdef LONGALLELE
    return is_carrayobject(obj) &&
      carray_type(obj) == 'B';
#else
    return is_carrayobject(obj) &&
      carray_type(obj) == 'I';
#endif
  }

  PyObject* Int_Vec_As_NumArray(int dim, int* buf, bool copyOver)
  {
    PyObject* res = newcarrayobjectfrommem('i', dim, reinterpret_cast<char*>(buf), copyOver);

    DBG_FAILIF(res==NULL, ValueError, "Can not convert buf to int num array");
    return res;
  }

  PyObject* Double_Vec_As_NumArray(int dim, double* buf, bool copyOver)
  {
    PyObject* res = newcarrayobjectfrommem('d', dim, reinterpret_cast<char*>(buf), copyOver);

    DBG_FAILIF(res==NULL, ValueError, "Can not convert buf to double num array");
    return res;
  }

  PyObject* Allele_Vec_As_NumArray(int dim, Allele* buf, bool copyOver)
  {
#ifdef LONGALLELE
    PyObject* res = newcarrayobjectfrommem('I',dim, reinterpret_cast<char*>(buf), copyOver);
#else
    PyObject* res = newcarrayobjectfrommem('B',dim, reinterpret_cast<char*>(buf), copyOver);
#endif
    DBG_FAILIF(res==NULL, ValueError, "Can not convert buf to Allele num array");
    return res;
  }

  int NumArray_Size(PyObject* obj)
  {
    // return PyArray_Size(obj);
    return carray_length(obj);
  }

  char* NumArray_Data(PyObject* obj)
  {
    // return reinterpret_cast<PyArrayObject*>(obj)->data;
    return carray_data(obj);
  }

  // //////////////////////////////////////////////////////

  // copy constructor
  SharedVariables::SharedVariables(const SharedVariables& rhs)
  {
    m_dict = PyDict_New();
    PyObject *key = 0, *value = 0;

    int i=0;
    while (PyDict_Next(rhs.m_dict, &i, &key, &value))
    {
      Py_INCREF(key);
      Py_INCREF(value);
      PyDict_SetItem(m_dict, key,value);
    }
  }

  SharedVariables::~SharedVariables()
  {
    if( m_ownVars)
    {
      PyDict_Clear(m_dict);
      Py_XDECREF(m_dict);
    }
  }

  /// check the ref count of this shared variable recursively.
  /// report anything that has refcount > 1
  void SharedVariables::checkRefCount()
  {
    checkRefCountRecursive(m_dict);
  }

  void SharedVariables::checkRefCountRecursive(PyObject * obj)
  {
    if(obj->ob_refcnt > 1)
    {
      PyObject* repr = PyObject_Repr(obj);
      cout << obj->ob_refcnt << "\t" << PyString_AsString(repr) << endl;
      Py_DECREF(repr);
    }
    if(PyList_Check(obj))
      for(int i=0;i<PyList_Size(obj);++i)
        checkRefCountRecursive(PyList_GetItem(obj,i));
    if(PyDict_Check(obj))
    {
      PyObject *key, *value;
      int pos = 0;

      while (PyDict_Next(obj, &pos, &key, &value))
        checkRefCountRecursive(value);
    }
  }

  /// setvars C++ ==> Python
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
    PyObject* curParent = m_dict;                 // should be always valid
    PyObject* curKey = PyString_FromString(const_cast<char*>(name.substr(0, i).c_str()));
    int curIdx = 0;
    PyObject * curChild = NULL;

    next:
    if( curType == 1 )
      curChild = PyDict_GetItem(curParent, curKey);
    else
      curChild = PyList_GetItem(curParent, curIdx);
    // get hold ofitem
    // Py_INCREF(curChild);

    if( name[i] == '[' )                          // we need to put an array in ...
    {
      // look for index
      s = ++i;
      for( ; name[i] != ']'; i++)
        if( !isdigit(name[i]) )
          throw ValueError("Expecting numbers: " + name);

      // get index
      size_t idx = atoi( name.substr(s,i-s).c_str());
                                                  // not exist
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
        if( curType == 1 )                        // dictionary?
        {
          PyDict_SetItem(curParent, curKey, curChild);
          Py_XDECREF(curChild);
        }
        else
          PyList_SetItem(curParent, curIdx, curChild);
      }
      else                                        // if exist, if it is a sequence?
      {
        size_t curSize = PyList_Size(curChild);
        /// if the size if enough, get the item
        if( curSize <= idx )
        {
          while(curSize++ < idx+1)
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
    else if( name[i] == '{' )                     // dictionary
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
                                                  // not exist
      if(curChild == NULL || !PyDict_Check(curChild) )
      {
        // Py_XDECREF(curChild);
        // create dictionary with given key, append
        curChild = PyDict_New();
        // keep this
        // Py_INCREF(curChild);

        if( curChild == NULL)
          throw  ValueError("Can not create dictionary" + name);

        // now set that one
        if( curType == 1 )                        // dictionary?
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
    else                                          // end here
    {
      // regardless curChild == NULL
      // set value in
      if( curType == 1 )                          // dictionary?
      {
        PyDict_SetItem(curParent, curKey, const_cast<PyObject*>(val));
        Py_XDECREF( const_cast<PyObject*>(val));
      }
      else
        PyList_SetItem(curParent, curIdx, const_cast<PyObject*>(val));
    }
    //Py_XDECREF(curParent);
    //Py_XDECREF(curChild);
    Py_XDECREF(curKey);
    return const_cast<PyObject*>(val);
  }

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
    PyObject* curParent = m_dict;                 // should be always valid
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

    if( name[i] == '[' )                          // we need to put an array in ...
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
    PyObject* curParent = m_dict;                 // should be always valid
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

    if( name[i] == '[' )                          // we need to put an array in ...
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
    else                                          // now the last piece, we know the parents
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
    Py_INCREF(obj);
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
    Py_INCREF(obj);
    return setVar(name, obj);
  }

  PyObject* SharedVariables::setDoubleNumArrayVar(const string& name, int dim, double * buf, bool copyOver)
  {
    PyObject * obj = Double_Vec_As_NumArray(dim, buf, copyOver);
    // do not know why, no need to inc reference. (checked with valgrind)
    // Py_INCREF(obj);
    return setVar(name, obj);
  }

  PyObject* SharedVariables::setIntNumArrayVar(const string& name, int dim, int * buf, bool copyOver)
  {
    PyObject * obj = Int_Vec_As_NumArray(dim, buf, copyOver);
    // do not know why, no need to inc reference. (checked with valgrind)
    // Py_INCREF(obj);
    return setVar(name, obj);
  }

  int SharedVariables::getVarAsDoubleNumArray(const string& name, double* & buf, bool nameError)
  {
    PyObject * obj = getVar(name, nameError);

    DBG_ASSERT( PyObj_Is_DoubleNumArray(obj), ValueError,
      name + " is not a Python Numeric double array");

    buf = reinterpret_cast<double*>(NumArray_Data(obj));
    return NumArray_Size(obj);
  }

  int SharedVariables::getVarAsIntNumArray(const string& name, int* & buf, bool nameError)
  {
    PyObject * obj = getVar(name, nameError);

    DBG_ASSERT( PyObj_Is_IntNumArray(obj), ValueError,
      name + " is not a Python Numeric int array");

    buf = reinterpret_cast<int*>(NumArray_Data(obj));
    return NumArray_Size(obj);
  }

  void save_none(string& str)
  {
    str += "n";
  }

  PyObject* load_none(const string& str, size_t& offset)
  {
    offset++;
    Py_INCREF(Py_None);
    return Py_None;
  }

  void save_int(string& str, PyObject * args)
  {
    long l = PyInt_AS_LONG((PyIntObject *)args);
    // type + string + ' '
    str += 'i' + toStr(l) + ' ';
  }

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

  void save_long(string& str, PyObject * args)
  {
    long l = PyInt_AS_LONG(args);
    // type +  string + ' '
    str += 'i' + toStr(l) + ' ';
  }

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

  void save_float(string& str, PyObject * args)
  {
    double d = PyFloat_AsDouble(args);
    // type + string
    str += 'f' + toStr(d) + ' ';
  }

  PyObject* load_float(const string& str, size_t& offset)
  {
    int len = 0;
    while( str[offset +len+1] != ' ') len++;
    double d = atof( const_cast<char*>(str.substr(offset+1, len).c_str()));
    offset += len + 2;
    return PyFloat_FromDouble(d);
  }

  void save_string(string& str, PyObject * args)
  {
    char * s = PyString_AsString(args);
    str += 's' + string(s)+'\0';
  }

  PyObject* load_string(const string& str, size_t& offset)
  {
    char * s = const_cast<char*>(str.c_str()) + offset + 1;
    size_t len = strlen(s);
    offset += len + 2;
    DBG_ASSERT( str[offset-1] == '\0', SystemError,
      "Error while loading string from str");

    return PyString_FromString(s);
  }

  void saveObj(string& str, PyObject* args);
  PyObject* loadObj(const string& vars, size_t& offset);

  void save_dict(string& str, PyObject* args)
  {
    PyObject *key = 0, *value = 0;

    str += 'd';                                   // dictionary
    int i = 0;
    while (PyDict_Next(args, &i, &key, &value))
    {
      saveObj(str, key);
      saveObj(str, value);
    }
    str += 'e';                                   // ending of a dictionary
  }

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
    offset++;                                     // skip 'e'
    return d;
  }

  void save_list(string& str, PyObject* args)
  {
    str += 'L';                                   // dictionary
    int len = PyList_Size(args);
    for (int i = 0; i < len; i++)
    {
      PyObject* elem = PyList_GET_ITEM((PyListObject *)args, i);
      saveObj(str, elem);
    }

    str += 'e';                                   // ending of a dictionary
  }

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
    offset++;                                     // skip 'e'
    return d;
  }

  void save_carray(string& str, PyObject* args)
  {
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
    PyObject* arr = newcarrayobjectfrommem(type, len, ptr, true);
    offset += len*carray_itemsize(arr) + 1;
    return arr;
  }

  void saveObj(string& str, PyObject* args)
  {
    PyTypeObject *type;

    if (args == Py_None)
    {
      save_none(str);
      return;
    }

    type = args->ob_type;

    switch (type->tp_name[0])
    {
      case 'i':
        if (type == &PyInt_Type)
        {
          save_int(str, args);
          return;
        }
        break;
      case 'd':
        if (type == &PyDict_Type)
        {
          save_dict(str, args);
          return;
        }
      case 's':
        if (type == &PyString_Type)
        {
          save_string(str, args);
          return;
        }
        break;
      case 'l':
        if (type == &PyLong_Type)
        {
          save_long(str, args);
          return;
        }
        else if (type == &PyList_Type)
        {
          save_list(str, args);
          return;
        }
        break;
      case 'a':
        if (type == &Arraytype)
        {
          save_carray(str, args);
          return;
        }
        break;
      case 'f':
        if (type == &PyFloat_Type)
        {
          save_float(str, args);
          return;
        }
        break;
    }
    // some other unknown type
    DBG_DO(DBG_UTILITY, cout << "Warning: object of type '" + toStr(type->tp_name) + "' cannot be saved. Use none.");
    save_none(str);
  }

  PyObject* loadObj(const string& vars, size_t& offset)
  {
    switch( vars[offset])
    {
      case 'd':                                   //
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
      case 'a':
        return load_carray(vars, offset);
      case 'n':
        return load_none(vars, offset);
      default:
      {
        DBG_DO(DBG_UTILITY, cout << endl);
        throw ValueError("Unknown type code at offset " + toStr(offset));
      }
    }
  }

  string SharedVariables::asString() const
  {
    // go through each variable and save
    string str;
    saveObj(str, m_dict);
    str += 'e';                                   // ending character
    return str;
  }

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
  swig_type_info* g_swigPopType, *g_swigIndType;

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
    return SWIG_NewPointerObj(p, g_swigIndType, 0);
  }

  // ////////////////////////////////////////////////////////////
  // Expression evaluation
  // ////////////////////////////////////////////////////////////
  // because of ref count, need to define copier

  Expression::Expression(const Expression& rhs)
    : m_expr(rhs.m_expr), m_stmts(rhs.m_stmts), m_locals(rhs.m_locals)
  {
    if(m_expr != NULL)
      Py_INCREF(m_expr);
    if(m_stmts != NULL)
      Py_INCREF(m_stmts);
  }

  Expression::~Expression()
  {
    // release compiled code
    if( m_expr != NULL )
      Py_DECREF(m_expr);
    if( m_stmts != NULL )
      Py_DECREF(m_stmts);
  }

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
      throw ValueError("Statement '" + stmts + "' is not valid.");
  }

  /// python expression
  PyObject* Expression::evaluate()
  {
    if( m_expr == NULL && m_stmts == NULL )
      return NULL;

    DBG_ASSERT( mainVars().dict() != NULL && m_locals != NULL,
      ValueError, "Can not evalulate. Dictionary is empty!");

    if( m_stmts != NULL)
    {
      if( PyEval_EvalCode((PyCodeObject*) m_stmts,
        mainVars().dict(), m_locals ) == NULL)
      {
        PyErr_Print();
        throw SystemError("Evalulation of statements failed");
      }
    }

    PyObject * res=NULL;
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

  /// returned values.
#define ExpressionValueAsType(T, TypeName, init) \
T Expression::valueAs##TypeName() \
{ \
  PyObject* res = evaluate(); \
  if( res==NULL ) \
  return init; \
  T val; \
  PyObj_As_##TypeName(res, val); \
  Py_XDECREF(res); \
  return val; \
}

    ExpressionValueAsType(bool, Bool, false);
  ExpressionValueAsType(int, Int, 0);
  ExpressionValueAsType(double, Double, 0.0);
  ExpressionValueAsType(string, String, "");
  ExpressionValueAsType(vectorf, Array, vectorf());
  ExpressionValueAsType(strDict, StrDict, strDict());
  ExpressionValueAsType(intDict, IntDict, intDict());

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
    else                                          // no string, usual file
    {
      m_append = realAppend;
      if( readable )
      {
        m_type = FSTREAM;
        if (realAppend )                          // readable , append
          m_stream = new fstream(name.c_str(),  std::ios::in | std::ios::out | std::ios::ate );
        else                                      // readable, ! append
          // existing file will be truncated...
          m_stream = new fstream(name.c_str(),  std::ios::in | std::ios::trunc | std::ios::out);
      }
      else
      {
        m_type = OFSTREAM;
        if( realAppend)                           // ! readable, append
          m_stream = new ofstream(name.c_str(), std::ios::out | std::ios::app );
        else                                      //  !readable, !append )
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
    closeAll(true);
  }

  ostream* OstreamManager::getOstream( const string& name, bool readable,  bool realAppend, bool useString)
  {

    ostreamMapIterator it = m_ostreams.find(name);

    if( it == m_ostreams.end() )                  // not found.
    {

      DBG_DO(DBG_UTILITY, cout << "Create new file " << name << endl);

      return m_ostreams.insert( ostreamMapValue(name,
        StreamElem( name, readable, realAppend, useString))).first->second.stream();
    }
    else                                          // already exist
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
  void OstreamManager::closeAll(bool closeAppend)
  {
    ostreamMap temp;

    // pick out those really persistent files
    for(ostreamMapIterator it = m_ostreams.begin(), itEnd = m_ostreams.end(); it != itEnd;  ++it)
    {
      if( !closeAppend  && it->second.append() )
      {

        DBG_DO(DBG_UTILITY, cout << "Do not close file: " << it->first << endl);

        temp.insert( *it );
        continue;
      }
    }
    m_ostreams.clear();
    m_ostreams.swap(temp);

    // the destructor of original ostreams will be called but
    // the original pointers have been set to NULL during the copy constructor to temp.
    // this is **very** tricky and error prone. I wish I can make it better later.

    // even if we do not close realAppend file, we flush them to make sure we can see them on disk
    for(ostreamMapIterator it = m_ostreams.begin(), itEnd = m_ostreams.end(); it != itEnd;  ++it)
      it->second.stream()->flush();

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
    else                                          // not in append mode, but check if this file is alreay there
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
    {                                             // no specification >
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
        if (strcmp (rng, (*t)->name) == 0)
        {
          // free current RNG
          if( m_RNG != NULL )
            gsl_rng_free(m_RNG);

          m_RNG = gsl_rng_alloc( *t );
          break;
        }
      }

      if ( *t == 0 )
        throw SystemError("GSL_RNG_TYPE=" + toStr(rng)
          + " not recognized\n c.f. listAllRNG()" );
    }
    else
    {
      // free current RNG
      if( m_RNG != NULL )
        gsl_rng_free(m_RNG);

      m_RNG = gsl_rng_alloc( gsl_rng_mt19937 );
    }

    // set seed
    if( seed == 0 )
      gsl_rng_set(m_RNG, static_cast<UINT>(time(0)));

  }

  ///////////// Weighted Sampler //////////////
  /// FIXME: consider adopting R's implementation.
  /// They may be quicker.
  void WeightedSampler::set(const vectorf& weight)
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
    ULONG * H = HL+m_N-1;                         // point to the end.

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

      L--;                                        // remove j from L
      if( m_q[k] < 1. )
      {
        *L++ = k;                                 // add k to L
        ++H;                                      // remove k from H
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

  BernulliTrials::BernulliTrials(RNG& rng)
    :m_RNG(&rng), m_N(0), m_prob(0), m_table(0),
    m_cur(0), m_bitSet(0)
  {
  }

  BernulliTrials::BernulliTrials(RNG& rng, const vectorf& prob, ULONG trials)
    :m_RNG(&rng), m_N(trials), m_prob(prob), m_table(prob.size()),
    m_cur(trials), m_bitSet(prob.size())
  {
    DBG_FAILIF( trials<=0 , ValueError, "trial number can not be zero.");
    DBG_FAILIF( prob.empty(), ValueError, "probability table can not be empty.");

    for(vector<BitSet>::iterator it=m_table.begin(), itEnd = m_table.end();
      it != itEnd;  ++it)
    {
      it->resize(trials);
    }
  }

  ///
  BernulliTrials::~BernulliTrials()
  {
  }

  void BernulliTrials::setParameter(const vectorf& prob, ULONG trials)
  {
    m_N = trials;
    m_prob = prob;
    m_table.resize(m_prob.size());
    m_bitSet.resize(m_prob.size());
    m_cur = m_N;                                  // will trigger doTrial.

    DBG_FAILIF( trials<=0 , ValueError, "trial number can not be zero.");
    DBG_FAILIF( prob.empty(), ValueError, "probability table can not be empty.");

    for(vector<BitSet>::iterator it=m_table.begin(), itEnd = m_table.end();
      it != itEnd;  ++it)
    {
      it->resize(trials);
    }
  }

  void BernulliTrials::doTrial()
  {
    DBG_FAILIF( m_N == 0, ValueError, "number of trials should be positive");

    DBG_DO(DBG_UTILITY, cout << "n=" << m_N << " doTrial, cur trial: " << m_cur << endl );

    // for each column
    for(size_t cl = 0, clEnd = m_prob.size(); cl < clEnd; ++cl)
    {
      // clear previous result
      BitSet& succ = m_table[cl];
      succ.reset();

      if( m_prob[cl] == 0.5)                      // random 0,1 bit
      {
        // treat a randInt as random bits and set them directly.
        // I.e., we will call 1/16 or 1/32 times of rng for this specifal case.
        unsigned long max = rng().max();
        int block;
        unsigned long ri;
        if( max == 0xFFFFFFFF)
        {
          block = 32;
          ri = rng().randGet();
        }
        else
        {
          block = 16;
          ri = rng().randInt(0xFFFF);
        }
        // loop
        int bit = 0;
        for(UINT i=0; i < m_N; ++i)
        {
          if(bit == block)                        // generate a block
          {
            if(block == 32)
              ri = rng().randGet();
            else
              ri = rng().randInt(0xFFFF);
            bit = 0;
          }
          // assign many bits a time.
          if( ((ri >> bit) & 0x1) == 1)
            succ.set(i);
          bit++;
        }
      }
      else if( m_N > 100)
        // it may make sense to limit the use of this method to low p,
        // but it rurns out that this may be a good idea.
      {
        // number of success trials
        // M should be 0...m_N
        ULONG  M = m_RNG->randBinomial( m_N, m_prob[cl] );

        ULONG loc;

        // randomly choose events
        for(UINT i=0; i<M; ++i)
        {
          // choose among [0,m_N)
          do
          {
            loc = m_RNG->randInt(m_N);
          }while( succ[loc] );
          succ.set(loc);
        }
      }
      else
      {
        double p = m_prob[cl];
        for(UINT i = 0; i < m_N; ++i)
          if( m_RNG->randUniform01() < p )
            succ.set(i);
      }
    }
    m_cur = 0;
    DBG_DO(DBG_UTILITY, cout << " doTrial done" << m_cur << endl );
  }

  UINT BernulliTrials::curTrial()
  {
    return m_cur;
  }

  const BitSet& BernulliTrials::trial()
  {
    if(m_cur == m_N )                             // reach the last trial
      doTrial();

    for(size_t cl = 0, clEnd = m_prob.size(); cl < clEnd; ++cl)
      m_bitSet.set( cl, m_table[cl][m_cur] );

    m_cur++;
    return m_bitSet;
  }

  const BitSet& BernulliTrials::succ(UINT index)
  {
    // since m_cur will not change when using succ
    // doTrial should be called by caller.
    // if(m_cur == m_N )                             // reach the last trial
    //  doTrial();

    DBG_FAILIF( index >= m_table.size(), ValueError, "succ: index out of range");

    return m_table[index];
  }

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
  void setRNG(const string r, unsigned long seed)
  {
    rng().setRNG(r.c_str(), seed);
  }

  /// list all available RNG.
  vectorstr listAllRNG()
  {
    vectorstr list;

    const gsl_rng_type **t, **t0;

    t0 = gsl_rng_types_setup();

    for(t=t0; *t !=0; t++)
      list.push_back((*t)->name);

    DBG_DO(DBG_GENERAL, cout << "Please refer to GSL (GNU Scientific Library) website:" << endl
      << "  http://www.gnu.org/software/gsl/manual/gsl-ref_17.html#SEC271 " << endl
      << "about the details of each algorithm. " << endl
      << "c.f. setRNG()"<<endl);

    return list;
  }

  /// show turned on bits
  vectorlu bitSet(const BitSet& bs)
  {
    vectorlu list;
    for(BitSet::size_type i=0, iEnd = bs.size(); i < iEnd;  i++)
      if( bs[i] )
        list.push_back(i);
    return list;
  }

  //////////////////////////////////////////////////////////////
  ///  Global debug and initialization related functions
  //////////////////////////////////////////////////////////////

  // new handler (not used by now.)
  //void my_new_handler()
  //{
  //  throw OutOfMemory("Not enough memory available. Please free some memory or simulate smaller population.");
  //}

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
  /** This file is used to initialize simuPOP when being load into
     python. The swig interface file will has a init% % entry to
     include this file. */
  bool initialize()
  {
    // tie python stdout to cout
    setLogOutput();

#ifndef OPTIMIZED
    // turn on some debug info
    TurnOnDebug(DBG_GENERAL);
#endif

    // SIMUPOP_MODULE is passed as name, but we need it to be quoted.
    // Note that under gcc, I could pass the macro from command line
    // using \" \" but this trick does not work under VC.
    // the following process is safer.
#define SimuPOP_Module_Name " ## SIMUPOP_MODULE ## "

    // set global dictionary/variable
    PyObject* mm = PyImport_AddModule(SimuPOP_Module_Name);
    g_module_vars = SharedVariables(PyModule_GetDict(mm), false);

    // main dictionary
    mm = PyImport_AddModule("__main__");
    g_main_vars = SharedVariables(PyModule_GetDict(mm), false);

    // get population and individual type pointer
    g_swigPopType = SWIG_TypeQuery(PopSWIGType);
    g_swigIndType = SWIG_TypeQuery(IndSWIGType);
    if( g_swigPopType == NULL || g_swigIndType == NULL)
      throw SystemError("Can not get population and individual type pointer, your SWIG version may be run.");

    /// load carray function and type
    initcarray();

    // set gsl error handler
    gsl_set_error_handler(& gsl_error_handler);

    return true;
  }
  /* !COMPILER */

  // record COMPILER, PY_VERSION and __DATE__ , these info will
  // be displayed when simuPOP is loaded.

#ifndef COMPILER
#ifdef __GNUC__
#define COMPILER "[GCC " __VERSION__ "]"
#endif
#endif                                          /* !COMPILER */

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
#ifndef SIMUPOP_VER
#define SIMUPOP_VER "snapshot"
#endif

#ifndef SIMUPOP_REV
#define SIMUPOP_REV "9999"
#endif

  int simuRev()
  {
    char * rev = SIMUPOP_REV;
    // can have form xx:xxM etc, or simply a number
    // we certainly like a single number but it is often the case that
    // svn local copy is not up to date.
    int num;
    // first try XX:XX or XX:XXM
    if( sscanf(rev, "%*d:%d", &num) == 1)         //success
      return num;
    else if( sscanf(rev, "%d", &num) == 1)        // XX or XXM
      return num;
    else
    {
      cout << "Can not extract revision information from " << SIMUPOP_REV << endl;
      return 0;
    }
  }

  string simuVer()
  {
    return SIMUPOP_VER;
  }

  bool optimized()
  {
#ifdef OPTIMIZED
    return(true);
#else
    return(false);
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

}
