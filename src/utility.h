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
* GappedAlleleIterator
* StreamElem
* StreamManager
* StreamProvider
* RNG
* BernulliTrials
*/

#include "simupop_cfg.h"

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

/// for random number generator
#include "gsl/gsl_sys.h"                          // for floating point comparison
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

/// for bernulli trials.
#include "boost/dynamic_bitset.hpp"
typedef boost::dynamic_bitset<> BitSet;

#include "Python.h"

namespace simuPOP
{

  // ////////////////////////////////////////////////////////////
  // / Debug and info functions
  // ////////////////////////////////////////////////////////////

  /// set debug code, default to turn all code on
  void TurnOnDebug(DBG_CODE code=DBG_ALL);

  /// turn off debug, default to turn all code off
  void TurnOffDebug(DBG_CODE code=DBG_ALL);

#ifndef OPTIMIZED
  /// test if one code is turned on
  bool debug(DBG_CODE code);
#endif

  /// show all dbg codes (print to cout)
  void ListDebugCode();

  /// dbg string for a code
  string dbgString(DBG_CODE code);

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

#define fcmp_lt(a,b) (gsl_fcmp(a,b,cmp_epsilon) == -1)
#define fcmp_le(a,b) (gsl_fcmp(a,b,cmp_epsilon) <= 0 )
#define fcmp_gt(a,b) (gsl_fcmp(a,b,cmp_epsilon) == 1 )
#define fcmp_ge(a,b) (gsl_fcmp(a,b,cmp_epsilon) >= 0 )
#define fcmp_eq(a,b) (gsl_fcmp(a,b,cmp_epsilon) == 0 )
#define fcmp_ne(a,b) (gsl_fcmp(a,b,cmp_epsilon) != 0 )
#define f_asBool(a)  (fcmp_ne(a, 0.))
}


namespace std
{
  /// CPPONLY how to output any std::pair
  template<class T1, class T2>
    ostream& operator<< (ostream& out, const pair<T1, T2> & pair)
  {
    out << "(" << setw(3) << pair.first << "," << setw(3) << pair.second << ")";
    return out;
  }

  /// CPPONLY how to output any vector.
  template<class T>
    ostream& operator<< (ostream& out, const vector<T> & vec)
  {
    if( ! vec.empty() )
    {
      typename vector<T>::const_iterator it = vec.begin();
      out <<  *it;
      for( ++it; it!= vec.end(); ++it)
        out << ", " << *it ;
    }
    return out;
  }

  /// CPPONLY how to output a dictionary
  ostream& operator<<(ostream& out, const strDict& dict);

  /// CPPONLY how to output a dictionary
  ostream& operator<<(ostream& out, const intDict& dict);

}


namespace simuPOP
{

  /// CPPONLY operator to tell the affected status of an individual
  template<class IndType>
    class isAffected
  {
    public:
      isAffected(){};

      bool operator() ( const IndType& ind)
      {
        return ind.affected();
      }
  };

  /** \brief CPPONLY iterator to access an allele across all ploidy and individuals.

  The reason why this iterator is called GappedAlleleIterator is that
  when alleles at the same locus are placed in equal distance ( a strategy
  currently employeed by Population ), this iterator can
  iterate through all alleles quickly by moving internal pointer
  \e m_step steps a time.

  e.g. the following code set alleles at the first locus all to 2.
  \code
  for(typename Pop::AlleleIterator it=pop.alleleBegin(0),
  itEnd = pop.alleleEnd(0); it < itEnd;  ++it)
  *it = 2;
  \endcode

  This iterator is designed to be a random access operator compatible to
  that of STL. All operators should work as expected so no comments are given.

  \note note that it+1 is not the next locus, but the next allele on the next
  chromosome.
  */
  class GappedAlleleIterator
  {
    public:

      typedef std::random_access_iterator_tag iterator_category;
      typedef Allele                     value_type;
      typedef int                        difference_type;
      typedef value_type&                reference;
      typedef value_type*                pointer;

      /// CPPONLY
      GappedAlleleIterator():m_step(1), m_ptr(NULL)
      {
      }

      /// CPPONLY
      GappedAlleleIterator( pointer p, difference_type s = 1)
        :m_step(s), m_ptr(p)
      {
      }

      /// CPPONLY
      GappedAlleleIterator(const GappedAlleleIterator& rhs)
      {
        m_ptr = rhs.m_ptr;
        m_step = rhs.m_step;
      }

      /// CPPONLY
      GappedAlleleIterator& operator= (const GappedAlleleIterator& rhs)
      {
        m_ptr = rhs.m_ptr;
        m_step = rhs.m_step;
        return *this;
      }

      /// CPPONLY
      ~GappedAlleleIterator()
      {
      }

      /// get pointer so that ptr()+1 will be the next locus.
      pointer ptr()
      {
        return m_ptr;
      }

      /// CPPONLY
      difference_type step()
      {
        return m_step;
      }

      /// CPPONLY
      value_type& operator*() const
      {
        return *m_ptr;
      }

      /// CPPONLY
      value_type& operator[] (difference_type diff) const
      {
        return *(m_ptr + diff*m_step);
      }

      /// CPPONLY
      GappedAlleleIterator operator++ (difference_type)
      {
        GappedAlleleIterator tmp = *this;
        m_ptr += m_step;
        return tmp;
      }

      /// CPPONLY
      GappedAlleleIterator & operator++ ()
      {
        m_ptr += m_step;
        return *this;
      }

      /// CPPONLY
      GappedAlleleIterator operator-- (difference_type)
      {
        GappedAlleleIterator tmp = *this;

        m_ptr -= m_step;
        return tmp;
      }

      /// CPPONLY
      GappedAlleleIterator & operator-- ()
      {
        m_ptr -= m_step;
        return *this;
      }

      /// CPPONLY
      GappedAlleleIterator & operator+= (difference_type diff)
      {
        m_ptr += diff * m_step;
        return *this;
      }

      /// CPPONLY
      GappedAlleleIterator & operator-=(difference_type diff)
      {
        m_ptr -= diff * m_step;
        return *this;
      }

      /// CPPONLY
      GappedAlleleIterator operator+ (difference_type diff)
      {
        GappedAlleleIterator tmp = *this;

        tmp += diff;
        return tmp;

      }

      /// CPPONLY
      GappedAlleleIterator operator- (difference_type diff)
      {
        GappedAlleleIterator tmp = *this;

        tmp -= diff;
        return tmp;
      }

      /// CPPONLY
      difference_type operator- (GappedAlleleIterator rhs)
      {
        return (m_ptr - rhs.m_ptr)/m_step;
      }

      /// CPPONLY
      bool operator== ( const GappedAlleleIterator & rhs)
      {
        return (m_ptr == rhs.m_ptr);
      }

      /// CPPONLY
      bool operator!= ( const GappedAlleleIterator&  rhs)
      {
        return (m_ptr != rhs.m_ptr);
      }

      /// CPPONLY
      bool operator>= ( const GappedAlleleIterator & rhs)
      {
        return (m_ptr >= rhs.m_ptr);
      }

      /// CPPONLY
      bool operator<= ( const GappedAlleleIterator & rhs)
      {
        return (m_ptr <= rhs.m_ptr);
      }

      /// CPPONLY
      bool operator> ( const GappedAlleleIterator & rhs)
      {
        return (m_ptr > rhs.m_ptr);
      }

      /// CPPONLY
      bool operator< ( const GappedAlleleIterator & rhs)
      {
        return (m_ptr < rhs.m_ptr);
      }

    private:
      /// gap between alleles. This is usally the total number of loci on all chromosoms.
      difference_type m_step;

      /// internal pointer
      pointer m_ptr;

  };

  // ////////////////////////////////////////////////////////////
  // / Shared variables
  // ////////////////////////////////////////////////////////////

  /// CPPONLY
  void PyObj_As_Bool(PyObject *obj, bool& val);

  /// CPPONLY
  void PyObj_As_Int(PyObject *obj, int& val);

  /// CPPONLY
  void PyObj_As_Double(PyObject *obj, double& val);

  /// CPPONLY
  void PyObj_As_String(PyObject *obj, string& val);

  /// CPPONLY
  void PyObj_As_Array(PyObject *obj, vectorf& val);

  /// CPPONLY
  void PyObj_As_StrDict(PyObject *obj, strDict& val);

  /// CPPONLY
  void PyObj_As_IntDict(PyObject *obj, intDict& val);

  /// CPPONLY
  bool PyObj_Is_IntNumArray(PyObject * obj);

  /// CPPONLY
  bool PyObj_Is_DoubleNumArray(PyObject * obj);

  /// CPPONLY
  bool PyObj_Is_AlleleNumArray(PyObject * obj);

  /// CPPONLY
  PyObject* Int_Vec_As_NumArray(int dim, int* buf, bool copyOver=true);

  /// CPPONLY
  PyObject* Double_Vec_As_NumArray(int dim, double* buf, bool copyOver=true);

  /// CPPONLY
  PyObject* Allele_Vec_As_NumArray(int dim, Allele* buf, bool copyOver=true);

  /// CPPONLY
  int NumArray_Size(PyObject* obj);

  /// CPPONLY
  char* NumArray_Data(PyObject* obj);

  // ///////////////////////////////////////////////////////
  /** shared variables.
   CPPONLY

  This class set and read Python variables using the given
  dictionary.

  */

  class SharedVariables
  {
    public:

      /// CPPONLY
      SharedVariables():m_dict(NULL),m_ownVars(false)
      {
      }

      /// CPPONLY
      SharedVariables(PyObject * dict, bool ownVars)
        :m_dict(dict),m_ownVars(ownVars)
      {
        // if not given a pre-existing dictionary, create one
        if(m_dict == NULL)
          m_dict = PyDict_New();

        DBG_ASSERT( m_dict != NULL, SystemError, "Can not create a new dictionary");

        if(! PyDict_Check( m_dict) )
          throw SystemError("Invaid dictionary. The local namespace may have been cleared.");
      }

      /// CPPONLY, copy and increase ref count
      SharedVariables(const SharedVariables& rhs);

      /// CPPONLY swap two ShatedVariables
      void swap(SharedVariables& rhs)
      {
        std::swap(m_dict, rhs.m_dict);
        std::swap(m_ownVars, rhs.m_ownVars);
      }

      /// CPPONLY  destructor
      /// I can not clear dict here since
      /// a resize of g_vars will copy this object and
      /// hence call this destructore.
      ~SharedVariables();

      void clear()
      {
        PyDict_Clear(m_dict);
      }

      void setDict(PyObject* dict)
      {
        m_dict = dict;
      }

      void checkRefCount();

      /// CPPONLY set arbitrary value (type)
      PyObject* setVar(const string& name, const PyObject* val);

      /// CPPONLY get variable, NULL is returned if nothing is found.
      PyObject* getVar(const string& name, bool nameError=true);

      bool hasVar(const string& name)
      {
        // not null = has variable
        return PyDict_GetItemString(m_dict, const_cast<char*>(name.c_str())) != NULL;
      }

      /// remove variable
      void removeVar(const string& name);

      ///CPPONLY
      PyObject* setBoolVar(const string& name, const bool val);

      ///CPPONLY
      PyObject* setIntVar(const string& name, const int val);

      ///CPPONLY
      PyObject* setDoubleVar(const string& name, const double val);

      ///CPPONLY
      PyObject* setStringVar(const string& name, const string& val);

      ///CPPONLY
      PyObject* setStrDictVar(const string& name, const strDict& val);

      ///CPPONLY
      PyObject* setIntDictVar(const string& name, const intDict& val);

      ///CPPONLY
      PyObject* setDoubleNumArrayVar(const string& name, int dim, double* buf, bool copyOver=true);

      ///CPPONLY
      PyObject* setIntNumArrayVar(const string& name, int dim, int* buf, bool copyOver=true);

      /// CPPONLY
      bool getVarAsBool(const string& name, bool nameError=true)
      {
        bool val;
        PyObj_As_Bool( getVar(name, nameError), val);
        return val;
      }

      /// CPPONLY
      int getVarAsInt(const string& name, bool nameError=true)
      {
        int val;
        PyObj_As_Int( getVar(name, nameError), val);
        return val;
      }

      /// CPPONLY
      double getVarAsDouble(const string& name, bool nameError=true)
      {
        double val;
        PyObj_As_Double( getVar(name, nameError), val);
        return val;
      }

      /// CPPONLY, obsolete
      // void numArrayOwnMem(const string& name);

      /// CPPONLY
      string getVarAsString(const string& name, bool nameError=true)
      {
        string val;
        PyObj_As_String( getVar(name, nameError), val);
        return val;
      }

      /// CPPONLY
      strDict getVarAsStrDict(const string& name, bool nameError=true)
      {
        strDict val;
        PyObj_As_StrDict( getVar(name, nameError), val);
        return val;
      }

      /// CPPONLY
      intDict getVarAsIntDict(const string& name, bool nameError=true)
      {
        intDict val;
        PyObj_As_IntDict( getVar(name, nameError), val);
        return val;
      }

      /// CPPONLY
      /// put pointer in bug
      /// return size
      int getVarAsDoubleNumArray(const string& name, double* &buf, bool nameError=true);

      /// CPPONLY
      int getVarAsIntNumArray(const string& name, int* &buf, bool nameError=true);

      PyObject*& dict()
      {
        return m_dict;
      }

      /// CPPONLY, save m_dist as string
      string asString() const;

      void fromString(const string& vars);

    private:

      ///
      void checkRefCountRecursive(PyObject* obj);

      /// the list
      PyObject * m_dict;

      /// whether or not the object owns the dictionary
      bool m_ownVars;

  };

  /// CPPONLY
  /// get main dictionary (user namespace)
  SharedVariables& mainVars();

  /// CPPONLY
  /// get module dictionary (it is different than mainDict!
  SharedVariables& moduleVars();

  // ////////////////////////////////////////////////////////////
  // Expression evaluation
  // ////////////////////////////////////////////////////////////

  /**  \brief expression CPPONLY

  This class evaluate python expressions.

  */
  class Expression
  {
    public:
      /// CPPONLY
      Expression(const string& expr="", const string& stmts="",
        PyObject* locals=NULL)
        : m_expr(NULL), m_stmts(NULL), m_locals(locals)
      {
        if(m_locals == NULL)
          m_locals = mainVars().dict();

        // empty expression
        if( expr.empty() && stmts.empty())
          return;

        compileExpr(expr);
        compileStmts(stmts);
      }

      /// CPPONLY
      ~Expression();

      /// Copy constructor, need to be defined because of ref count issue.
      Expression(const Expression& rhs);

      /// CPPONLY
      /// set local dictionary
      void setLocalDict(PyObject* dict)
      {
        m_locals = dict;
      }

      /// CPPONLY
      bool empty()
      {
        return m_expr==NULL && m_stmts==NULL;
      }

      /// CPPONLY
      void setExpr(const string& expr="")
      {
        compileExpr(expr);
      }

      /// CPPONLY
      void setStmts(const string& stmts="")
      {
        compileStmts(stmts);
      }

      /// CPPONLY evaluate with string output
      PyObject* evaluate();

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
      void compileExpr(const string& expr);

      /// compile statements into byte code
      void compileStmts(const string& stmts);

      /// compiled code
      PyObject* m_expr, *m_stmts;

      /// local namespace
      PyObject* m_locals;
  };

  //////////////////////////////////////////////////////////////
  /// Stream element, can be of different types
  //////////////////////////////////////////////////////////////

  /// CPPONLY element of persistent stream
  class StreamElem
  {
    public:
      /// type of persistent string:
      /// 1. output only file stream
      /// 2. read/write file stream
      /// 3. string stream for maximum performance.
      enum streamType{ OFSTREAM, FSTREAM, SSTREAM };

      /// CPPONLY create a stream
      /// \param name filename
      /// \param readable iostream or just ostream
      /// \param realAppend whether or not keep old content when open an existing file
      /// \param useString use a stringstream rather than a file.
      StreamElem(const string& name, bool readable, bool realAppend, bool useString);

      /// CPPONLY copy constructor, called by map generator
      StreamElem(const StreamElem& rhs);

      /// CPPONLY destructor. Close stream and delete m_stream pointer.
      ~StreamElem();

      /// CPPONLY if the file is opened as write-only, chage to read-write mode using this function
      void makeReadable();

      /// CPPONLY set realAppend status. clear existing file if necessary
      void makeAppend(bool append);

      /// CPPONLY return the stream
      ostream* stream()
      {
        return m_stream;
      }

      /// CPPONLY return type of the stream
      streamType type()
      {
        return m_type;
      }

      /// CPPONLY output info (filename, type and file status)
      /// mostly for debug purposes
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
      ostream* m_stream;

      /// file name.
      string m_filename;
  };

  //////////////////////////////////////////////////////////////
  /// Stream Manager
  //////////////////////////////////////////////////////////////

  /** CPPONLY stream manager:
  This class keeps track of persistent ostreams, aceesible by
  their names.
  */
  class OstreamManager
  {
    public:

      /// CPPONLY
      OstreamManager();

      /// CPPONLY
      ~OstreamManager();

      /// CPPONLY get an ostream pointer from a name.
      /// if the stream does not exist, create one and return.
      ostream* getOstream( const string& name, bool readable,  bool realAppend, bool useString);

      /// CPPONLY if persistant ostream exist for a filename
      /// this is mostly for debug purposes
      bool hasOstream(const string& filename );

      /// CPPONLY list all registered ostreams
      void listAll();

      /// CPPONLY close all files and clean the map
      /// \param closeAppend if close all files including those persistent after evolution
      void closeAll(bool closeAppend);

    private:

      typedef map<string, StreamElem> ostreamMap;
      typedef map<string, StreamElem>::iterator ostreamMapIterator;
      typedef map<string, StreamElem>::value_type ostreamMapValue;

      /// collection of ostreams
      ostreamMap m_ostreams;
  };

  /// CPPONLY return a global Ostream Manager.
  OstreamManager& ostreamManager();

  //////////////////////////////////////////////////////////////
  /// OStream wrapper with our name extensions: >, >> etc
  //////////////////////////////////////////////////////////////

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
      StreamProvider(const string& output, const string& outputExpr);

      /// CPPONLY destructor
      ~StreamProvider(){}

      /// CPPONLY reset format string
      void setOutput(const string& output, const string& outputExpr);

      /// CPPONLY m_flag: if type ''
      bool noOutput(){ return ISSETFLAG(m_flags, m_flagNoOutput); }

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
      \param group group index of calling replicate

      these three parameters may be used to determine filename (when %grp etc
      are used in filename specification.

      */
      ostream& getOstream(PyObject*dict = NULL, bool readable = false);

      /// CPPONLY close ostream and delete ostream pointer.. if it is a ofstream.
      void closeOstream();

    private:

      /// we need to set
      ///  m_flagNoOutput
      ///  m_flagUseDefault (when no file is specified.)
      ///  m_flagAppend  (>> )
      ///  m_flagRealAppend (>>>)
      ///  m_flagUseString (|)
      void analyzeOutputString(const string& output);

      /// m_flag: if use and close type
      bool closeAfterUse(){ return ISSETFLAG(m_flags, m_flagCloseAfterUse); }

    private:
      /// plain text filename
      string m_filename;

      /// output filename parser
      Expression m_filenameExpr;

      /// internal m_flags of the operator. They are set during initialization for
      /// performance considerations.
      static const size_t  m_flagNoOutput      = 1;
      static const size_t  m_flagUseDefault    = 2;
      static const size_t  m_flagAppend        = 4;
      static const size_t  m_flagRealAppend    = 8;
      static const size_t  m_flagCloseAfterUse = 16;
      static const size_t  m_flagUseString     = 32;
      static const size_t  m_flagReadable      = 64;

      /// m_flags
      unsigned char m_flags;

      /// file stream pointer, used ONLY if we create a new file stream
      ostream* m_filePtr;
  };

  //////////////////////////////////////////////////////////////
  /// Random number generator
  //////////////////////////////////////////////////////////////

  /// for output purpose only.
  vectorlu bitSet(const BitSet&);

  /* \brief random number generator, functor, initializer etc,

    They can be individual functions but I will need a
    destructor to clear RNG structure.

    This implementation uses GSL rather than boost:;random because
    GSL has more distributions and easier/wide choice of RNGs.

    Right now, I can using the gsl library and will copy individual files
    to simupop directory later for easy distribution.
  */
  class RNG
  {

    public:
      ///
      RNG(const char* rng=NULL, unsigned long seed=0);

      ///
      ~RNG();

      ///
      /// choose an random number generator.
      /// This can be done by setting GSL_RNG_TYPE as well.
      void setRNG(const char * rng=NULL, unsigned long seed=0);

      /// return RNG name
      const char* name()
      {
        return gsl_rng_name(m_RNG);
      }

      unsigned long int max()
      {
        return gsl_rng_max(m_RNG);
      }

      string __repr__()
      {
        return toStr("<simuPOP::RNG ") + gsl_rng_name(m_RNG) + ">";
      }

      /// return min, 1, 2, ... max
      unsigned long int randGet()
      {
        return gsl_rng_get(m_RNG);
      }

      /// return [0, 1, 2, ... n-1]
      unsigned long int randInt(unsigned long int n)
      {
        return gsl_rng_uniform_int(m_RNG, n);
      }

      /// return [0, 1, 2, ... n-1]
      void randIntArray(ULONG n, ULONG size, ULONG* vec)
      {

        DBG_FAILIF( n <= 0 , ValueError, "RandInt: n should be positive.");

        for(size_t i=0; i<size; ++size)
          vec[i] =  gsl_rng_uniform_int(m_RNG, n);
      }

      /// return geometric with parameter p (k>=1)
      int randGeometric(double p)
      {
        return gsl_ran_geometric(m_RNG, p);
      }

      /// return double [0,1)
      double randUniform01()
      {
        return gsl_rng_uniform(m_RNG);
      }

      /// return double -inf, inf, v is standard deviation
      double randNormal(double m, double v)
      {
        return gsl_ran_gaussian(m_RNG, v) + m;
      }

      /// return double [0,1)
      void randUniform01Array(ULONG size, double * vec)
      {
        for(size_t i=0; i<size; ++size)
          vec[i] =  gsl_rng_uniform(m_RNG);
      }

      /// binomial distribution
      UINT randBinomial(UINT n, double p)
      {
        DBG_FAILIF( n <= 0 , ValueError, "RandInt: n should be positive.");

        return gsl_ran_binomial(m_RNG, p, n);
      }

      /// Poisson distribution
      UINT randPoisson( double p)
      {
        return gsl_ran_poisson(m_RNG, p);
      }

    private:

      /// global random number generator
      gsl_rng * m_RNG;
  };

  // weighted sampling using Walker's alias algorithm
  class WeightedSampler
  {
    public:
      // set up AliasMethod table
      WeightedSampler(RNG& rng, const vectorf& weight=vectorf(), bool fast=true )
        :m_RNG(&rng), m_q(0), m_a(0)              // , m_fast(fast)
      {
        set(weight);
      };

      ~WeightedSampler(){};

      void set(const vectorf& weight);

      ULONG biSearch(double a)
      {
        unsigned long left=0, right=m_N-1;

        // a reasonable initial guess will improve the performance a little bit
        unsigned long mid = static_cast<ULONG>(a*(m_N-1));

        // at most loop 4 times for n=31
        while( left <= right )
        {

          if( m_q[mid] > a )
          {
            if( mid==0 || m_q[mid-1] <= a)
              return mid;
            else
              right = mid - 1;
          }
          else
          {
            if( m_q[mid+1] > a )
              return mid+1;
            else
              left = mid + 1;
          }
          mid = (right+left)/2;
        }
        // this should not be reached.
        return 0;
      }

      // sample without replacement from 0,...,n-1,
      // with weight freq
      ULONG get()
      {
        //        if(m_fast)
        //      {
        double rN = m_RNG->randUniform01() * m_N;

        size_t K = static_cast<size_t>(rN);

        rN -= K;

        if( rN < m_q[K] )
          return K;
        else
          return  m_a[K];
        //      }
        //    else
        //   {
        //     return biSearch( m_RNG->randUniform01());
        //   }
      }

      // sample without replacement from 0,...,n-1,
      // with weight freq
      ULONG get(vectorlu& res, ULONG shift=0)
      {
        //        if(m_fast)
        //    {
        double rN;
        size_t K;

        for( vectorlu::iterator it=res.begin(); it != res.end(); ++it)
        {
          rN = m_RNG->randUniform01() * m_N;
          K = static_cast<size_t>(rN);

          rN -= K;

          if( rN < m_q[K] )
            *it = K+shift;
          else
            *it = m_a[K]+shift;
        }
        return 0;
        //  }
        //  else
        //  {
        //   for( vectorlu::iterator it=res.begin(); it != res.end(); ++it)
        //     *it = biSearch(m_RNG->randUniform01())+shift;
        //  }
      }

      template<class Iterator>
        ULONG get( Iterator beg, Iterator end,  ULONG shift=0)
      {
        //        if(m_fast)
        //       {
        double rN;
        size_t K;

        for( Iterator it=beg; it != end; ++it)
        {
          rN = m_RNG->randUniform01() * m_N;
          K = static_cast<size_t>(rN);

          rN -= K;

          if( rN < m_q[K] )
            *it = K+shift;
          else
            *it = m_a[K]+shift;
        }
        return 0;
        //     }
        //   else
        //   {
        //     for( Iterator it=beg; it != end; ++it)
        //       *it = biSearch(m_RNG->randUniform01())+shift;
        //   }
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
      RNG* m_RNG;

      /// length of weight.
      size_t m_N;

      /// internal table.
      vectorf m_q;

      vectorlu m_a;

      /// whether or not use fast (but use more RAM) method
      // bool m_fast;
  };

  /// CPPONLY
  /// this class encapsulate behavior of a sequence of Bernulli trial.
  /// the main idea is that when doing a sequence of Bernulli trials
  /// of the same probability, we can first use a binomial distribution
  /// to calculate the number of trials and then choose events randomly.
  ///
  /// For example, when N=10000, p=0.001. The usual way to do N Bin(p)
  /// trials is to do N randUnif(0,1)<p comparison.
  ///
  /// using the new method, we can generate number of success M~100
  /// and then find M locations.
  ///
  /// This class maintain a two dimensional table:
  /// a vector of probabilities cross expected number of trials
  ///
  ///           p1 p2 p3 p4 p5
  /// trial 1
  /// trial 2
  /// ...
  /// trial N
  ///
  /// We expect that N is big (usually populaiton size) and p_i are small
  ///
  /// using fast BernulliTrial method for fix p,
  /// we can fill up this table very quickly column by column
  ///
  /// This class will provide easy access to row (each trial) or column
  /// of this table.
  ///
  /// if this table is accessed row by row (each trial), a internal index
  /// is used.
  ///
  /// if index exceeds N, trials will be generated all again.
  /// if trial will be called, e.g., N+2 times all the time,
  /// this treatment might not be very efficient.
  ///
  class BernulliTrials
  {
    public:

      /// CPPONLY
      BernulliTrials(RNG& rng);

      ///
      /// CPPONLY
      BernulliTrials(RNG& rng, const vectorf& prob, ULONG trials);

      ///
      /// CPPONLY
      ~BernulliTrials();

      /// CPPONLY
      /// return size of trial
      ULONG size() const
      {
        return m_N;
      }

      vectorf prob() const
      {
        return m_prob;
      }

      /// CPPONLY
      void setParameter(const vectorf& prob, ULONG trials);

      /// generate the trial table, reset m_cur
      /// CPPONLY
      void doTrial();

      /// CPPONLY
      UINT curTrial();

      /// return the result of m_cur trial.
      /// if necessary, do trail again.
      /// CPPONLY
      const BitSet& trial();

      /// return succeed trials for p[index]
      /// fail when m_cur is not 0. (i.e., has retrieve the table through trial()
      /// CPPONLY
      const BitSet& succ(UINT index);

      /// CPPONLY
      vectorf probabilities()
      {
        return m_prob;
      }

      /// CPPONLY
      ULONG numTrials()
      {
        return m_N;
      }

    private:
      /// pointer to a random number generator.
      /// this is in preparation for multiple thread/RNG.
      RNG* m_RNG;

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

      /// current trial. Used when user want to access the table row by row
      ULONG m_cur;

      /// result to be returned.
      /// if we do need a result in bitset format. (each row)
      /// return this guy.
      BitSet m_bitSet;

  };

  /// currently, return a global RNG.
  RNG& rng();

  /// set random number generator
  void setRNG(const string rng="", unsigned long seed=0);

  /// list all available RNG.
  vectorstr listAllRNG();

  //////////////////////////////////////////////////////////////
  ///  Global debug and initialization related functions
  //////////////////////////////////////////////////////////////

  /// return svn revision
  int revision();

  /// return version infomation of simuPOP
  string simuVer();

  /// CPPONLY
  /// initialize module simuPOP when using "import simuPOP"
  bool initialize();

  bool optimized();

  

  /// CPPONLY
  /// get a null stream that discard everything
  ostream& cnull();

  /// CPPONLY
  /// set standard output to (default standard Python output)
  void setLogOutput(const string filename="");
}


/// CPPONLY
/// combine a pointer with its type obj and return a new PyObject*
/// this function will be appended to simuPOP_wrap.cpp
extern "C" PyObject* pointer2pyObj(void* p, char* typeString);
#endif
