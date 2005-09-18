/**************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu                                                        *
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

#ifndef _STATOR_H
#define _STATOR_H
/**
\file
\brief head file of class Stator:public Operator
*/
#include "utility.h"
#include "population.h"
#include "operator.h"

#include <string>
#include <numeric>
using std::string;
using std::count_if;
using std::bind2nd;
using std::equal;
using std::greater;
using std::min;
using std::max;

#include <iomanip>
using std::setprecision;

// for eval of string
#include "Python.h"

namespace simuPOP
{
  /// NOTE: the default output for stator is "", i.e., no output
  /// i.e., stator will write to shared variables and
  /// unless specified by output=">" etc, no output will be generated.
  ///
  /// this class will also list ALL statistics and its names?
  ///
  template<class Pop>
    class Stator: public Operator<Pop>
  {
    public:
      /// constructor. default to be always active.
      /// default to have NO output (shared variables will be set.)
      /// phase: if we treat Aa!=aA, default is false, i.e, Aa=aA.
      Stator(string output="", string outputExpr="",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL, string sep="\t"):
      Operator<Pop>(output, outputExpr, stage, begin, end, step, at, rep, grp, sep)
      {
      };

      /// destructor
      virtual ~Stator()
      {
      }

      /// this function is very important
      virtual Operator<Pop>* clone() const
      {
        return new Stator<Pop>(*this);
      }

  };

  /// evaluate an expression.
  ///
  template<class Pop>
    class PyEval: public Stator<Pop>
  {
    public:
      /** \brief evaluate expr/statments in local replicate namespace

      \param expr expression to be evaluated. Its result will be sent to output.
      \param stmts statement that will be executed before the expression.
      \param preStmts statement that will be executed when the operaot is
      constructed.
      \param postStmts statement that will be executed when the operator is
      destroyed.
      \param exposePop if true, expose current pop as variable "pop"
      \param name used to let pure python operator to identify themselves
      \param output default to ">" . i.e., output to standard output.
      for usage of other parameters, see help(baseOperator.__init__)
      */
      PyEval(const string& expr="", const string& stmts="", const string& preStmts="",
        const string& postStmts="", bool exposePop=false, const string& name="",
        string output=">", string outputExpr="",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL,  int grp=GRP_ALL, string sep="\t")
        :Stator<Pop>(output, outputExpr, stage, begin, end, step, at, rep, grp, sep),
        m_expr(expr, stmts), m_postExpr("", postStmts), m_exposePop(exposePop), m_name(name)
      {
        if( preStmts != "")
          Expression("", preStmts).evaluate();
      }

      ~PyEval()
      {
        m_postExpr.evaluate();
      }

      /// this function is very important
      virtual Operator<Pop>* clone() const
      {
        return new PyEval<Pop>(*this);
      }

      // check all alleles in vector allele if they are fixed.
      virtual bool apply(Pop& pop)
      {

        if(m_exposePop)
        {
          PyObject* popObj=pointer2pyObj((void*)(&pop),
            PopSWIGType);
          if( popObj == NULL)
            throw SystemError("Could not expose population pointer. Compiled with the wrong version of SWIG? ");

          Py_INCREF(popObj);
          // set dictionary variable pop to this object
          pop.setVar("pop", popObj);
        }

        m_expr.setLocalDict(pop.dict());
        string res = m_expr.valueAsString();

        if( ! this->noOutput() )
        {
          ostream & out = this->getOstream(pop.dict());
          out << res;
          this->closeOstream();
        }
        return true;
      }

      virtual string __repr__()
      {
        return "<simuPOP::PyEval " + m_name + ">";
      }

      string& name()
      {
        return m_name;
      }

    private:
      /// expression to evaluate
      Expression m_expr, m_postExpr;

      vectorstr m_headers;

      /// if expose pop
      bool m_exposePop;

      string m_name;
  };

  /// evaluate an expression.
  ///
  template<class Pop>
    class PyExec: public PyEval<Pop>
  {
    public:
      /** \brief evaluate statments in local replicate namespace, no return value

      \param stmts statements (a single or multi-line string) that will be evaluated before the expression.
      \param preStmts statement that will be executed when the operaot is
      constructed.
      \param postStmts statement that will be executed when the operator is
      destroyed.
      \param exposePop if true, expose current pop as variable "pop"
      \param output default to ">" . i.e., output to standard output.
      for usage of other parameters, see help(baseOperator.__init__)
      */
      PyExec( const string& stmts="", const string & preStmts="", const string& postStmts="",
        bool exposePop=false, const string& name="",
        string output=">", string outputExpr="",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL,  int grp=GRP_ALL, string sep="\t")
        :PyEval<Pop>("", stmts, preStmts, postStmts, exposePop, name, "", "",
        stage, begin, end, step, at, rep, grp, sep)
      {
      }

      /// this function is very important
      virtual Operator<Pop>* clone() const
      {
        return new PyExec<Pop>(*this);
      }

      virtual string __repr__()
      {
        return "<simuPOP::PyExec " + this->name() + ">";
      }
  };

  ///
  /// The following classes apply various statistics
  /// and basicStat class will provide an opearator interface
  /// to all of them.

  ///
  /// each class defines how to apply the statistics and
  /// provide interface to allow others to retrieve the value.
  /// NOTE: the values are population dependent so these
  /// values are only meaningful within the same basicStat.apply
  /// call;
  ///
  /// The design separate the calculation of statistics
  /// as much as possible and allow easier debug and writting
  /// new statistics.

  /// CPPONLY
  /// return {'a-b-b'} for a b c
  string haploKey(const vectori& seq)
  {
    string key = "{'"+toStr(seq[0]);

    for(size_t i=1; i< seq.size(); ++i)
      key += toStr("-") + toStr(seq[i]);

    return key + "'}";
  }

  /// CPPONLY
  /// post population sizes etc.
  template<class Pop>
    class statPopSize
  {
    private:

#define  numSubPop_String   "numSubPop"
#define  popSize_String     "popSize"
#define  subPopSize_String  "subPopSize"

    public:
      statPopSize(bool popSize=false)
        :m_isActive(popSize)
      {
      }

      void activate()
      {
        m_isActive = true;
      }

      void apply(Pop& pop)
      {
        UINT numSP = pop.numSubPop();
        ULONG popSize = pop.popSize();

        pop.setIntVar(numSubPop_String, numSP);
        pop.setIntVar(popSize_String, popSize);

        if( !m_isActive)
          return;

        // type mismatch, can not use subPopSizes() directly.
        vectori spSize(numSP);
        for(size_t sp=0; sp < numSP; ++sp)
          spSize[sp] = pop.subPopSize(sp);

        pop.setIntNumArrayVar(subPopSize_String, numSP, &spSize[0] );

        for( size_t sp=0; sp < numSP; ++sp)
          pop.setIntVar(subPopVar_String(sp, popSize_String), spSize[sp]);
      }

    private:
      bool m_isActive;
  };

  /// CPPONLY operator to tell the sex of an individual
  template<class IndType>
    class isMale
  {
    public:
      isMale(){};

      bool operator() ( const IndType& ind)
      {
        return ind.sex() == Male;
      }
  };

  /// CPPONLY
  template<class Pop>
    class statNumOfMale
  {
    private:

#define  numOfMale_String    "numOfMale"
#define  numOfFemale_String  "numOfFemale"

    public:

      statNumOfMale(bool numOfMale=false)
        :m_numOfMale(numOfMale?1:0), m_numOfFemale(0)
      {
      }

      void activate(bool yes=true)
      {
        m_numOfMale.resize(yes?1:0);
        m_numOfFemale.resize(yes?1:0);
      }

      ULONG numOfMale()
      {
        DBG_ASSERT( m_numOfMale.size() > 1, ValueError,
          "num of male has not been counted.");

        return m_numOfMale[ m_numOfMale.size() -1 ];
      }

      ULONG numOfFemale()
      {
        DBG_ASSERT( m_numOfFemale.size() > 1, ValueError,
          "num of female has not been counted.");

        return m_numOfFemale[ m_numOfFemale.size() -1 ];
      }

      ULONG numOfMale(UINT subPop)
      {
        DBG_ASSERT( m_numOfMale.size() > 1, ValueError,
          "num of male has not been counted.");

        DBG_ASSERT( subPop >= m_numOfMale.size() - 1, ValueError,
          "subPop index out of range.");

        return m_numOfMale[ subPop ];
      }

      ULONG numOfFemale(UINT subPop)
      {
        DBG_ASSERT( m_numOfFemale.size() > 1, ValueError,
          "num of male has not been counted.");

        DBG_ASSERT( subPop >= m_numOfFemale.size() - 1, ValueError,
          "subPop index out of range.");

        return m_numOfFemale[ subPop ];
      }

      void apply(Pop& pop)
      {
        if(m_numOfMale.empty())
          return;

        UINT numSP = pop.numSubPop();
        m_numOfMale.resize(numSP+1);
        m_numOfFemale.resize(numSP+1);

        ULONG numOfMale=0;

        for( size_t sp=0; sp < numSP; ++sp)
        {
          ULONG n = count_if(pop.indBegin( sp ), pop.indEnd( sp ),
            isMale<typename Pop::IndType>());
          numOfMale += n;
          m_numOfMale[sp] = n;

          pop.setIntVar(subPopVar_String(sp, numOfMale_String), n);

          n = pop.subPopSize(sp) - n;
          m_numOfFemale[sp] = n;

          pop.setIntVar( subPopVar_String(sp, numOfFemale_String), n);
        }
        pop.setIntVar( numOfMale_String, numOfMale);
        pop.setIntVar( numOfFemale_String, pop.popSize() - numOfMale);
        m_numOfMale[numSP] = numOfMale;
        m_numOfFemale[numSP] = pop.popSize() - numOfMale;
      }

    private:
      /// whether or not to apply number of male/female
      vectorlu m_numOfMale, m_numOfFemale;
  };

  /// CPPONLY
  template<class Pop>
    class statNumOfAffected
  {
    private:

#define  numOfAffected_String    "numOfAffected"
#define  propOfAffected_String   "propOfAffected"
#define  numOfUnaffected_String  "numOfUnaffected"
#define  propOfUnaffected_String  "propOfUnaffected"

    public:

      statNumOfAffected(bool numOfAffected=false)
        : m_numOfAffected(numOfAffected?1:0), m_numOfUnaffected(0)
      {
      }

      void activate(bool yes=true)
      {
        m_numOfAffected.resize(yes?1:0);
        m_numOfUnaffected.resize(yes?1:0);
      }

      ULONG numOfAffected()
      {
        DBG_ASSERT( m_numOfAffected.size() > 1, ValueError,
          "num of affected has not been counted.");

        return m_numOfAffected[ m_numOfAffected.size() -1 ];
      }

      ULONG numOfUnaffected()
      {
        DBG_ASSERT( m_numOfUnaffected.size() > 1, ValueError,
          "num of unaffected has not been counted.");

        return m_numOfUnaffected[ m_numOfUnaffected.size() -1 ];
      }

      ULONG numOfAffected(UINT subPop)
      {
        DBG_ASSERT( m_numOfAffected.size() > 1, ValueError,
          "num of affected has not been counted.");

        DBG_ASSERT( subPop >= m_numOfAffected.size() - 1, ValueError,
          "subPop index out of range.");

        return m_numOfAffected[ subPop ];
      }

      ULONG numOfUnaffected(UINT subPop)
      {
        DBG_ASSERT( m_numOfUnaffected.size() > 1, ValueError,
          "num of unaffected has not been counted.");

        DBG_ASSERT( subPop >= m_numOfUnaffected.size() - 1, ValueError,
          "subPop index out of range.");

        return m_numOfUnaffected[ subPop ];
      }

      void apply(Pop& pop)
      {
        if( m_numOfAffected.empty() )
          return;

        ULONG numOfAffected=0;
        UINT numSP = pop.numSubPop();
        m_numOfAffected.resize(numSP+1);
        m_numOfUnaffected.resize(numSP+1);

        for( size_t sp=0; sp < numSP; ++sp)
        {
          ULONG n = count_if(pop.indBegin( sp ), pop.indEnd( sp ),
            isAffected<typename Pop::IndType>());
          numOfAffected += n;
          m_numOfAffected[sp] = n;
          pop.setIntVar(subPopVar_String(sp, numOfAffected_String), n);
          pop.setDoubleVar( subPopVar_String(sp, propOfAffected_String),
            (double)(n)/pop.subPopSize(sp));

          n = pop.subPopSize(sp) - n;
          m_numOfUnaffected[sp] = n;

          pop.setIntVar(subPopVar_String(sp, numOfUnaffected_String), n);
          pop.setDoubleVar( subPopVar_String(sp, propOfUnaffected_String),
            (double)(n)/pop.subPopSize(sp));
        }
        pop.setIntVar( numOfAffected_String, numOfAffected);
        pop.setIntVar( numOfUnaffected_String, pop.popSize() - numOfAffected);
        pop.setDoubleVar( propOfAffected_String, (double)(numOfAffected)/pop.popSize());
        pop.setDoubleVar( propOfUnaffected_String, (double)(pop.popSize() - numOfAffected)/pop.popSize());
        m_numOfAffected[numSP] = numOfAffected;
        m_numOfUnaffected[numSP] = pop.popSize() - numOfAffected;
      }

    private:
      /// record the result.
      vectorlu m_numOfAffected, m_numOfUnaffected;
  };

  /// CPPONLY
  template<class Pop>
    class statAlleleFreq
  {
    private:

#define  NumOfAlleles_String  "numOfAlleles"
#define  AlleleNum_String     "alleleNum"
#define  AlleleFreq_String    "alleleFreq"

    public:

      statAlleleFreq( const vectori& atLoci = vectori())
        :m_atLoci(atLoci), m_ifPost(atLoci.size()), m_numOfAlleles(0),
        m_alleleNum(0), m_alleleFreq(0)
      {
        for(size_t i=0; i<atLoci.size(); ++i)
          m_ifPost[i] = 1;                        // true, post result
      }

      void addLocus(int locus, bool post=false)
      {
        if( find(m_atLoci.begin(), m_atLoci.end(), locus) == m_atLoci.end() )
        {
          m_atLoci.push_back( locus );
          // do not post result.
          m_ifPost.push_back(static_cast<int>(post));
        }
      }

      intMatrix& alleleNumAll()
      {
        return m_alleleNum.back();
      }

      vectori& alleleNumVec(int loc)
      {
        return  m_alleleNum.back()[loc];
      }

      int alleleNum(UINT allele, int loc)
      {
        // test for size?
        vectori& an= m_alleleNum.back()[loc];
        return allele < an.size() ? an[allele] : 0;
      }

      matrix& alleleFreqAll()
      {
        return m_alleleFreq.back();
      }

      vectorf& alleleFreqVec(int loc)
      {
        return m_alleleFreq.back()[loc];
      }

      double alleleFreq(UINT allele, int loc)
      {
        vectorf& af = m_alleleFreq.back()[loc];
        return allele < af.size() ? af[allele] : 0;
      }

      intMatrix& alleleNumAll(UINT subPop)
      {
        return m_alleleNum[subPop];
      }

      vectori& alleleNumVec(int loc, UINT subPop)
      {
        return m_alleleNum[subPop][loc];
      }

      int alleleNum(UINT allele, int loc, UINT subPop)
      {
        vectori& an = m_alleleNum[subPop][loc];
        return allele < an.size() ? an[allele] : 0;
      }

      matrix& alleleFreqAll(UINT subPop)
      {
        return m_alleleFreq[subPop];
      }

      vectorf& alleleFreqVec(int loc, UINT subPop)
      {
        return m_alleleFreq[subPop][loc];
      }

      double alleleFreq(UINT allele, int loc, UINT subPop)
      {
        vectorf& af = m_alleleFreq[subPop][loc];
        return allele < af.size() ? af[allele] : 0.;
      }

      vectori& numOfAlleles()
      {
        UINT idx=m_numOfAlleles.size()==2?0:(m_numOfAlleles.size()-1);
        return m_numOfAlleles[idx];
      }

      vectori& numOfAlleles(UINT subPop)
      {
        DBG_ASSERT( subPop < m_numOfAlleles.size()-1, IndexError,
          "Subpop index " + toStr(subPop) + " out of range of 0 ~ "
          +toStr(m_numOfAlleles.size()-2));
        return m_numOfAlleles[subPop];
      }

      vectori alleles(int loc)
      {
        vectori al;

        for(size_t j=1; j < m_alleleNum.back()[loc].size(); ++j)
          if( m_alleleNum.back()[loc][j] > 0)
            al.push_back(j);

        DBG_ASSERT( al.size() == static_cast<UINT>(numOfAlleles()[loc]),
          SystemError, "Number of alleles at locus " + toStr(loc)
          + " does not match.Observed "
          + toStr( al.size()) + " previous count: "
          + toStr( numOfAlleles()[loc] ));

        return al;
      }

      vectori alleles(int loc, UINT subPop)
      {
        vectori al;

        for(size_t j=1; j < m_alleleNum[subPop][loc].size(); ++j)
          if( m_alleleNum[subPop][loc][j] != 0)
            al.push_back(j);

        DBG_ASSERT( al.size() == static_cast<UINT>(numOfAlleles(subPop)[loc]),
          SystemError, "Number of alleles at locus " + toStr(loc)
          + " at subpop " + toStr(subPop) + " does not match. Observed "
          + toStr( al.size()) + " previous count: "
          + toStr( numOfAlleles(subPop)[loc] ));

        return al;
      }

      void apply(Pop& pop)
      {
        if( m_atLoci.empty())
          return;

        UINT numSP = pop.numSubPop();
        UINT numLoci = m_atLoci.size();

        UINT len = numSP==1?1:(numSP+1);
        // if not initialized or m_atLoci/numSP changes
        if( m_alleleNum.size() != len )
        {
          for( size_t i=0; i < numLoci;  ++i)
            DBG_FAILIF( static_cast<UINT>(m_atLoci[i]) >= pop.totNumLoci(),
              IndexError, "locus index (" + toStr(m_atLoci[i])
              + ") out of range of 0 - " + toStr(pop.totNumLoci()-1));

          m_alleleNum.resize(len);
          m_alleleFreq.resize(len);

          m_numOfAlleles.resize(len);
          for(size_t i = 0; i < len; ++i)
          {
            m_alleleNum[i].resize(pop.totNumLoci());
            m_alleleFreq[i].resize(pop.totNumLoci());
            m_numOfAlleles[i].resize(pop.totNumLoci());
            fill(m_numOfAlleles[i].begin(), m_numOfAlleles[i].end(), 0);
          }
        }

        string varname;

        for(size_t i = 0; i < numLoci; ++i)
        {
          UINT loc = m_atLoci[i];

          vectori& sum = m_alleleNum.back()[loc];
          fill(sum.begin(), sum.end(), 0);

          // for each subPopulation
          for( UINT sp=0; sp < numSP;  ++sp)
          {
            vectori& num = m_alleleNum[sp][loc];
            // clear all current values
            fill(num.begin(), num.end(),0);

            // go through all alleles
            for(typename Pop::AlleleIterator a=pop.alleleBegin(loc, sp),
              aEnd=pop.alleleEnd(loc, sp); a != aEnd; ++a)
            {
              if( *a >= num.size() )
                num.resize(*a+1, 0);
              num[*a]++;
            }

            // add this number to overall num
            // calculate frequency
            if(numSP > 1 )
            {
              if(sum.size() < num.size())
                sum.resize(num.size(), 0);
              for(size_t e=0, eEnd=num.size(); e < eEnd; ++e)
                sum[e] += num[e];
            }

            vectorf& freq = m_alleleFreq[sp][loc];
            freq.resize( num.size(), 0.);
            double dy = pop.subPopSize(sp)*pop.ploidy();
            for(size_t e=0, eEnd=num.size(); e < eEnd; ++e)
              freq[e] = static_cast<double>(num[e])/dy;

            if(m_ifPost[i])
            {
              varname = subPopVar_String(sp, AlleleNum_String) + "[" + toStr(loc) + "]";
              PyObject * d = pop.setIntNumArrayVar(varname, num.size(), &num[0]);

              if(numSP == 1)
              {
                varname = toStr(AlleleNum_String) + "[" + toStr(loc) + "]";
                Py_INCREF(d);
                pop.setVar(varname, d);
              }

              varname = subPopVar_String(sp, AlleleFreq_String) + "[" + toStr(loc) + "]";
              d = pop.setDoubleNumArrayVar(varname, freq.size(), &freq[0]);

              if(numSP == 1)
              {
                varname = toStr(AlleleFreq_String) + "[" + toStr(loc) + "]";
                Py_INCREF(d);
                pop.setVar(varname, d);
              }
            }                                     // post

            // set numOfAlleles if necessary
            m_numOfAlleles[sp][loc] = count_if( num.begin(), num.end(),
              bind2nd(std::greater<int>(),0));
          }                                       // subpop

          if(numSP > 1 )                          // calculate sum and post overall result
          {
            // summary?
            vectorf& freq = m_alleleFreq.back()[loc];
            freq.resize( sum.size());
            double dy = pop.popSize() * pop.ploidy();
            for(size_t e=0, eEnd=sum.size(); e < eEnd; ++e)
              freq[e] = static_cast<double>(sum[e])/dy;

            if(m_ifPost[i])
            {
              varname = string(AlleleNum_String) + "[" + toStr(loc) + "]";
              pop.setIntNumArrayVar(varname, sum.size(), &sum[0]);

              varname = string(AlleleFreq_String) + "[" + toStr(loc) + "]";
              pop.setDoubleNumArrayVar(varname, freq.size(), &freq[0]);
            }

            // set numOfAlleles if necessary
            m_numOfAlleles.back()[loc] = count_if( sum.begin(), sum.end(),
              bind2nd(std::greater<int>(),0));
          }
        }                                         // all loci

        if( accumulate(m_ifPost.begin(), m_ifPost.end(), 0) > 0 )
        {
          // post number of alleles
          for(UINT sp = 0; sp < numSP; ++sp)
          {
            PyObject* d = pop.setIntNumArrayVar( subPopVar_String(sp, NumOfAlleles_String),
              pop.totNumLoci(), &m_numOfAlleles[sp][0]);
            if(numSP == 1)
            {
              Py_INCREF(d);
              pop.setVar(NumOfAlleles_String, d);
            }
          }
          if(numSP>1)
          {
            pop.setIntNumArrayVar( NumOfAlleles_String, pop.totNumLoci(),
              &m_numOfAlleles.back()[0]);
          }
        }
      }

    private:

      /// which alleles?
      vectori m_atLoci;

      /// whether or not post result
      vectori m_ifPost;

      /// count number of alleles
      intMatrix m_numOfAlleles;

      /// allele counter! use this matrix to avoid frequent
      /// allocation of memory
      vector<intMatrix> m_alleleNum;

      /// allele Freq
      vector<matrix> m_alleleFreq;
  };

  /// CPPONLY
  /// use alleleFreq to get number of alleles
  /// so statNumOfAlleles is just a proxy class
  ///
  template<class Pop>
    class statNumOfAlleles
  {
    public:
      statNumOfAlleles(statAlleleFreq<Pop>& calc, const vectori& atLoci = vectori())
        :m_calc(calc)
      {
        for(vectori::const_iterator it = atLoci.begin(); it != atLoci.end(); ++it)
          m_calc.addLocus(*it, true);
      }

      // do nothing. m_calc.spply will be called by basicStat.
      void apply(Pop& pop)
      {
      }

    private:

      /// a reference to an existing allelefreq calculator
      statAlleleFreq<Pop>& m_calc;
  };

  /// CPPONLY
  template<class Pop>
    class statHeteroFreq
  {
    private:

#define HeteroNum_String        "heteroNum"
#define HeteroFreq_String       "heteroFreq"
#define HomoNum_String          "homoNum"
#define HomoFreq_String         "homoFreq"

      int locusIdx(int loc)
      {
        UINT idx=0;
        while( m_atLoci[idx] != loc && idx < m_atLoci.size() )
          idx++;
        DBG_ASSERT( m_atLoci[idx] == loc, ValueError,
          "Can not find allele freq for locus " + toStr(loc));
        return idx;
      }

      // the result layout is
      //  loc 1, sup1
      //  loc 2, sup1
      //  ...
      //  loc 1 , sup2
      //  loc 2, sup2
      //  ...
      //  loc 1, summary
      //  loc 2, summary
      //  ...

      int resIdx(int idx)                         // for summary
      {
        UINT numSP = m_heteroNum.size()/m_atLoci.size() - 1;
        if(numSP==1)
          return idx;
        else
          return idx + numSP*m_atLoci.size();
      }

      int resIdx(int idx, UINT subPop)
      {
        DBG_ASSERT( subPop < m_heteroNum.size()/m_atLoci.size() - 1,
          IndexError,
          "Subpop index " + toStr(subPop) + " out of range of 0 ~ "
          +toStr( m_heteroNum.size()/m_atLoci.size() -2));

        return idx + subPop*m_atLoci.size();
      }

    public:
      statHeteroFreq(const vectori& heteroFreq=vectori(),
        const vectori& homoFreq=vectori())
        : m_atLoci(heteroFreq), m_ifPost(0),
        m_postHetero(!heteroFreq.empty()), m_postHomo(!homoFreq.empty()),
        m_heteroNum(0), m_heteroFreq(0), m_homoNum(0), m_homoFreq(0)
      {
        // add homofreq to m_atLoci
        for(size_t i=0; i < homoFreq.size(); ++i)
          if( find(m_atLoci.begin(), m_atLoci.end(), homoFreq[i]) == m_atLoci.end())
            m_atLoci.push_back( homoFreq[i]);

        m_ifPost.resize(m_atLoci.size(), 1);      // post result
      }

      void addLocus(int locus, bool post=false)
      {
        if( find(m_atLoci.begin(), m_atLoci.end(), locus) == m_atLoci.end() )
        {
          m_atLoci.push_back( locus );
          // default not post result
          m_ifPost.push_back(static_cast<int>(post));
        }
        m_postHetero = true;
      }

      int  heteroNum(UINT allele, int loc)
      {
        UINT idx=locusIdx(loc);

        vectori& hn =  m_heteroNum[resIdx(idx)];
        return allele < hn.size() ? hn[allele] : 0;
      }

      double  heteroFreq(UINT allele, int loc)
      {
        UINT idx=locusIdx(loc);
        vectorf& hf =  m_heteroFreq[resIdx(idx)];
        return allele < hf.size() ? hf[allele] : 0.;
      }

      int  heteroNum(UINT allele, int loc, UINT subPop)
      {
        UINT idx=locusIdx(loc);
        vectori& hn =  m_heteroNum[resIdx(idx,subPop)];
        return allele < hn.size() ? hn[allele] : 0;
      }

      double  heteroFreq(UINT allele, int loc, UINT subPop)
      {
        UINT idx=locusIdx(loc);
        vectorf& hf =  m_heteroFreq[resIdx(idx,subPop)];
        return allele < hf.size() ? hf[allele] : 0;
      }

      void apply(Pop& pop)
      {
        if( m_atLoci.empty())
          return;

        UINT numSP = pop.numSubPop();
        UINT numLoci = m_atLoci.size();

        // may be resizing for different replicate of populations.
        // if not initialized or m_atLoci/numSP changes
        if(m_heteroNum.size() != (numSP+1)*numLoci)
        {
          m_heteroNum.resize((numSP+1)*numLoci);
          m_heteroFreq.resize((numSP+1)*numLoci);
          m_homoNum.resize( numSP+1 );
          m_homoFreq.resize( numSP+1 );
        }

        string varname;
        ULONG popSize = pop.popSize();

        for( size_t i = 0; i < numLoci; ++i)
        {
          UINT loc = m_atLoci[i];
          DBG_DO(DBG_STAT, cout << "Counting heterozygotes at locus " << loc << endl);

          vectori& sum = m_heteroNum[resIdx(i)];
          fill(sum.begin(), sum.end(), 0);

          // for each subPopulation
          for( UINT sp=0; sp < numSP;  ++sp)
          {
            vectori& num = m_heteroNum[resIdx(i,sp)];
            fill(num.begin(), num.end(), 0 );

            // go through all alleles
            //?>> \todo here we assume diploid population
            for(typename Pop::AlleleIterator a=pop.alleleBegin(loc, sp),
              aEnd=pop.alleleEnd(loc, sp);
              a != aEnd; a+=2)
            {
              if( *a >= num.size() )
                num.resize(*a+1);

              if( *(a+1) >= num.size() )
                num.resize(*(a+1)+1);

              if(  *a != *(a+1) )                 // heterozygote
              {
                num[0] ++;                        // overall x != y
                num[*a]++;
                num[*(a+1)]++;
              }
            }

            // add this number to overall num
            // calculate frequency
            if(numSP > 1)
            {
              if(sum.size() < num.size())
                sum.resize(num.size());
              for(size_t e=0, eEnd=num.size(); e < eEnd; ++e)
                sum[e] += num[e];
            }

            vectorf& freq = m_heteroFreq[resIdx(i,sp)];
            freq.resize( num.size());
            for(size_t e=0, eEnd=num.size(); e < eEnd; ++e)
              freq[e] = static_cast<double>(num[e])/pop.subPopSize(sp);

            // set variable.
            if(m_ifPost[i] && m_postHetero)
            {
              varname =  subPopVar_String(sp, HeteroNum_String) + "[" + toStr(loc) + "]";
              PyObject* d = pop.setIntNumArrayVar(varname, num.size(), &num[0]);
              if(numSP == 1)
              {
                Py_INCREF(d);
                varname =  toStr(HeteroNum_String) + "[" + toStr(loc) + "]";
                pop.setVar(varname, d);
              }

              varname = subPopVar_String(sp, HeteroFreq_String) + "[" + toStr(loc) + "]";
              d = pop.setDoubleNumArrayVar(varname, freq.size(), &freq[0]);

              if(numSP == 1)
              {
                Py_INCREF(d);
                varname =  toStr(HeteroFreq_String) + "[" + toStr(loc) + "]";
                pop.setVar(varname, d);
              }
            }
          }

          if(numSP > 1 && m_postHetero)
          {
            vectorf& freq = m_heteroFreq[resIdx(i)];
            freq.resize( sum.size());
            for(size_t e=0, eEnd=sum.size(); e < eEnd; ++e)
              freq[e] = static_cast<double>(sum[e])/popSize;

            if(m_ifPost[i])
            {                                     //
              varname = string(HeteroNum_String) + "[" + toStr(loc) + "]";
              pop.setIntNumArrayVar(varname, sum.size(), &sum[0]);

              varname = string(HeteroFreq_String) + "[" + toStr(loc) + "]";
              pop.setDoubleNumArrayVar(varname, freq.size(), &freq[0]);
            }
          }                                       // whole population
        }                                         // for all loci

        if( m_postHomo )
        {
          for( size_t i = 0; i < numLoci; ++i)
          {
            UINT loc = m_atLoci[i];

            if(loc+1 >= m_homoFreq[0].size())
            {
              for( UINT sp=0; sp < numSP+1;  ++sp)
              {
                m_homoFreq[sp].resize(loc+1, 0.0);
                m_homoNum[sp].resize(loc+1, 0);
              }
            }

            // calculate homoNum
            for( UINT sp=0; sp < numSP;  ++sp)
            {
              m_homoNum[sp][loc] = pop.subPopSize(sp) - m_heteroNum[resIdx(i, sp)][0];
              m_homoFreq[sp][loc] = (double)(m_homoNum[sp][loc])/pop.subPopSize(sp);
            }
            m_homoNum[numSP][loc] = pop.popSize() - m_heteroNum[resIdx(i)][0];
            m_homoFreq[numSP][loc] = (double)(m_homoNum[numSP][loc])/pop.popSize();

          }                                       // all loci
          // post result
          for( UINT sp=0; sp < numSP; ++sp)
          {
            pop.setIntNumArrayVar(subPopVar_String(sp,HomoNum_String),
              m_homoNum[sp].size(), &m_homoNum[sp][0]);
            pop.setDoubleNumArrayVar(subPopVar_String(sp,HomoFreq_String),
              m_homoFreq[sp].size(), &m_homoFreq[sp][0]);
          }
          pop.setIntNumArrayVar(HomoNum_String,
            m_homoNum[numSP].size(), &m_homoNum[numSP][0]);
          pop.setDoubleNumArrayVar(HomoFreq_String,
            m_homoFreq[numSP].size(), &m_homoFreq[numSP][0]);
        }
      }

    private:

      /// heteroFreq
      vectori m_atLoci;

      ///
      vectori m_ifPost;

      bool m_postHetero, m_postHomo;

      /// hetero counter
      intMatrix m_heteroNum;

      /// hetero Freq
      matrix m_heteroFreq;

      /// expected heterozygosity
      matrix m_expHetero;

      /// homozygosity number and freq
      intMatrix m_homoNum;

      matrix m_homoFreq;
  };

  /// CPPONLY
  template<class Pop>
    class statExpHetero
  {
    private:

#define ExpHetero_String "expHetero"

    public:
      statExpHetero(statAlleleFreq<Pop>& alleleFreq, const vectori& expHetero=vectori())
        : m_alleleFreq(alleleFreq), m_atLoci(expHetero), m_expHetero(0)
      {
        // add expected hetero to m_alleleFreq
        for(size_t i=0; i < expHetero.size(); ++i)
          m_alleleFreq.addLocus( expHetero[i]);
      }

      void apply(Pop& pop)
      {
        if( m_atLoci.empty())
          return;

        UINT numSP = pop.numSubPop();
        UINT numLoci = m_atLoci.size();

        if(m_expHetero.size() != numSP+1)
          m_expHetero.resize(numSP+1 );

        for( size_t i = 0; i < numLoci; ++i)
        {
          UINT loc = m_atLoci[i];

          if(loc+1 >= m_expHetero[0].size())
          {
            for( UINT sp=0; sp < numSP+1;  ++sp)
              m_expHetero[sp].resize(loc+1, 0.0);
          }

          // for each subPopulation
          for( UINT sp=0; sp < numSP;  ++sp)
          {
            // calculate expected heterozygosity
            // get allele frequency
            vectorf& af = m_alleleFreq.alleleFreqVec(loc, sp);
            double expHeter=1;
            // 1-sum pi^2
            for(int al = 0, alEnd=af.size() ; al < alEnd; al++)
              expHeter -= af[al]*af[al];

            m_expHetero[sp][loc] = expHeter;
          }

          vectorf& af = m_alleleFreq.alleleFreqVec(loc);
          double expHeter=1;
          // 1-sum pi^2
          for(int al = 0, alEnd = af.size(); al < alEnd; al++)
            expHeter -= af[al]*af[al];

          m_expHetero[numSP][loc] = expHeter;
        }

        // post result
        for( UINT sp=0; sp < numSP; ++sp)
          pop.setDoubleNumArrayVar(subPopVar_String(sp,ExpHetero_String),
            m_expHetero[sp].size(), &m_expHetero[sp][0]);
        pop.setDoubleNumArrayVar(ExpHetero_String,
          m_expHetero[numSP].size(), &m_expHetero[numSP][0]);
      }

    private:

      /// need this to apply alleleFreq
      statAlleleFreq<Pop>& m_alleleFreq;

      /// heteroFreq
      vectori m_atLoci;

      /// expected heterozygosity
      matrix m_expHetero;
  };

  /// CPPONLY
  /// currently there is no need to expose the result.
  /// may add that later.
  template<class Pop>
    class statGenoFreq
  {
    private:

#define  GenotypeNum_String   "genoNum"
#define  GenotypeFreq_String  "genoFreq"

    public:
      statGenoFreq(const vectori& genoFreq = vectori(), bool phase=false )
        : m_atLoci(genoFreq), m_phase(phase)
      {
      }

      void apply(Pop& pop)
      {
        if( m_atLoci.empty())
          return;

        UINT numSP = pop.numSubPop();
        ULONG popSize = pop.popSize();

        for( size_t i=0, iEnd = m_atLoci.size(); i < iEnd;  ++i)
        {
          if(static_cast<UINT>(m_atLoci[i]) >= pop.totNumLoci() )
            throw IndexError("Absolute locus index "
              + toStr( m_atLoci[i]) + " is out of range of 0 ~ "
              + toStr( pop.totNumLoci() - 1));
        }

        // first remove genoNum that may be set by previous count.
        // for example if genoNum[a][10] was set but this time there
        // is no allele 10, we will not reset genoNum[a][10] ...
        pop.removeVar( GenotypeNum_String);
        pop.removeVar( GenotypeFreq_String);
        for( UINT sp=0; sp < numSP;  ++sp)
        {
          // remove genoNum, genoFreq that may be set by previous run.
          pop.removeVar( subPopVar_String(sp, GenotypeNum_String));
          pop.removeVar( subPopVar_String(sp, GenotypeFreq_String));
        }

        string varname;

        // deal with genotype
        for( size_t i=0, iEnd = m_atLoci.size(); i < iEnd;  ++i)
        {
          // for each locus, we need to use a vector of dictionaries.
          vector< intDict> sum;

          int loc = m_atLoci[i];

          Allele a, b;

          // for each subPopulation
          for( UINT sp=0; sp < numSP;  ++sp)
          {
            DBG_DO(DBG_STAT, cout << "Counting genotypes at locus " <<
              loc << " subPop " << sp << endl);

            vector< intDict> num;

            /// go through a single allele for all individual, all diploid
            for( typename Pop::AlleleIterator it = pop.alleleBegin(loc, sp),
              itEnd = pop.alleleEnd(loc, sp); it != itEnd;  it+=2 )
            {
              a = *it;
              b = *(it+1);

              if( !m_phase && a > b )
                std::swap(a,b);

              // count:
              // DBG_ASSERT( a < num.size(), SystemError,
              //  "Allele number " + toStr(int(a)) + " is greater than population maxAllele() (" + toStr(pop.maxAllele()) +")");

              if( a >= num.size() )
                num.resize(a+1);

              num[a][b]++;

              if( a >= sum.size() )
                sum.resize(a+1);

              sum[a][b]++;
            }

            // register values for this subpopulation
            for(a=1; a < num.size(); ++a)
            {
              /// need to replace previous values
              // if( num[a].empty() )
              //   continue;

              // empty dictionary should be allowed
              varname = subPopVar_String(sp, GenotypeNum_String) +
                + "[" + toStr(loc) + "][" + toStr(int(a)) + "]";
              pop.setIntDictVar( varname, num[a] );

              // apply frequency
              for(intDict::iterator it=num[a].begin(), itEnd=num[a].end(); it!=itEnd; ++it)
                it->second = it->second/ pop.subPopSize(sp);

              varname =  subPopVar_String(sp, GenotypeFreq_String) +
                + "[" + toStr(loc) + "][" + toStr(int(a)) + "]";
              pop.setIntDictVar( varname, num[a] );
            }
          }

          for(a=1; a<sum.size(); ++a)
          {
            // if( sum[a].empty() )
            //  continue;

            // empty dictionary should be allowed
            varname = toStr(GenotypeNum_String) + "[" + toStr(loc) + "][" + toStr(int(a)) + "]";
            pop.setIntDictVar( varname, sum[a] );

            // apply frequency
            for(intDict::iterator it=sum[a].begin(), itEnd=sum[a].end(); it!=itEnd; ++it)
              it->second = it->second/ popSize;

            varname =  toStr(GenotypeFreq_String) + "[" + toStr(loc) + "][" + toStr(int(a)) + "]";
            pop.setIntDictVar( varname, sum[a] );
          }
        }
      }

    private:
      /// which genotypes
      vectori m_atLoci;

      /// phase
      bool m_phase;
  };

  /// CPPONLY
  template<class Pop>
    class statHaploFreq
  {
    private:

#define  HaplotypeNum_String    "haploNum"
#define  HaplotypeFreq_String   "haploFreq"
      int haploIndex(const vectori& haplo)
      {
        // first locate haplo
        UINT idx = 0;
        while( m_haplotypes[idx] != haplo && idx < m_haplotypes.size())
          idx ++;

        DBG_ASSERT( idx != m_haplotypes.size(), ValueError,
          "Can not find haplotype.");

        return idx;
      }

    public:

      statHaploFreq(const intMatrix& haploFreq=intMatrix())
        :m_haplotypes(haploFreq), m_ifPost(haploFreq.size())
      {
        for(size_t i=0; i<haploFreq.size(); ++i)
          m_ifPost[i] = 1;
      }

      void addHaplotype(const vectori& haplo, bool post=false)
      {
        if( find(m_haplotypes.begin(), m_haplotypes.end(), haplo) == m_haplotypes.end() )
        {
          m_haplotypes.push_back( haplo );
          m_ifPost.push_back(static_cast<int>(post));
        }
      }

      int numOfHaplotypes(const vectori& haplo)
      {
        // first locate haplo
        UINT idx = haploIndex(haplo);

        return m_haploNum[idx + m_haploNum.size() -
          m_haplotypes.size()].size();
      }

      int numOfHaplotypes(const vectori& haplo, UINT subPop)
      {
        // first locate haplo
        UINT idx = haploIndex(haplo);

        return m_haploNum[idx + subPop* m_haplotypes.size()].size();
      }

      map<vectori, UINT>& haploNum(const vectori& haplo)
      {
        UINT idx = haploIndex(haplo);

        return m_haploNum[idx + m_haploNum.size() - m_haplotypes.size() ];
      }

      map<vectori, double>& haploFreq(const vectori& haplo)
      {
        UINT idx = haploIndex(haplo);

        return m_haploFreq[idx + m_haploNum.size() - m_haplotypes.size() ];
      }

      map<vectori, UINT>& haploNum(const vectori& haplo, UINT subPop)
      {
        UINT idx = haploIndex(haplo);

        return m_haploNum[idx + subPop*m_haplotypes.size() ];
      }

      map<vectori, double>& haploFreq(const vectori& haplo, UINT subPop)
      {
        UINT idx = haploIndex(haplo);

        return m_haploFreq[idx + subPop*m_haplotypes.size() ];
      }

      void apply(Pop& pop)
      {
        if( m_haplotypes.empty())
          return;

        UINT nHap = m_haplotypes.size();

        // first time?
        if( m_haploNum.size() != (pop.numSubPop()+1)*nHap)
        {
          for( size_t i = 0; i < nHap;  ++i)
          {
            vectori& haplotype = m_haplotypes[i];
            size_t sz = haplotype.size();

            if( sz == 0 )
              throw ValueError("has to specify some haplotype.");
            else if( sz == 1 )
              throw ValueError("Haplotype must contain alleles from more than one locus.");
          }
          m_haploNum.resize( (pop.numSubPop()+1)*nHap);
          m_haploFreq.resize( (pop.numSubPop()+1)*nHap);
        }

        UINT numSP = pop.numSubPop();

        DBG_DO(DBG_STAT, cout << "Counting haplotypes" << endl);

        // clear all statistics
        for(size_t h = 0; h < nHap*(numSP+1); ++h)
        {
          m_haploNum[h].clear();
          m_haploFreq[h].clear();
        }

        // for each subPopulation
        for( UINT sp=0; sp < numSP;  ++sp)
        {
          for(size_t h = 0; h < nHap; ++h)
          {
            vectori& haplotype = m_haplotypes[h];
            map< vectori, UINT>& count = m_haploNum[h + sp*nHap];
            map< vectori, UINT>& sum = m_haploNum[h + numSP*nHap];

            size_t sz = haplotype.size();

            vectori sampleHap(sz);

            for( typename Pop::AlleleIterator it = pop.alleleBegin(0, sp),
              itEnd = pop.alleleEnd(0, sp); it != itEnd;  ++it )
            {
              for( size_t hap=0; hap < sz; ++hap)
                sampleHap[hap] = *(it.ptr()+haplotype[hap]);

              // add sampleHap count
              count[sampleHap] ++;
              sum[sampleHap] ++;
            }
          }
        }
        // finish count

        // calculate haploFreq,
        // for each subPopulation
        for( UINT sp=0; sp < numSP;  ++sp)
        {
          // record both num and freq
          string varNumName =  subPopVar_String(sp, HaplotypeNum_String);
          string varFreqName =  subPopVar_String(sp, HaplotypeFreq_String);
          for(size_t h = 0; h < nHap; ++h)
          {
            vectori& haplotype = m_haplotypes[h];
            map< vectori, UINT>& count = m_haploNum[h + sp*nHap];
            map< vectori, double>& freq = m_haploFreq[h + sp*nHap];

            for(map<vectori, UINT>::iterator it = count.begin(); it!=count.end(); ++it)
            {
              freq[ it->first] = double(it->second)/(pop.subPopSize(sp)*pop.ploidy());
              if(m_ifPost[h])
              {
                pop.setIntVar( varNumName + haploKey(haplotype) + haploKey( it->first) , it->second);
                pop.setDoubleVar( varFreqName + haploKey(haplotype) + haploKey( it->first) ,
                  double(it->second)/(pop.subPopSize(sp)*pop.ploidy()));
              }
            }
          }
        }
        // whole population
        for(size_t h = 0; h < nHap; ++h)
        {
          vectori& haplotype = m_haplotypes[h];
          map< vectori, UINT>& count = m_haploNum[h + numSP*nHap];
          map< vectori, double>& freq = m_haploFreq[h+numSP*nHap];

          for(map<vectori, UINT>::iterator it = count.begin(); it!=count.end(); ++it)
          {
            freq[ it->first] = double(it->second)/(pop.popSize()*pop.ploidy());
            if(m_ifPost[h])
            {
              pop.setIntVar( HaplotypeNum_String + haploKey(haplotype) + haploKey( it->first) , it->second);
              pop.setDoubleVar( HaplotypeFreq_String + haploKey(haplotype) + haploKey( it->first) ,
                double(it->second)/(pop.popSize()*pop.ploidy()));
            }
          }
        }
      }

    private:
      /// which haplotypes
      intMatrix m_haplotypes;

      vectori m_ifPost;

      /// keep results
      vector< map< vectori, UINT> > m_haploNum;

      /// keep result
      vector< map< vectori, double> > m_haploFreq;
  };

  /// CPPONLY
  template<class Pop>
    class statLD
  {
    private:

#define   LD_String           "ld"
#define   LDPRIME_String      "ld_prime"
#define   R2_String           "r2"
#define   AvgLD_String        "LD"
#define   AvgLDPRIME_String   "LD_prime"
#define   AvgR2_String        "R2"

    public:

      statLD(statAlleleFreq<Pop>& alleleFreq, statHaploFreq<Pop>& haploFreq,
        const intMatrix& LD=intMatrix(), bool midValues=false)
        :m_alleleFreq(alleleFreq), m_haploFreq(haploFreq),
        m_LD(LD), m_midValues(midValues)
      {
        for( size_t i=0, iEnd = m_LD.size(); i < iEnd; ++i)
        {
          if( m_LD[i].size() != 2 && m_LD[i].size() != 4 )
            throw ValueError("Expecting [locus locus [allele allele ]] items");

          // what we will need?
          // 1. allele freq
          m_alleleFreq.addLocus( m_LD[i][0], midValues );
          m_alleleFreq.addLocus( m_LD[i][1], midValues );
          // 2. haplotype
          vectori hap(2);
          hap[0] = m_LD[i][0];
          hap[1] = m_LD[i][1];
          m_haploFreq.addHaplotype( hap, midValues );
        }
      }

      // calculate, right now,  do not tempt to save values
      void apply(Pop& pop)
      {
        if( m_LD.empty())
          return;

        UINT numSP = pop.numSubPop();
        UINT nLD = m_LD.size();

        // remove previous values.
        pop.removeVar(LD_String);
        pop.removeVar(LDPRIME_String);
        pop.removeVar(R2_String);
        pop.removeVar(AvgLDPRIME_String);
        pop.removeVar(AvgR2_String);
        for( UINT sp=0; sp < numSP;  ++sp)
        {
          pop.removeVar( subPopVar_String(sp, LD_String));
          pop.removeVar( subPopVar_String(sp, LDPRIME_String));
          pop.removeVar( subPopVar_String(sp, R2_String));
          pop.removeVar( subPopVar_String(sp, AvgLDPRIME_String));
          pop.removeVar( subPopVar_String(sp, AvgR2_String));
        }

        for(size_t i=0; i < nLD; ++i)
        {
          // specifying alleles
          if( m_LD[i].size() == 4)
          {
            vectori hapLoci(2);
            vectori hapAlleles(2);

            hapLoci[0] = m_LD[i][0];
            hapLoci[1] = m_LD[i][1];
            hapAlleles[0] = m_LD[i][2];
            hapAlleles[1] = m_LD[i][3];

            for( UINT sp=0; sp < numSP;  ++sp)
            {                                     // all necessary numbers have been applyd
              // get haplotype freq
              double P_AB = m_haploFreq.haploFreq(hapLoci, sp)[hapAlleles];
              double P_A = m_alleleFreq.alleleFreq(m_LD[i][2], m_LD[i][0], sp);
              double P_B = m_alleleFreq.alleleFreq(m_LD[i][3], m_LD[i][1],sp);

              // apply LD
              double D = P_AB - P_A * P_B;
              // apply LD'
              double D_max = D > 0 ? std::min(P_A*(1-P_B), (1-P_A)*P_B):std::min(P_A*P_B,(1-P_A)*(1-P_B));
              double D_prime = fcmp_eq(D_max, 0.)?0.:D/D_max;
              double r2 = (fcmp_eq(P_A,0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1))?0.:D*D/P_A/(1-P_A)/P_B/(1-P_B);

              DBG_DO(DBG_STAT, cout << "LD: subpop " << sp << " : P_AB: " << P_AB
                << " P_A: " << P_A << " P_B: " << P_B << " D_max: " << D_max <<
                " LD: " << D << " LD': " << D_prime << " r2: " << r2 << endl);

              // m_D[sp][key] = D;
              // m_Dprime[sp][key] = D_prime;
              // m_r2[sp][key] = r2;
              pop.setDoubleVar( subPopVar_String(sp, LD_String) +
                haploKey( hapLoci) + haploKey(hapAlleles), D);
              pop.setDoubleVar( subPopVar_String(sp, LDPRIME_String) +
                haploKey( hapLoci) + haploKey(hapAlleles), D_prime);
              pop.setDoubleVar( subPopVar_String(sp, R2_String) +
                haploKey( hapLoci) + haploKey(hapAlleles), r2);
            }

            // whole population
            // get haplotype freq
            double P_AB = m_haploFreq.haploFreq(hapLoci)[hapAlleles];
            double P_A = m_alleleFreq.alleleFreq(m_LD[i][2], m_LD[i][0]);
            double P_B = m_alleleFreq.alleleFreq(m_LD[i][3], m_LD[i][1]);

            // apply LD
            double D = P_AB - P_A * P_B;
            // apply LD'
            double D_max = D > 0 ? std::min(P_A*(1-P_B), (1-P_A)*P_B):std::max(P_A*P_B,(1-P_A)*(1-P_B));
            double D_prime = fcmp_eq(D_max, 0)?0:D/D_max;
            double r2 = (fcmp_eq(P_A,0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1))?0:D*D/P_A/(1-P_A)/P_B/(1-P_B);

            DBG_DO(DBG_STAT, cout << "LD: P_AB: " << P_AB
              << " P_A: " << P_A << " P_B: " << P_B << " D_max: " << D_max <<
              " LD: " << D << " LD': " << D_prime << " r2: " << r2 << endl);

            // m_D[numSP][key] = D;
            // m_Dprime[numSP][key] = D_prime;
            // m_r2[numSP][key] = r2;
            pop.setDoubleVar(LD_String + haploKey( hapLoci) + haploKey(hapAlleles), D);
            pop.setDoubleVar(LDPRIME_String + haploKey( hapLoci) + haploKey(hapAlleles), D_prime);
            pop.setDoubleVar(R2_String + haploKey( hapLoci) + haploKey(hapAlleles), r2);
          }
          else                                    // No alleles, using averages
          {
            vectori hapLoci(2);
            vectori hapAlleles(2);

            hapLoci[0] = m_LD[i][0];
            hapLoci[1] = m_LD[i][1];
            string hapLociStr = toStr("[") + toStr(hapLoci[0])
              + toStr("][") + toStr(hapLoci[1]) + toStr("]");

            if( hapLoci[0] == hapLoci[1])         // meaningless?
            {
              for( UINT sp=0; sp < numSP;  ++sp)
              {
                pop.setDoubleVar( subPopVar_String(sp, AvgLD_String) + hapLociStr, 0);
                pop.setDoubleVar( subPopVar_String(sp, AvgLDPRIME_String) + hapLociStr, 0);
                pop.setDoubleVar( subPopVar_String(sp, AvgR2_String) + hapLociStr, 0);
              }
              pop.setDoubleVar( AvgLD_String + hapLociStr, 0);
              pop.setDoubleVar( AvgLDPRIME_String + hapLociStr, 0);
              pop.setDoubleVar( AvgR2_String + hapLociStr, 0);

              continue;
            }

            // find out all alleles
            vectori A_alleles = m_alleleFreq.alleles(hapLoci[0]);
            vectori B_alleles = m_alleleFreq.alleles(hapLoci[1]);

            map<vectori, double>::iterator hapKey;

            for( UINT sp=0; sp < numSP;  ++sp)
            {
              map<vectori, double>& hapFreq = m_haploFreq.haploFreq(hapLoci, sp);

              // all necessary numbers have been applyd
              double D=0.0, D_prime = 0.0, r2 = 0.0;
              for(vectori::iterator A_ale = A_alleles.begin();
                A_ale != A_alleles.end(); ++A_ale)
              {
                for(vectori::iterator B_ale = B_alleles.begin();
                  B_ale != B_alleles.end(); ++B_ale)
                {
                  hapAlleles[0] = *A_ale;
                  hapAlleles[1] = *B_ale;

                  // get haplotype freq
                  double P_AB = (hapKey=hapFreq.find(hapAlleles))==hapFreq.end()?0.:hapKey->second;
                  // alleleFreq(allele, locus, subPop)
                  double P_A = m_alleleFreq.alleleFreq(hapAlleles[0], m_LD[i][0], sp);
                  double P_B = m_alleleFreq.alleleFreq(hapAlleles[1], m_LD[i][1], sp);

                  // apply LD
                  double D_ = P_AB - P_A * P_B;

                  if( m_midValues)
                  {
                    pop.setDoubleVar( subPopVar_String(sp, LD_String) + haploKey(hapLoci) +
                      haploKey(hapAlleles), D_);
                    if( numSP == 1)
                      pop.setDoubleVar( LD_String + haploKey(hapLoci) +
                        haploKey(hapAlleles), D_);
                  }

                  double D_max = D_ > 0 ? std::min(P_A*(1-P_B), (1-P_A)*P_B):std::min(P_A*P_B,(1-P_A)*(1-P_B));

                  double D_prime_ = fcmp_eq(D_max, 0.)?0.:D_/D_max;

                  if( m_midValues)
                  {
                    pop.setDoubleVar( subPopVar_String(sp, LDPRIME_String) + haploKey(hapLoci) +
                      haploKey(hapAlleles), D_prime_);
                    if( numSP == 1)
                      pop.setDoubleVar( LDPRIME_String + haploKey(hapLoci) +
                        haploKey(hapAlleles), D_prime_);
                  }

                  double r2_ = (fcmp_eq(P_A,0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1))?0.:D_*D_/P_A/(1-P_A)/P_B/(1-P_B);

                  if( m_midValues)
                  {
                    pop.setDoubleVar( subPopVar_String(sp, R2_String) + haploKey(hapLoci) +
                      haploKey(hapAlleles), r2_);
                    if( numSP == 1)
                      pop.setDoubleVar( R2_String + haploKey(hapLoci) +
                        haploKey(hapAlleles), r2_);
                  }

                  D += P_A*P_B*fabs(D_);
                  D_prime += P_A*P_B*fabs(D_prime_);
                  r2 += P_A*P_B*r2_;

                  DBG_DO(DBG_STAT, cout << "LD: subpop " << sp
                    << m_LD[i][0] << m_LD[i][1] << " " << hapAlleles[0] << hapAlleles[1]
                    << " : P_AB: " << P_AB
                    << " P_A: " << P_A << " P_B: " << P_B << " D_max: " << D_max
                    << " LD: " << D_ << " LD': " << D_prime_ << " r2: " << r2_ << endl);
                }
              }

              // average take average.

              // m_D[sp][key] = D;
              // m_Dprime[sp][key] = D_prime;
              // m_r2[sp][key] = r2;
              // fixme: room for performance impovement.
              // setStrDictVar instead of setting them one by one.
              // also may need to save result.
              pop.setDoubleVar( subPopVar_String(sp, AvgLD_String) + hapLociStr, D);
              pop.setDoubleVar( subPopVar_String(sp, AvgLDPRIME_String) + hapLociStr, D_prime);
              pop.setDoubleVar( subPopVar_String(sp, AvgR2_String) + hapLociStr, r2);
              if( numSP == 1)
              {
                pop.setDoubleVar( AvgLD_String + hapLociStr, D);
                pop.setDoubleVar( AvgLDPRIME_String + hapLociStr, D_prime);
                pop.setDoubleVar( AvgR2_String + hapLociStr, r2);
              }
            }

            if(numSP > 1 )
            {
              map< vectori, double>& hapFreq = m_haploFreq.haploFreq(hapLoci);

              double D = 0.0, D_prime = 0.0, r2 = 0.0;
              for(vectori::iterator A_ale = A_alleles.begin();
                A_ale != A_alleles.end(); ++A_ale)
              {
                for(vectori::iterator B_ale = B_alleles.begin();
                  B_ale != B_alleles.end(); ++B_ale)
                {
                  hapAlleles[0] = *A_ale;
                  hapAlleles[1] = *B_ale;

                  // whole population
                  // get haplotype freq
                  double P_AB = (hapKey=hapFreq.find(hapAlleles))==hapFreq.end()?0.:hapKey->second;
                  double P_A = m_alleleFreq.alleleFreq(hapAlleles[0], m_LD[i][0]);
                  double P_B = m_alleleFreq.alleleFreq(hapAlleles[1], m_LD[i][1]);

                  // apply LD
                  double D_ = P_AB - P_A * P_B;

                  if( m_midValues)
                    pop.setDoubleVar(LD_String + haploKey(hapLoci) + haploKey(hapAlleles), D_);

                  double D_max = D_ > 0 ? std::min(P_A*(1-P_B), (1-P_A)*P_B):std::min(P_A*P_B,(1-P_A)*(1-P_B));

                  double D_prime_ = fcmp_eq(D_max, 0.)?0.:D_/D_max;

                  if( m_midValues)
                    pop.setDoubleVar(LDPRIME_String + haploKey(hapLoci) + haploKey(hapAlleles), D_prime_);

                  double r2_ = (fcmp_eq(P_A,0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1))?0.:D_*D_/P_A/(1-P_A)/P_B/(1-P_B);

                  if( m_midValues)
                    pop.setDoubleVar(R2_String + haploKey(hapLoci) + haploKey(hapAlleles), r2_);

                  D += P_A * P_B * fabs(D_);
                  D_prime += P_A * P_B * fabs(D_prime_);
                  r2 += P_A * P_B * r2_;

                  DBG_DO(DBG_STAT, cout <<  m_LD[i][0] << m_LD[i][1] << hapAlleles[0] << hapAlleles[1]
                    << "LD: P_AB: " << P_AB  << " P_A: " << P_A << " P_B: " << P_B << " D_max: " << D_max
                    << " LD: " << D_ << " LD': " << D_prime_ << " r2: " << r2_ << endl);

                }                                 // all haplotypes
              }
              pop.setDoubleVar( AvgLD_String + hapLociStr, D);
              pop.setDoubleVar( AvgLDPRIME_String + hapLociStr, D_prime);
              pop.setDoubleVar( AvgR2_String + hapLociStr, r2);
            }                                     // length 2
          }
        }                                         // for all LD
      }

    private:

      /// need to get allele freq
      statAlleleFreq<Pop>& m_alleleFreq;

      /// need to get haplofreq
      statHaploFreq<Pop>& m_haploFreq;

      /// LD
      intMatrix m_LD;

      ///
      bool m_midValues;
  };

  /// CPPONLY
  /// currently there is no need to retrieve calculated value
  template<class Pop>
    class statFst
  {

    private:

#define  Fst_String  "Fst"
#define  Fis_String  "Fis"
#define  Fit_String  "Fit"
#define  AvgFst_String  "AvgFst"
#define  AvgFis_String  "AvgFis"
#define  AvgFit_String  "AvgFit"

    public:
      statFst(statAlleleFreq<Pop>& alleleFreq, statHeteroFreq<Pop>& heteroFreq,
        const vectori& Fst=vectori(), bool midValues=false)
        : m_alleleFreq(alleleFreq), m_heteroFreq(heteroFreq), m_atLoci(Fst)
      {
        for(size_t i=0; i < m_atLoci.size(); ++i)
        {
          // need to get allele frequency at this locus
          m_alleleFreq.addLocus( m_atLoci[i], midValues);

          // need to get heterozygous proportion  at this locus
          m_heteroFreq.addLocus( m_atLoci[i], midValues);
        }
      }

      double Fst()
      {
        return m_avgFst;
      }

      double Fis()
      {
        return m_avgFis;
      }

      double Fit()
      {
        return m_avgFit;
      }

      double Fst(UINT loc)
      {
        return m_Fst[loc];
      }

      double Fis(UINT loc)
      {
        return m_Fis[loc];
      }

      double Fit(UINT loc)
      {
        return m_Fit[loc];
      }

      void apply(Pop& pop)
      {
        if( m_atLoci.empty())
          return;

        m_Fst.clear();
        m_Fit.clear();
        m_Fis.clear();

        // dicitonary to save values.
        UINT numSP = pop.numSubPop();
        ULONG popSize = pop.popSize();

        // do not save these values now
        double aa = 0., bb = 0., cc = 0.;

        // vector to store p[i]
        vectorf p_i = vectorf(numSP);

        // calculate Fst for each locus
        for(size_t st = 0; st < m_atLoci.size(); ++st)
        {
          int loc = m_atLoci[st];

          DBG_ASSERT( static_cast<size_t>(loc) < pop.totNumLoci(), IndexError,
            "Index out of range of 0 ~ " + toStr(pop.totNumLoci()-1));

          // get all available alleles
          vectori alleles = m_alleleFreq.alleles(loc);

          DBG_DO(DBG_STAT, cout << "Using alleles " << alleles << endl);

          // n_bar
          double r = numSP;
          double n = popSize;
          double n_bar = n/r;
          vectorlu n_i = pop.subPopSizes();

          // n_c
          double n_c = n;
          for( int i = 0; i < r; ++i)
            n_c  -= n_i[i]*n_i[i]/n;
          n_c /= (r - 1);

          double a=0.0, b=0.0, c=0.0;

          for( vectori::iterator ale = alleles.begin(); ale != alleles.end(); ++ale)
          {
            // p_i
            for(int sp=0; sp < r; ++sp)
              p_i[sp] = m_alleleFreq.alleleFreq(*ale, loc, sp);

            // p_bar
            double p_bar = 0;
            for(int sp = 0; sp < r; ++sp)
              p_bar += n_i[sp] * p_i[sp];
            p_bar /= n;

            // s^2
            double s_2 = 0;
            for(int sp=0; sp < r; ++sp)
              s_2 += n_i[sp]*(p_i[sp] - p_bar)*(p_i[sp]-p_bar);
            s_2 /= (r-1)*n_bar;

            // h_bar
            double h_bar = 0;
            for(int sp=0; sp < r; ++sp)
              h_bar += m_heteroFreq.heteroFreq(*ale, loc, sp) * n_i[sp];
            h_bar /= n;

            // a, b, c
            a += n_bar/n_c *( s_2 - ( p_bar*(1-p_bar) - (r-1.)/r*s_2 - h_bar/4.)/(n_bar -1. ) );
            b += n_bar / (n_bar -1 )* ( p_bar *(1-p_bar) - (r-1)/r * s_2 - (2 * n_bar -1 )/(4.*n_bar)* h_bar );
            c += h_bar /2.;

            DBG_DO(DBG_STAT, cout << "allele " << *ale << "\tn_c: " << n_c
              << "\tp_i: " << p_i << "\tp_bar: " << p_bar << "\ts^2: " << s_2 << "\th_bar:"
              << h_bar << "\ta: " << a << "\tb: " << b << "\tc: " << c << endl);
          }                                       // each allele

          DBG_DO(DBG_STAT, cout << "Fst= " << a/(a+b+c) << endl);

          if( static_cast<size_t>(loc) >= m_Fst.size())
          {
            m_Fst.resize(loc+1,0.);
            m_Fit.resize(loc+1,0.);
            m_Fis.resize(loc+1,0.);
          }
          m_Fst[loc] = fcmp_eq(a+b+c,0.)?0.:(a / ( a+b+c));
          m_Fit[loc] = fcmp_eq(a+b+c,0.)?1.:(1 - c / (a+b+c));
          m_Fis[loc] = fcmp_eq(b+c,0.)?1.:(1 - c /( b+c));

          aa += a;
          bb += b;
          cc += c;
        }
        m_avgFst = fcmp_eq(aa+bb+cc,0.)?0:(aa / ( aa+bb+cc));
        m_avgFit = fcmp_eq(aa+bb+cc,0.)?1.:(1 - cc / (aa+bb+cc));
        m_avgFis = fcmp_eq(aa+bb+cc,0)?1.:(1 - cc /( bb+cc));

        // post results
        pop.setDoubleNumArrayVar(Fst_String, m_Fst.size(), &m_Fst[0]);
        pop.setDoubleNumArrayVar(Fit_String, m_Fit.size(), &m_Fit[0]);
        pop.setDoubleNumArrayVar(Fis_String, m_Fis.size(), &m_Fis[0]);
        pop.setDoubleVar(AvgFst_String, m_avgFst);
        pop.setDoubleVar(AvgFit_String, m_avgFit);
        pop.setDoubleVar(AvgFis_String, m_avgFis);
      }

    private:

      statAlleleFreq<Pop>& m_alleleFreq;
      statHeteroFreq<Pop>& m_heteroFreq;

      /// Fst
      vectori m_atLoci;

      /// result
      vectorf m_Fst, m_Fit, m_Fis;

      double m_avgFst, m_avgFit, m_avgFis;
  };

#define REL_Queller             1
#define REL_Lynch               2
#define REL_IR                  3
#define REL_D2                  4
#define REL_Rel                 5

  /// the relatedness measure between two individuals/families
  /// using Queller and Goodnight or Lynch's method.
  /// or internal relatedness values
  template<class Pop>
    class statRelatedness
  {

    public:

#define Rel_Queller_String "relQueller"
#define Rel_Lynch_String   "relLynch"
#define Rel_IR_String      "relIR"
#define Rel_D2_String      "relD2"
#define Rel_Rel_String     "relRel"

    private:

      typedef std::pair<double, double> fraction;

    public:

      /// \brief calculate relatedness measures between elements in groups
      /**
      \param groups can be [ [1,2,3],[4,5,6],[7,8,9]] as three groups of
      individuals; or [ 1 3 4] as three subpopulations. To specify between
      individual relatedness, use [[1],[2],[3]] (the first form). If this
      parameter is ignored, this operator calculate relatedness between
      all subpopulations.

      \param method can be REL_Queller, REL_Lynch, REL_IR, REL_D2
      or REL_Rel. Please refer to the manual for details.
      */
      statRelatedness(statAlleleFreq<Pop>& alleleFreq, const intMatrix& groups=intMatrix(),
        bool useSubPop=false, const vectori& loci=vectori(), vectori method=vectori(),
        int minScored=10, bool midValues=false):
      m_alleleFreq(alleleFreq), m_groups(groups), m_useSubPop(useSubPop),
        m_atLoci(loci), m_method(method), m_minScored(minScored)
      {
        if( m_groups.empty() || m_groups[0].empty() )
          return;

        DBG_FAILIF(  method.empty(), ValueError, "Please specify relatedness method");

        DBG_FAILIF( m_atLoci.empty(), ValueError, "Please specify parameter relLoci.");

        for(size_t i=0; i<m_atLoci.size(); ++i)
          m_alleleFreq.addLocus(m_atLoci[i], midValues);
      }

      // relatedness between individuals
      fraction relQueller(typename Pop::IndType ind1,
        typename Pop::IndType ind2)
      {

        matrix& af = m_alleleFreq.alleleFreqAll();
        fraction res(0.,0.);
        int numScored = 0;
        for(vectori::iterator locus=m_atLoci.begin();
          locus != m_atLoci.end(); ++locus)
        {
          Allele a = ind1.allele(*locus, 0);
          Allele b = ind1.allele(*locus, 1);
          Allele c = ind2.allele(*locus, 0);
          Allele d = ind2.allele(*locus, 1);

          if( a==0 || b==0 || c==0 || d==0 ) continue;

          double s = af[*locus][a] + af[*locus][b] + af[*locus][c] + af[*locus][d];
          double r = (a==c) + (b==c) + (a==d) + (b==d);
          double h = 2 + (a==b) + (c==d);

          res.first += (r - s);
          res.second += (h - s);
          numScored ++;
        }
        // chromosome
        // cout << numScored <<  " " << res.value() << endl;
        if( res.second == 0 || numScored <= m_minScored)
          cout << "Warning: Not enough alleles available to calculate relatedness. Num scored = "
            << numScored << endl;

        return res;
      }

      fraction relLynch(typename Pop::IndType ind1,
        typename Pop::IndType ind2)
      {
        matrix& af = m_alleleFreq.alleleFreqAll();

        fraction res(0.,0.);
        int numScored = 0;

        double rel_xy = 0.0;
        double rel_yx = 0.0;
        double weight_xy = 0.0;
        double weight_yx = 0.0;
        for(vectori::iterator locus=m_atLoci.begin();
          locus != m_atLoci.end(); ++locus)
        {
          Allele a = ind1.allele(*locus, 0);
          Allele b = ind1.allele(*locus, 1);
          Allele c = ind2.allele(*locus, 0);
          Allele d = ind2.allele(*locus, 1);

          if( a==0 || b==0 || c==0 || d==0 ) continue;

          double pa = af[*locus][a];
          double pb = af[*locus][b];
          double pc = af[*locus][c];
          double pd = af[*locus][d];

          if( pa < 1e-8 || pb < 1e-8 || pc < 1e-8 || pd < 1e-8) continue;

          double r_xy = ( pa*((b==c)+(b==d)) + pb*((a==c)+(a==d)) - 4*pa*pb ) /
            ( (1+(a==b))*(pa+pb) - 4*pa*pb);
          double r_yx = ( pc*((d==a)+(d==b)) + pd*((c==a)+(c==b)) - 4*pc*pd ) /
            ( (1+(c==d))*(pc+pd) - 4*pc*pd);
          double w_xy = ( (1+(a==b))*(pa+pb) - 4*pa*pb) / (2*pa*pb);
          double w_yx = ( (1+(c==d))*(pc+pd) - 4*pc*pd) / (2*pc*pd);

          rel_xy += w_xy*r_xy;
          rel_yx += w_yx*r_yx;
          weight_xy += w_xy;
          weight_yx += w_yx;

          numScored ++;
        }                                         // all loci
        DBG_FAILIF( numScored <= m_minScored || weight_xy == 0.0 || weight_yx == 0.0,
          ValueError, "Not enough allels to calculated relatedness.");
        res.first = (rel_xy/weight_xy + rel_yx/weight_yx)/2.;
        res.second = 1.;
        return res;
      }

      // IR measure for individual ind at specified locus
      fraction relIR(typename Pop::IndType ind1, int locus)
      {
        matrix& af = m_alleleFreq.alleleFreqAll();
        fraction res(0.,0.);
        Allele a = ind1.allele(locus, 0);
        Allele b = ind1.allele(locus, 1);
        double pa = af[locus][a];
        double pb = af[locus][b];

        res.first = 2*(a==b) - pa - pb;
        res.second = 2 - pa - pb;

        if( fcmp_eq( res.second, 0) )
          cout << "Warning: IR value: pa = pb =1. NA will be returned." << endl;
        return res;
      }

      // D2 measure for individual ind at specified locus
      fraction relD2(typename Pop::IndType ind1, int locus)
      {
        fraction res(0.,1.);
        Allele a = ind1.allele(locus, 0);
        Allele b = ind1.allele(locus, 1);

        UINT mx = m_alleleFreq.numOfAlleles()[locus];
        res.first = ((a - b)/(mx-2.0))*((a - b)/(mx-2.0));

        return res;
      }

      // REL measure for individual ind at specified locus
      fraction relRel(typename Pop::IndType ind1,
        typename Pop::IndType ind2,  int locus)
      {
        matrix& af = m_alleleFreq.alleleFreqAll();

        fraction res(0.,0.);
        Allele a = ind1.allele(locus, 0);
        Allele b = ind1.allele(locus, 1);
        Allele c = ind2.allele(locus, 0);
        Allele d = ind2.allele(locus, 1);

        if( a==0 || b==0 || c==0 || d==0 )
          return res;

        double s = af[locus][a] + af[locus][b] + af[locus][c] + af[locus][d];
        double r = (a==c) + (b==c) + (a==d) + (b==d);
        double h = 2 + (a==b) + (c==d);

        res.first += (r - s);
        res.second += (h - s);

        return res;
      }

      // between group i and j if method=REL_Queller and REL_Lynch
      /// for group i and locus j otherwise
      double groupRelatedness(Pop& pop, int i, int j, int method)
      {
        fraction res(0., 0.);

        if( m_useSubPop)
        {
          UINT sp1, sp2;
          if( m_groups[0].empty())
          {
            sp1 = i;
            sp2 = j;
          }
          else
          {
            sp1 = m_groups[0][i];
            sp2 = m_groups[0][j];
          }

          switch(method)
          {
            case REL_Queller:
              // from subpop i and j
              for(typename Pop::IndIterator ind1 = pop.indBegin(sp1);
                ind1 != pop.indEnd(sp1); ++ind1)
              {
                for(typename Pop::IndIterator ind2 = pop.indBegin(sp2);
                  ind2 != pop.indEnd(sp2); ++ind2)
                {
                  fraction tmp= relQueller(*ind1, *ind2);
                  res.first += tmp.first/tmp.second;
                  res.second += 1.;
                }
              }
              return res.first/res.second;
            case REL_Lynch:
              // from subpop i and j
              for(typename Pop::IndIterator ind1 = pop.indBegin(sp1);
                ind1 != pop.indEnd(sp1); ++ind1)
              {
                for(typename Pop::IndIterator ind2 = pop.indBegin(sp2);
                  ind2 != pop.indEnd(sp2); ++ind2)
                {
                  fraction tmp= relLynch(*ind1, *ind2);
                  res.first += tmp.first/tmp.second;
                  res.second += 1.;
                }
              }                                   // lynch
              return res.first/res.second;
            case REL_IR:
              for(typename Pop::IndIterator ind1 = pop.indBegin(sp1);
                ind1 != pop.indEnd(sp1); ++ind1)
              {
                fraction tmp= relIR(*ind1, j);
                res.first += tmp.first;
                res.second += tmp.second;
              }                                   // lynch
              return res.first/res.second;
            case REL_D2:
              for(typename Pop::IndIterator ind1 = pop.indBegin(sp1);
                ind1 != pop.indEnd(sp1); ++ind1)
              {
                fraction tmp= relD2(*ind1, j);
                res.first += tmp.first/tmp.second;
                res.second += 1.;
              }                                   // lynch
              return res.first/res.second;
            case REL_Rel:
              for(typename Pop::IndIterator ind1 = pop.indBegin(sp1);
                ind1 != pop.indEnd(sp1); ++ind1)
              {
                for(typename Pop::IndIterator ind2 = ind1 + 1;
                  ind2 != pop.indEnd(sp1); ++ind2)
                {
                  fraction tmp= relRel(*ind1, *ind2, j);
                  res.first += tmp.first;
                  res.second += tmp.second;
                }
              }                                   // lynch
              return res.first/res.second;
          }                                       // switch
        }                                         // m_useSubPop
        else
        {
          switch(method)
          {
            case REL_Queller:
              // from specified group i and j
              for(vectori::iterator ind1 = m_groups[i].begin();
                ind1 != m_groups[i].end(); ++ind1)
              {
                for(vectori::iterator ind2 = m_groups[j].begin();
                  ind2 != m_groups[j].end(); ++ind2)
                {
                  fraction tmp = relQueller(pop.individual(*ind1), pop.individual(*ind2));
                  res.first += tmp.first/tmp.second;
                  res.second += 1.;
                }
              }
              return res.first/res.second;
            case REL_Lynch:
              for(vectori::iterator ind1 = m_groups[i].begin();
                ind1 != m_groups[i].end(); ++ind1)
              {
                for(vectori::iterator ind2 = m_groups[j].begin();
                  ind2 != m_groups[j].end(); ++ind2)
                {
                  fraction tmp = relLynch(pop.individual(*ind1), pop.individual(*ind2));
                  res.first += tmp.first/tmp.second;
                  res.second += 1.;
                }
              }
              return res.first/res.second;
            case REL_IR:
              for(vectori::iterator ind1 = m_groups[i].begin();
                ind1 != m_groups[i].end(); ++ind1)
              {
                fraction tmp = relIR(pop.individual(*ind1), j);
                res.first += tmp.first;
                res.second += tmp.second;
              }
              return res.first/res.second;
            case REL_D2:
              for(vectori::iterator ind1 = m_groups[i].begin();
                ind1 != m_groups[i].end(); ++ind1)
              {
                fraction tmp = relD2(pop.individual(*ind1), j);
                res.first += tmp.first/tmp.second;
                res.second += 1.;
              }
              return res.first/res.second;
            case REL_Rel:
              for(vectori::iterator ind1 = m_groups[i].begin();
                ind1 != m_groups[i].end(); ++ind1)
              {
                for(vectori::iterator ind2 = ind1 + 1;
                  ind2 != m_groups[i].end(); ++ind2)
                {
                  fraction tmp = relRel(pop.individual(*ind1),
                    pop.individual(*ind2), j);
                  res.first += tmp.first;
                  res.second += tmp.second;
                }
              }
              return res.first/res.second;
          }                                       // switch
        }                                         // m_useSubPop
        // will never reach here.
        return 0;
      }

      void apply(Pop& pop)
      {
        if(m_groups.empty())
          return;

        // calculate relatendness values between two groups
        UINT nGroups;
        if(m_useSubPop)
        {
          if( m_groups[0].empty() )
            nGroups = pop.numSubPop();
          else
            nGroups = m_groups[0].size();
        }
        else
          nGroups = m_groups.size();

        UINT nLoci = m_atLoci.size();

        for(size_t m=0; m < m_method.size(); ++m)
        {
          switch(m_method[m])
          {
            case REL_Queller:
              m_relQueller.resize(nGroups);
              for(UINT i=0; i<nGroups; ++i)
                m_relQueller[i].resize(nGroups);

              for(UINT i=0; i<nGroups; ++i)
                for(UINT j=i+1; j<nGroups; ++j)
              {
                m_relQueller[i][j] = groupRelatedness(pop, i, j, m_method[m]);
                m_relQueller[j][i] = m_relQueller[i][j];
              }

              for(UINT i=0; i<nGroups; ++i)
                pop.setDoubleNumArrayVar(Rel_Queller_String + toStr("[")
                  + toStr(i) + "]", nGroups, &(m_relQueller[i][0]));
              break;
            case REL_Lynch:
              m_relLynch.resize(nGroups);
              for(UINT i=0; i<nGroups; ++i)
                m_relLynch[i].resize(nGroups);

              for(UINT i=0; i<nGroups; ++i)
                for(UINT j=i+1; j<nGroups; ++j)
              {
                m_relLynch[i][j] = groupRelatedness(pop, i, j, m_method[m]);
                m_relLynch[j][i] = m_relLynch[i][j];
              }

              for(UINT i=0; i<nGroups; ++i)
                pop.setDoubleNumArrayVar(Rel_Lynch_String + toStr("[")
                  + toStr(i) + "]", nGroups, &(m_relLynch[i][0]));
              break;
            case REL_IR:
              m_relIR.resize(nGroups);
              for(UINT i=0; i<nGroups; ++i)
                m_relIR[i].resize(nLoci);

              for(UINT i=0; i<nGroups; ++i)
                for(UINT j=0; j<nLoci; ++j)
                  m_relIR[i][j] = groupRelatedness(pop, i, m_atLoci[j], m_method[m]);

              for(UINT i=0; i<nGroups; ++i)
                pop.setDoubleNumArrayVar(Rel_IR_String + toStr("[")
                  + toStr(i) + "]", nLoci, &(m_relIR[i][0]));
              break;
            case REL_D2:
              m_relD2.resize(nGroups);
              for(UINT i=0; i<nGroups; ++i)
                m_relD2[i].resize(nLoci);

              for(UINT i=0; i<nGroups; ++i)
                for(UINT j=0; j<nLoci; ++j)
                  m_relD2[i][j] = groupRelatedness(pop, i, m_atLoci[j], m_method[m]);

              for(UINT i=0; i<nGroups; ++i)
                pop.setDoubleNumArrayVar(Rel_D2_String + toStr("[")
                  + toStr(i) + "]", nLoci, &(m_relD2[i][0]));
              break;
            case REL_Rel:
              m_relRel.resize(nGroups);
              for(UINT i=0; i<nGroups; ++i)
                m_relRel[i].resize(nLoci);

              for(UINT i=0; i<nGroups; ++i)
                for(UINT j=0; j<nLoci; ++j)
                  m_relRel[i][j] = groupRelatedness(pop, i, m_atLoci[j], m_method[m]);

              for(UINT i=0; i<nGroups; ++i)
                pop.setDoubleNumArrayVar(Rel_Rel_String + toStr("[")
                  + toStr(i) + "]", nLoci, &(m_relRel[i][0]));
              break;
          }
        }                                         // for all method
      }

    private:

      /// need to get allele freq
      statAlleleFreq<Pop>& m_alleleFreq;

      /// method to use
      intMatrix m_groups;

      ///
      bool m_useSubPop;

      /// loci used
      vectori m_atLoci;

      /// method
      vectori m_method;

      /// control number of scored loci
      int m_minScored;

      // save result
      matrix m_relQueller, m_relLynch, m_relIR, m_relD2, m_relRel;
  };

  template<class Pop>
    class BasicStat: public Stator<Pop>
  {
    public:

      /// create an basicStat
      /**
      \param popSize whether or not calculate population sizes. will set numSubPop, subPopSize,
        popSize, subPop[sp]['popSize']

      \param numOfMale whether or not count number of male and female, will set numOfMale and numOfFemale for
        all population/subpopulations

      \param numOfAffected whether or not count number of affected individuals. Will set numOfAffected, and numOfUnaffected.

      \param numOfAlleles an array of loci at which number of alleles will be counted (0 is excluded). Note that
        number of alleles will be automatically set if alleleFreq is counted.

      \param alleleFreq an array of loci at which all alleles will be counted.

      \param genoFreq an array of loci at which all genotype will be counted

      each item is the locus index followed by allele pairs.

      \param heteroFreq an array of loci at which the observed proportion of individuausl heterozygous
      will be applyd for each allele. Expected heterozygosity will also be calculuate and put
      in heteroFreq[locus][0] (since allele 0 is not used.)

      \param homoFreq an array of loci at which homozygosity number and frequcies will be calculated

      \param expHetero an array of loci at which expected heterozygosity will be calculated.

      \param haploFreq a matrix of haplotypes (allele sequence) to count

      format:  haploFreq = [ [ 0,1,2 ], [1,2] ]

      All haplotypes on loci 012, 12 will be counted. If only one haplotype is
      specified, the outer [] can be ommited. I.e., haploFreq=[0,1] is acceptable.

      \param LD apply LD, LD' and r2, given
      LD=[ [locus1 locus2], [ locus1 locus2 allele1 allele2], ,...]
      If two numbers are given, D, D' and r2 overall possible allele pairs will
      be calculated and saved as AvgLD etc. If four numbers are given, D, D' and r2
      using specified alleles are provided. If only one item is specified, the
      outer [] can be ignored. I.e., LD=[locus1 locus2] is acceptable.

      \param Fst calculate Fst. Fis and Fit will be given as a side product.

      format
      Fst = [ 0, 1 ] Fst calculated at locus 0 and 1. Since any allele can be used
      to calculate Fst, Fst[0] will be the average over all alleles as suggested
      by Weir and Cockerham.

      \param relGroups  calculated pairwise relatedness between groups. relGroups
      can be in the form of either [[1,2,3],[4,5],[7,8]] (gourps of individuals)
      or [1,3,4] (use subpopulations).

      \param relMethod method used to calculate relatedness. Can be
      either REL_Queller or REL_Lynch.

      \param relLoci loci on which relatedness calues are calculated.

      \param hasPhase if a/b and b/a are the same genotype. default is false.

      \param midValues whether or not post intermediate results. Default to false.
      For example, Fst will need to calculate allele frequency, if midValues is set
      to true, allele frequencies will be posted as well. This will help debuging and
      sometimes derived statistics.

      \param others there is NO output for this operator. other parameters
      please see help(baseOperator.__init__)

      \note previous provide 'search for all allele/genotype' option but this
      has proven to be troublesome. In this version, everything should be
      explicitly specified.
      **/
      BasicStat(
        bool popSize=false,
        bool numOfMale=false,
        bool numOfAffected=false,
        vectori numOfAlleles=vectori(),
        vectori alleleFreq=vectori(),
        vectori heteroFreq=vectori(),
        vectori expHetero=vectori(),
        vectori homoFreq=vectori(),
        vectori genoFreq=vectori(),
        intMatrix haploFreq=intMatrix(),
        intMatrix LD=intMatrix(),
        vectori Fst=vectori(),
        intMatrix relGroups=intMatrix(),
        vectori relLoci=vectori(),
        bool relBySubPop=false,                   // internal use
        vectori relMethod=vectori(),
        int relMinScored=10,                      // minimal number of loci required.
        bool hasPhase=false,
        bool midValues=false,
      // regular parameters
        string output="", string outputExpr="",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL,  int grp=GRP_ALL, string sep="\t")
        :Stator<Pop>("", outputExpr, stage, begin, end, step, at, rep, grp, sep),
      // the order of initialization is meaningful since they may depend on each other
        m_popSize(popSize),
        m_numOfMale(numOfMale),
        m_numOfAffected(numOfAffected),
        m_alleleFreq(alleleFreq),
        m_numOfAlleles(m_alleleFreq, numOfAlleles),
        m_heteroFreq(heteroFreq, homoFreq),
        m_expHetero(m_alleleFreq, expHetero),
        m_genoFreq(genoFreq, hasPhase),
        m_haploFreq(haploFreq),
        m_LD(m_alleleFreq, m_haploFreq, LD, midValues),
        m_Fst(m_alleleFreq, m_heteroFreq, Fst, midValues),
        m_relatedness(m_alleleFreq, relGroups, relBySubPop, relLoci,
        relMethod, relMinScored, midValues)
      {
      }

      ~BasicStat()
      {
      }

      /// this function is very important
      virtual Operator<Pop>* clone() const
      {
        return new BasicStat<Pop>(*this);
      }

      /// count various statistics.
      /// use m_alleles etc to save (potentially) time to
      /// resize all these variables.
      virtual bool apply(Pop & pop)
      {
        m_popSize.apply(pop);
        m_numOfMale.apply(pop);
        m_numOfAffected.apply(pop);
        m_alleleFreq.apply(pop);
        m_heteroFreq.apply(pop);
        m_expHetero.apply(pop);
        m_genoFreq.apply(pop);
        m_haploFreq.apply(pop);
        m_LD.apply(pop);
        m_Fst.apply(pop);
        m_relatedness.apply(pop);
        return true;
      }

      virtual string __repr__()
      {
        return "<simuPOP::basic statistics>";
      }

    private:
      statPopSize<Pop> m_popSize;
      statNumOfMale<Pop> m_numOfMale;
      statNumOfAffected<Pop> m_numOfAffected;
      statAlleleFreq<Pop> m_alleleFreq;
      statNumOfAlleles<Pop> m_numOfAlleles;
      statHeteroFreq<Pop> m_heteroFreq;
      statExpHetero<Pop> m_expHetero;
      statGenoFreq<Pop> m_genoFreq;
      statHaploFreq<Pop> m_haploFreq;
      statLD<Pop> m_LD;
      statFst<Pop> m_Fst;
      statRelatedness<Pop> m_relatedness;
  };
}
#endif
