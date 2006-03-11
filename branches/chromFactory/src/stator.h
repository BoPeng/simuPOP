/**************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu
 *                                                                         *
 *   $LastChangedDate$
 *   $Rev$                                                      *
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
\brief head file of class stator:public Operator
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

  class stator: public Operator
  {
    public:
      /// constructor. default to be always active.
      /// default to have NO output (shared variables will be set.)
      /// phase: if we treat Aa!=aA, default is false, i.e, Aa=aA.
      stator(string output="", string outputExpr="",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL):
      Operator(output, outputExpr, stage, begin, end, step, at, rep, grp)
      {
      };

      /// destructor
      virtual ~stator()
      {
      }

      /// this function is very important
      virtual Operator* clone() const
      {
        return new stator(*this);
      }
  };

  /// evaluate an expression.
  ///

  class pyEval: public stator
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
      pyEval(const string& expr="", const string& stmts="", const string& preStmts="",
        const string& postStmts="", bool exposePop=false, const string& name="",
        string output=">", string outputExpr="",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL,  int grp=GRP_ALL)
        :stator(output, outputExpr, stage, begin, end, step, at, rep, grp),
        m_expr(expr, stmts), m_postExpr("", postStmts), m_exposePop(exposePop), m_name(name)
      {
        if( preStmts != "")
          Expression("", preStmts).evaluate();
      }

      ~pyEval()
      {
        m_postExpr.evaluate();
      }

      /// this function is very important
      virtual Operator* clone() const
      {
        return new pyEval(*this);
      }

      // check all alleles in vector allele if they are fixed.
      virtual bool apply(population& pop);

      virtual string __repr__()
      {
        return "<simuPOP::pyEval " + m_name + ">";
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

  class pyExec: public pyEval
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
      pyExec( const string& stmts="", const string & preStmts="", const string& postStmts="",
        bool exposePop=false, const string& name="",
        string output=">", string outputExpr="",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL,  int grp=GRP_ALL)
        :pyEval("", stmts, preStmts, postStmts, exposePop, name, "", "",
        stage, begin, end, step, at, rep, grp)
      {
      }

      /// this function is very important
      virtual Operator* clone() const
      {
        return new pyExec(*this);
      }

      virtual string __repr__()
      {
        return "<simuPOP::pyExec " + this->name() + ">";
      }
  };

  ///
  /// The following classes apply various statistics
  /// and stat class will provide an opearator interface
  /// to all of them.

  ///
  /// each class defines how to apply the statistics and
  /// provide interface to allow others to retrieve the value.
  /// NOTE: the values are population dependent so these
  /// values are only meaningful within the same stat.apply
  /// call;
  ///
  /// The design separate the calculation of statistics
  /// as much as possible and allow easier debug and writting
  /// new statistics.

  /// CPPONLY
  /// return {'a-b-b'} for a b c
  string haploKey(const vectori& seq);

  /// CPPONLY
  /// post population sizes etc.

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

      void apply(population& pop)
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

        pop.setIntVectorVar(subPopSize_String, spSize );

        for( size_t sp=0; sp < numSP; ++sp)
          pop.setIntVar(subPopVar_String(sp, popSize_String), spSize[sp]);
      }

    private:
      bool m_isActive;
  };

  /// CPPONLY

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

      void apply(population& pop);

    private:
      /// whether or not to apply number of male/female
      vectorlu m_numOfMale, m_numOfFemale;
  };

  /// CPPONLY

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

      void apply(population& pop);

    private:
      /// record the result.
      vectorlu m_numOfAffected, m_numOfUnaffected;
  };

  /// CPPONLY

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

        for(size_t j=0; j < m_alleleNum.back()[loc].size(); ++j)
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

        // whether or not count 0 allele?
        // consider NA as an allele is reasonable.
        for(size_t j=0; j < m_alleleNum[subPop][loc].size(); ++j)
          if( m_alleleNum[subPop][loc][j] != 0)
            al.push_back(j);

#ifndef BINARYALLELE
        DBG_WARNING(m_alleleNum[subPop][loc][0] != 0,
          "Having zero (NA) allele, counted as one allele.");
#endif

        DBG_ASSERT( al.size() == static_cast<UINT>(numOfAlleles(subPop)[loc]),
          SystemError, "Number of alleles at locus " + toStr(loc)
          + " at subpop " + toStr(subPop) + " does not match. Observed "
          + toStr( al.size()) + " previous count: "
          + toStr( numOfAlleles(subPop)[loc] ));

        return al;
      }

      void apply(population& pop);

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

  class statNumOfAlleles
  {
    public:
      statNumOfAlleles(statAlleleFreq& calc, const vectori& atLoci = vectori())
        :m_calc(calc)
      {
        for(vectori::const_iterator it = atLoci.begin(); it != atLoci.end(); ++it)
          m_calc.addLocus(*it, true);
      }

      // do nothing. m_calc.spply will be called by stat.
      void apply(population& pop)
      {
      }

    private:

      /// a reference to an existing allelefreq calculator
      statAlleleFreq& m_calc;
  };

  /// CPPONLY

  class statHeteroFreq
  {
    private:

#define HeteroNum_String        "heteroNum"
#define HeteroFreq_String       "heteroFreq"
#define AllHeteroNum_String     "HeteroNum"
#define AllHeteroFreq_String    "HeteroFreq"
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

      int heteroNum(UINT allele, int loc)
      {
        UINT idx=locusIdx(loc);

        vectori& hn =  m_heteroNum[resIdx(idx)];
        return allele < hn.size() ? hn[allele] : 0;
      }

      double heteroFreq(UINT allele, int loc)
      {
        UINT idx=locusIdx(loc);
        vectorf& hf =  m_heteroFreq[resIdx(idx)];
        return allele < hf.size() ? hf[allele] : 0.;
      }

      int heteroNum(UINT allele, int loc, UINT subPop)
      {
        UINT idx=locusIdx(loc);
        vectori& hn =  m_heteroNum[resIdx(idx,subPop)];
        return allele < hn.size() ? hn[allele] : 0;
      }

      double heteroFreq(UINT allele, int loc, UINT subPop)
      {
        UINT idx=locusIdx(loc);
        vectorf& hf =  m_heteroFreq[resIdx(idx,subPop)];
        return allele < hf.size() ? hf[allele] : 0;
      }

      void apply(population& pop);

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

  class statExpHetero
  {
    private:

#define ExpHetero_String "expHetero"

    public:
      statExpHetero(statAlleleFreq& alleleFreq, const vectori& expHetero=vectori())
        : m_alleleFreq(alleleFreq), m_atLoci(expHetero), m_expHetero(0)
      {
        // add expected hetero to m_alleleFreq
        for(size_t i=0; i < expHetero.size(); ++i)
          m_alleleFreq.addLocus( expHetero[i]);
      }

      void apply(population& pop);

    private:

      /// need this to apply alleleFreq
      statAlleleFreq& m_alleleFreq;

      /// heteroFreq
      vectori m_atLoci;

      /// expected heterozygosity
      matrix m_expHetero;
  };

  /// CPPONLY
  /// currently there is no need to expose the result.
  /// may add that later.

  class statGenoFreq
  {
    private:

#define  GenotypeNum_String   "genoNum"
#define  GenotypeFreq_String  "genoFreq"

    public:
      statGenoFreq(const vectori& genoFreq = vectori(), bool phase=true )
        : m_atLoci(genoFreq), m_phase(phase)
      {
      }

      void apply(population& pop);

    private:
      /// which genotypes
      vectori m_atLoci;

      /// phase
      bool m_phase;
  };

  /// CPPONLY

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

      void apply(population& pop);

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

      statLD(statAlleleFreq& alleleFreq, statHaploFreq& haploFreq,
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
      void apply(population& pop);

    private:

      /// need to get allele freq
      statAlleleFreq& m_alleleFreq;

      /// need to get haplofreq
      statHaploFreq& m_haploFreq;

      /// LD
      intMatrix m_LD;

      ///
      bool m_midValues;
  };

  /// CPPONLY
  /// currently there is no need to retrieve calculated value

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
      statFst(statAlleleFreq& alleleFreq, statHeteroFreq& heteroFreq,
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

      void apply(population& pop);

    private:

      statAlleleFreq& m_alleleFreq;
      statHeteroFreq& m_heteroFreq;

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

  class statRelatedness
  {

    public:

#define Rel_Queller_String "relQueller"
#define Rel_Lynch_String   "relLynch"
#define Rel_IR_String      "relIR"
#define Rel_D2_String      "relD2"
#define Rel_Rel_String     "relRel"

    public:

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
      statRelatedness(statAlleleFreq& alleleFreq, const intMatrix& groups=intMatrix(),
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
      fraction relQueller(individual ind1,
        individual ind2);

      fraction relLynch(individual ind1,
        individual ind2);

      // IR measure for individual ind at specified locus
      fraction relIR(individual ind1, int locus);

      // D2 measure for individual ind at specified locus
      fraction relD2(individual ind1, int locus);

      // REL measure for individual ind at specified locus
      fraction relRel(individual ind1,
        individual ind2,  int locus);

      // between group i and j if method=REL_Queller and REL_Lynch
      /// for group i and locus j otherwise
      double groupRelatedness(population& pop, int i, int j, int method);

      void apply(population& pop);

    private:

      /// need to get allele freq
      statAlleleFreq& m_alleleFreq;

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

  class stat: public stator
  {
    public:

      /// create an stat
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
      stat(
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
        int rep=REP_ALL,  int grp=GRP_ALL)
        :stator("", outputExpr, stage, begin, end, step, at, rep, grp),
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

      ~stat()
      {
      }

      /// this function is very important
      virtual Operator* clone() const
      {
        return new stat(*this);
      }

      /// count various statistics.
      /// use m_alleles etc to save (potentially) time to
      /// resize all these variables.
      virtual bool apply(population& pop)
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
        return "<simuPOP::statistics>";
      }

    private:
      statPopSize m_popSize;
      statNumOfMale m_numOfMale;
      statNumOfAffected m_numOfAffected;
      statAlleleFreq m_alleleFreq;
      statNumOfAlleles m_numOfAlleles;
      statHeteroFreq m_heteroFreq;
      statExpHetero m_expHetero;
      statGenoFreq m_genoFreq;
      statHaploFreq m_haploFreq;
      statLD m_LD;
      statFst m_Fst;
      statRelatedness m_relatedness;
  };
}
#endif
