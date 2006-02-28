/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu
 *                                                                         *
 *   $LastChangedDate$
 *   $Rev$                                                       *
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

#ifndef _MATING_H
#define _MATING_H
/**
\file
\brief head file of class mating and its subclasses
*/
#include "simupop_cfg.h"
#include "utility.h"
#include "population.h"
#include "operator.h"

#include <string>
#include <algorithm>
using std::max_element;
using std::string;

namespace simuPOP
{
  /**
  The mating classes describe various mating scheme --- a required parameter
  of simulator.

  */

  class mating
  {
    public:

      /// numOffspring: constant, numOffspringFunc: call each time before mating
#define MATE_NumOffspring           1
      /// call numOffspringFunc each time during mating.
#define MATE_NumOffspringEachFamily 2
      /// numOffspring and numOffsrpingsFunc call each time before mating is
      /// the p for a geometric distribution
#define MATE_GeometricDistribution   3
#define MATE_PoissonDistribution     4
#define MATE_BinomialDistribution    5

    public:
      /// check if the mating type is compatible with population structure
      /**  possible things to check:
       -  need certain types of individual (age, sex etc)
       -  need resizeable population...
      */
      virtual bool isCompatible( population& pop)
      {
        return true;
      }

      /// constructor
      mating(double numOffspring=1.0,
        PyObject* numOffspringFunc=NULL,
        UINT maxNumOffspring = 0,
        UINT mode=MATE_NumOffspring,
        vectorlu newSubPopSize=vectorlu(),
        string newSubPopSizeExpr="",
        PyObject* newSubPopSizeFunc=NULL);

      mating(const mating& rhs)
        : m_numOffspring(rhs.m_numOffspring),
        m_numOffspringFunc(rhs.m_numOffspringFunc),
        m_maxNumOffspring(rhs.m_maxNumOffspring),
        m_mode(rhs.m_mode), m_firstOffspring(true),
        m_subPopSize(rhs.m_subPopSize),
        m_subPopSizeExpr(rhs.m_subPopSizeExpr),
        m_subPopSizeFunc(rhs.m_subPopSizeFunc)
      {
        if( m_subPopSizeFunc != NULL)
          Py_INCREF(m_subPopSizeFunc);
        if( m_numOffspringFunc != NULL)
          Py_INCREF(m_numOffspringFunc);
      }

      /// destructor
      virtual ~mating()
      {
        if( m_subPopSizeFunc != NULL)
          Py_DECREF(m_subPopSizeFunc);
        if( m_numOffspringFunc != NULL)
          Py_DECREF(m_numOffspringFunc);
      }

      /// clone() const. Generate a copy of itself and return a pointer
      /**
      This function is important because Python automatically
      release an object after it is used.

      For example:
      \code
        foo( mate = mating() )
        // in C++ implementation, foo keeps a pointer to mating()
        // mating * p = mate;
        foo.p->fun()
      \endcode
      will fail since mating() is released after the first line
      being executed.

      With the help of clone() const, the C++ implementation can avoid this problem by
      \code
      foo( mating* mate = &mating() )
      {
      mating * p = mate.clone() const;
      }
      \endcode
      */
      virtual mating* clone() const
      {
        return new mating(*this);
      }

      /// return name of the mating type
      /// used primarily in logging.
      virtual string __repr__()
      {
        return "<simuPOP::generic mating scheme>";
      }

      virtual void submitScratch(population& pop, population& scratch)
      {
      }

      /// mate: This is not supposed to be called for base mating class.
      /**
      \param pop population
      \param scratch scratch population
      \param ops during mating operators
      \return return false when mating fail.
      */
      virtual bool mate( population& pop, population& scratch, vector<Operator* >& ops, bool submit)
      {
        throw SystemError("You are not supposed to call base mating scheme.");
      }

      UINT numOffspring(int gen );

      void resetNumOffspring()
      {
        m_firstOffspring =  true;
      }

    public:
      /// whether or not to generate offspring genotype
      /// this is true when none of the during-mating operator can do this.
      bool formOffGenotype(const vector<Operator* >& ops);

      /// dealing with pop/subPop size change, copy of structure etc.
      void prepareScratchPop(population& pop, population& scratch);

    protected:

      /// number of offspring each mate
      double m_numOffspring;

      /// number of offspring func
      PyObject* m_numOffspringFunc;

      ///
      UINT m_maxNumOffspring;

      /// whether or not call m_numOffspringFunc each time
      UINT m_mode;

      ///
      bool m_firstOffspring;

      /// new subpopulation size. mostly used to 'keep' subPopsize
      /// after migration.
      vectorlu m_subPopSize;

      /// expression to evaluate subPopSize
      /// population size can change as a result of this
      /// e.g. "%popSize*1.3" whereas %popSize is predefined
      Expression m_subPopSizeExpr;

      /// the function version of the parameter.
      /// the python function should take one parameter (gen) and
      /// return a vector of subpop size.
      PyObject* m_subPopSizeFunc;
  };

  /**
     No mating. No subpopulation change.
     During mating operator will be applied, but
     the return values are not checked.
  */

  class noMating: public mating
  {
    public:

      /// constructor, no new subPopsize parameter
      noMating()
        :mating(1, NULL, 0, MATE_NumOffspring,
        vectorlu(),"",NULL)
      {

      }

      /// destructor
      ~noMating(){}

      /// clone() const. The same as mating::clone() const.
      /** \sa mating::clone() const
       */
      virtual mating* clone() const
      {
        return new noMating(*this);
      }

      /// return name of the mating type
      virtual string __repr__()
      {
        return "<simuPOP::no mating>";
      }

      virtual void submitScratch(population& pop, population& scratch)
      {
      }

      /**
       \brief do the mating. --- no mating :-)

       All individuals will be passed to during mating operators but
       no one will die (ignore during mating failing signal).
      */
      virtual bool mate( population& pop, population& scratch, vector<Operator *>& ops, bool submit);
  };

  /**
    binomial random selection

    No sex. Choose one individual from last generation.

    1. numOffspring protocol is honored
    2. population size changes are allowed
    3. selection is possible.

    So this works just like a sexless random mating.
    If ploidy is one, this is chromosomal mating.
  */

  class binomialSelection: public mating
  {
    public:

      /// constructor
      binomialSelection(double numOffspring=1.,
        PyObject* numOffspringFunc=NULL,
        UINT maxNumOffspring=0,
        UINT mode=MATE_NumOffspring,
        vectorlu newSubPopSize=vectorlu(),
        string newSubPopSizeExpr="",
        PyObject* newSubPopSizeFunc=NULL)
        :mating(numOffspring,
        numOffspringFunc, maxNumOffspring, mode,
        newSubPopSize, newSubPopSizeExpr,
        newSubPopSizeFunc),
        m_sampler(rng())
        {}

      /// destructor
      ~binomialSelection(){}

      /// clone() const. The same as mating::clone() const.
      /** \sa mating::clone() const
       */
      virtual mating* clone() const
      {
        return new binomialSelection(*this);
      }

      /// return name of the mating type
      virtual string __repr__()
      {
        return "<simuPOP::binomial random selection>";
      }

      virtual void submitScratch(population& pop, population& scratch)
      {
        pop.fitness().clear();
        // use scratch population,
        pop.pushAndDiscard(scratch);
        DBG_DO(DBG_MATING, pop.setIntVectorVar("famSizes", m_famSize));
      }

      /**
       \brief do the mating.
       \param pop population
       \param scratch scratch population, will be used in this mating scheme.
       \param ops during mating operators
       \return return false when mating fails.
      */
      virtual bool mate( population& pop, population& scratch, vector<Operator *>& ops, bool submit);

    protected:
      /// accumulative fitness
      Weightedsampler m_sampler;

#ifndef OPTIMIZED
      ///
      vectori m_famSize;
#endif
  };

  /**
  basic sexual random mating.

  Within each subpopulation, choose male and female randomly
  randmly get one copy of chromosome from father/mother.

  require: sexed individual; ploidy == 2

  apply during mating operators and put into the next generation.

  if ignoreParentsSex is set, parents will be chosen regardless of sex.

  Otherwise, male and female will be collected and be chosen randomly.

  If there is no male or female in a subpopulation,
  if m_UseSameSexIfUniSex is true, an warning will be generated and same
  sex mating (?) will be used
  otherwise, randomMating will return false.

  if there is no during mating operator to copy alleles, a direct copy
  will be used.
  */

  class randomMating : public mating
  {
    public:

      /// create a random mating scheme
      /**
      \param numOffspring, number of offspring or p in some modes
      \param numOffspringFunc, a python function that determine
       number of offspring or p depending on mode
      \param maxNumOffspring used when mode=MATE_BinomialDistribution
      \param mode one of  MATE_NumOffspring , MATE_NumOffspringEachFamily,
       MATE_GeometricDistribution, MATE_PoissonDistribution, MATE_BinomialDistribution
      \param newSubPopSize an array of subpopulation sizes, should have the same
      number of subpopulations as current population
      \param newSubPopSizeExpr an expression that will be evaluated as an array of subpop sizes
      \param newSubPopSizeFunc an function that have parameter gen and oldSize (current subpop
      size).
      \param contWhenUniSex continue when there is only one sex in the population, default to true

      */
      randomMating( double numOffspring=1.,
        PyObject* numOffspringFunc=NULL,
        UINT maxNumOffspring=0,
        UINT mode=MATE_NumOffspring,
        vectorlu newSubPopSize=vectorlu(),
        PyObject* newSubPopSizeFunc=NULL,
        string newSubPopSizeExpr="",
        bool contWhenUniSex=true)
        :mating(numOffspring,
        numOffspringFunc, maxNumOffspring, mode,
        newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc),
        m_contWhenUniSex(contWhenUniSex),
        m_maleIndex(0), m_femaleIndex(0),
        m_maleFitness(0), m_femaleFitness(0),
        m_malesampler(rng()), m_femalesampler(rng())
      {
      }

      /// destructor
      ~randomMating(){}

      /// clone() const. Generate a copy of itself and return pointer
      /// this is to make sure the object is persistent and
      /// will not be freed by python.
      virtual mating* clone() const
      {
        return new randomMating(*this);
      }

      virtual bool isCompatible( population& pop)
      {
        // test if individual has sex
        // if not, will yield compile time error.
        pop.indBegin()->sex();

#ifndef OPTIMIZED
        if( pop.ploidy() != 2 )
          cout << "Warning: This mating type only works with diploid population." << endl;
#endif

        return true;
      }

      /// return name of the mating type
      virtual string __repr__()
      {
        return "<simuPOP::sexual random mating>";
      }

      virtual void submitScratch(population& pop, population& scratch)
      {
        pop.fitness().clear();
        // use scratch population,
        pop.pushAndDiscard(scratch);
        DBG_DO(DBG_MATING, pop.setIntVectorVar("famSizes", m_famSize));
      }

      /// do the mating. parameters see mating::mate .
      /**
      Within each subpopulation, choose male and female randomly
      randmly get one copy of chromosome from father/mother.

      require: sexed individual; ploidy == 2

      apply during mating operators and put into the next generation.

      Otherwise, male and female will be collected and be chosen randomly.

      - If there is no male or female in a subpopulation,
      - if m_contWhenUniSex is true, an warning will be generated and same
      sex mating (?) will be used
      - otherwise, randomMating will return false.

      */
      virtual bool mate( population& pop, population& scratch, vector<Operator *>& ops, bool submit);

    protected:

      /// if no other sex exist in a subpopulation,
      /// same sex mating will occur if m_contWhenUniSex is set.
      /// otherwise, an exception will be thrown.
      bool m_contWhenUniSex;

      /// internal index to female/males.
      vectorlu m_maleIndex, m_femaleIndex;

      vectorf m_maleFitness, m_femaleFitness;

      // weighted sampler
      Weightedsampler m_malesampler, m_femalesampler;

#ifndef OPTIMIZED
      ///
      vectori m_famSize;
#endif
  };

  //
  // simulate trajectories of disease susceptibility loci using an extension of
  // the backward method described in Slatkin 2001.
  //
  // Parameters:
  //
  //    N      constant population size N
  //    NtFunc a python function that returns population size at each generation.
  //           gen is defined in reversed order. NtFunc(0) should be current
  //           generation number.
  //    freq   expected allele frequency of allele a.
  //    s      fitness for [AA, Aa, aa] or [AA, Aa, aa, BB, Bb, bb] for the multi-locus
  //           case. Assume constant selection pressure. s is default to [1,1,1].
  //           I.e, a neutral process.
  //    sFunc  a python function that returns selection pressure at each generation
  //           the function expects a single parameter gen which is defined
  //           in reversed order.
  //    T      maximum generation number. The process will terminate even if it
  //           can not reach allele zero after T generations. Default to 100,000,
  //           roughly 2,000,000 years which is longer than human history.
  //
  // Of course, you should specify only one of N/NtFunc and one of s/sFunc
  //
  // Tracking the allele frequency of allele a.
  //
  //
  vectorf FreqTrajectoryStoch( double freq, long N=0,
    PyObject* NtFunc=NULL, vectorf s=vectorf(), PyObject* sFunc=NULL,
    ULONG T=100000);

  class trajectory
  {
    public:
      trajectory(size_t nTraj=1):m_freqs(nTraj)
      {
      }

      int numTraj()
      {
        return m_freqs.size();
      }

      size_t maxLen();

      void setTraj(const vectorf& freq, size_t idx=0);

      vectorf traj(size_t idx=0)
      {
        DBG_FAILIF( idx >= m_freqs.size(), IndexError,
          "Index out of range");

        return m_freqs[idx];
      }

    private:
      // all frequencies
      vector<vector<double> > m_freqs;
  };

  //
  // simulate trajectories of disease susceptibility loci using an extension of
  // the backward method described in Slatkin 2001.
  //
  // Parameters:
  //
  //    N      constant population size N
  //    NtFunc a python function that returns population size at each generation.
  //           gen is defined in reversed order. NtFunc(0) should be current
  //           generation number.
  //    freq   expected allele frequencies of alleles of multiple unlinked loci
  //    s      constant fitness for [AA, Aa, aa, BB, Bb, bb ...]
  //    sFunc  a python function that returns selection pressure at each generation
  //           the function expects parameters gen and freq. gen is current generation
  //           number and freq is the allele frequency at all loci. This allows
  //           frequency dependent selection. gen is defined in reversed order.
  //    T      maximum generation number. The process will return if it can not
  //           reach allele zero after T generations, even if alleles are not fixed.
  //           Default to 100,000, roughly
  //           2,000,000 years, longer than human history.
  //
  // Of course, you should specify only one of N/NtFunc and one of s/sFunc
  //
  // Tracking the allele frequency of allele a.
  //
  //
  trajectory FreqTrajectoryMultiStoch( vectorf freq=vectorf(), long N=0,
    PyObject* NtFunc=NULL, vectorf s=vectorf(), PyObject* sFunc=NULL,
    ULONG T=100000);

  // simulate trajectory
  vectorf FreqTrajectorySelSim(
    double sel,                                   // strength of selection coef  ::8
    long Ne,                                      // effective population size ::9
    double freq,                                  // initial freq ::10
    double dom_h,                                 // strength of dominance ::27
    int selection                                 // selection ::5
    );

  /*simulate the sample path of the frequency of disease allele,
  conditional on non-extinction and non-fixation
  array of the mutant allele frequency is backword, start from present-day,
  end at the founding time of the disease
  But, simulation is forward, start from the past when only one copy
  of disease allele, until the
  present-day (if #disease allele <5, poisson dis. conditional
  on #disease >0, else, normal dis.)*/
  vectorf FreqTrajectoryForward(double lowbound, double highbound,
    int disAge, double grate, long N0, double seleCo);

  /**
    controlled mating
  */

  class controlledMating : public mating
  {
    public:

      /// Controlled mating,
      /// control allele frequency at a locus.
      /**
      \param mating a mating scheme.
      \param loci  loci at which allele frequency is controlled. Note that
        controlling the allele frequencies at several loci may take a long
        time.
      \param alleles alleles to control at each loci. Should have the
        same length as loci
      \param freqFunc frequency boundaries. If the length of the return
        value equals the size of loci, the range for loci will be
        [value0, value0+range], [value1, value1+range] etc.
        If the length of the return value is 2 times size of loci, it will
      be interpreted as [low1, high1, low2, high2 ...]
      */
      controlledMating(
        mating& matingScheme,
        vectori loci,
        vectori alleles,
        PyObject* freqFunc,
        double range = 0.01
        )
        : m_loci(loci), m_alleles(alleles),
        m_freqFunc(freqFunc),
        m_range(range)
      {
        DBG_FAILIF( m_loci.empty(), ValueError, "Have to specify a locus (or a loci) to control");

        DBG_FAILIF( m_alleles.empty(), ValueError, "Have to specify allele at each locus");

        DBG_FAILIF( m_loci.size() != m_alleles.size(), ValueError, "Should specify allele for each locus");

        if(m_freqFunc == NULL || !PyCallable_Check(m_freqFunc))
          throw ValueError("Please specify a valid frequency function");
        else
          Py_INCREF(m_freqFunc);

        m_matingScheme = matingScheme.clone();
      }

      /// CPPONLY
      controlledMating(const controlledMating& rhs)
        : m_loci(rhs.m_loci),
        m_alleles(rhs.m_alleles),
        m_freqFunc(rhs.m_freqFunc),
        m_range(rhs.m_range)
      {
        Py_INCREF(m_freqFunc);
        m_matingScheme = rhs.m_matingScheme->clone();
      }

      /// destructor
      ~controlledMating()
      {
        if( m_freqFunc != NULL)
          Py_DECREF(m_freqFunc);
        delete m_matingScheme;
      }

      /// clone() const. Generate a copy of itself and return pointer
      /// this is to make sure the object is persistent and
      /// will not be freed by python.
      virtual mating* clone() const
      {
        return new controlledMating(*this);
      }

      virtual bool isCompatible( population& pop)
      {
        return m_matingScheme->isCompatible(pop);
      }

      /// return name of the mating type
      virtual string __repr__()
      {
        return "<simuPOP::controlled mating>";
      }

      vectorlu countAlleles(population& pop, const vectori& loci, const vectori& alleles);

      virtual bool mate( population& pop, population& scratch, vector<Operator *>& ops, bool submit);

    private:

      /// mating scheme
      mating* m_matingScheme;

      /// loci at which mating is controlled.
      vectori m_loci;

      /// allele to be controlled at each locus
      vectori m_alleles;

      /// function that return an array of frquency range
      PyObject * m_freqFunc;

      /// range, used when m_freqFunc returns a vector of the same length as m_loci
      double m_range;
  };

  /**
    binomial random selection

    No sex. Choose one individual from last generation.

    1. numOffspring protocol is honored
    2. population size changes are allowed
    3. selection is possible.

    So this works just like a sexless random mating.
    If ploidy is one, this is chromosomal mating.
  */

  class controlledBinomialSelection: public binomialSelection
  {
    public:

      /// constructor
      controlledBinomialSelection(
        int locus,
        Allele allele,
        PyObject* freqFunc,
        double numOffspring=1.,
        PyObject* numOffspringFunc=NULL,
        UINT maxNumOffspring=0,
        UINT mode=MATE_NumOffspring,
        vectorlu newSubPopSize=vectorlu(),
        string newSubPopSizeExpr="",
        PyObject* newSubPopSizeFunc=NULL)
        :binomialSelection(numOffspring,
        numOffspringFunc, maxNumOffspring, mode,
        newSubPopSize, newSubPopSizeExpr,
        newSubPopSizeFunc),
        m_locus(locus),
        m_allele(allele),
        m_freqFunc(freqFunc)
      {
        if(m_freqFunc == NULL || !PyCallable_Check(m_freqFunc))
          throw ValueError("Please specify a valid frequency function");
        else
          Py_INCREF(m_freqFunc);
      }

      /// CPPONLY
      controlledBinomialSelection(const controlledBinomialSelection& rhs)
        : binomialSelection(rhs),
        m_locus(rhs.m_locus),
        m_allele(rhs.m_allele),
        m_freqFunc(rhs.m_freqFunc)
      {
        Py_INCREF(m_freqFunc);
      }

      /// destructor
      ~controlledBinomialSelection()
      {
        if( m_freqFunc != NULL)
          Py_DECREF(m_freqFunc);
      }

      /// clone() const. The same as mating::clone() const.
      /** \sa mating::clone() const
       */
      virtual mating* clone() const
      {
        return new controlledBinomialSelection(*this);
      }

      /// return name of the mating type
      virtual string __repr__()
      {
        return "<simuPOP::binomial random selection>";
      }

      virtual void submitScratch(population& pop, population& scratch)
      {
        pop.fitness().clear();
        // use scratch population,
        pop.pushAndDiscard(scratch);
        DBG_DO(DBG_MATING, pop.setIntVectorVar("famSizes", m_famSize));
      }

      /**
       \brief do the mating.
       \param pop population
       \param scratch scratch population, will be used in this mating scheme.
       \param ops during mating operators
       \return return false when mating fails.
      */
      virtual bool mate( population& pop, population& scratch, vector<Operator *>& ops, bool submit);

    private:
      /// locus at which mating is controlled.
      int m_locus;

      /// allele to be controlled at each locus
      Allele m_allele;

      /// function that return an array of frquency range
      PyObject * m_freqFunc;
  };

  /**
  basic sexual random mating.

  Within each subpopulation, choose male and female randomly
  randmly get one copy of chromosome from father/mother.

  require: sexed individual; ploidy == 2

  apply during mating operators and put into the next generation.

  if ignoreParentsSex is set, parents will be chosen regardless of sex.

  Otherwise, male and female will be collected and be chosen randomly.

  If there is no male or female in a subpopulation,
  if m_UseSameSexIfUniSex is true, an warning will be generated and same
  sex mating (?) will be used
  otherwise, randomMating will return false.

  if there is no during mating operator to copy alleles, a direct copy
  will be used.
  */

  class controlledRandomMating : public randomMating
  {
    public:

      /// create a random mating scheme
      /**
      \param numOffspring, number of offspring or p in some modes
      \param numOffspringFunc, a python function that determine
       number of offspring or p depending on mode
      \param maxNumOffspring used when mode=MATE_BinomialDistribution
      \param mode one of  MATE_NumOffspring , MATE_NumOffspringEachFamily,
       MATE_GeometricDistribution, MATE_PoissonDistribution, MATE_BinomialDistribution
      \param newSubPopSize an array of subpopulation sizes, should have the same
      number of subpopulations as current population
      \param newSubPopSizeExpr an expression that will be evaluated as an array of subpop sizes
      \param newSubPopSizeFunc an function that have parameter gen and oldSize (current subpop
      size).
      \param contWhenUniSex continue when there is only one sex in the population, default to true
      */
      controlledRandomMating(
        int locus,
        Allele allele,
        PyObject* freqFunc,
        double numOffspring=1.,
        PyObject* numOffspringFunc=NULL,
        UINT maxNumOffspring=0,
        UINT mode=MATE_NumOffspring,
        vectorlu newSubPopSize=vectorlu(),
        PyObject* newSubPopSizeFunc=NULL,
        string newSubPopSizeExpr="",
        bool contWhenUniSex=true)
        :randomMating(numOffspring,
        numOffspringFunc, maxNumOffspring, mode,
        newSubPopSize,
        newSubPopSizeFunc,
        newSubPopSizeExpr,
        contWhenUniSex),
        m_locus(locus),
        m_allele(allele),
        m_freqFunc(freqFunc)
      {
        if(m_freqFunc == NULL || !PyCallable_Check(m_freqFunc))
          throw ValueError("Please specify a valid frequency function");
        else
          Py_INCREF(m_freqFunc);
      }

      /// CPPONLY
      controlledRandomMating(const controlledRandomMating& rhs)
        : randomMating(rhs),
        m_locus(rhs.m_locus),
        m_allele(rhs.m_allele),
        m_freqFunc(rhs.m_freqFunc)
      {
        Py_INCREF(m_freqFunc);
      }

      /// destructor
      ~controlledRandomMating()
      {
        if( m_freqFunc != NULL)
          Py_DECREF(m_freqFunc);
      }

      /// clone() const. Generate a copy of itself and return pointer
      /// this is to make sure the object is persistent and
      /// will not be freed by python.
      virtual mating* clone() const
      {
        return new controlledRandomMating(*this);
      }

      virtual bool isCompatible( population& pop)
      {
        // test if individual has sex
        // if not, will yield compile time error.
        pop.indBegin()->sex();

#ifndef OPTIMIZED
        if( pop.ploidy() != 2 )
          cout << "Warning: This mating type only works with diploid population." << endl;
#endif

        return true;
      }

      /// return name of the mating type
      virtual string __repr__()
      {
        return "<simuPOP::sexual random mating>";
      }

      virtual void submitScratch(population& pop, population& scratch)
      {
        pop.fitness().clear();
        // use scratch population,
        pop.pushAndDiscard(scratch);
        DBG_DO(DBG_MATING, pop.setIntVectorVar("famSizes", m_famSize));
      }

      /// do the mating. parameters see mating::mate .
      /**
      Within each subpopulation, choose male and female randomly
      randmly get one copy of chromosome from father/mother.

      require: sexed individual; ploidy == 2

      apply during mating operators and put into the next generation.

      Otherwise, male and female will be collected and be chosen randomly.

      - If there is no male or female in a subpopulation,
      - if m_contWhenUniSex is true, an warning will be generated and same
      sex mating (?) will be used
      - otherwise, controlledRandomMating will return false.

      */
      virtual bool mate( population& pop, population& scratch, vector<Operator *>& ops, bool submit);

    private:

      /// locus at which mating is controlled.
      int m_locus;

      /// allele to be controlled at each locus
      Allele m_allele;

      /// function that return an array of frquency range
      PyObject * m_freqFunc;
  };
}
#endif
