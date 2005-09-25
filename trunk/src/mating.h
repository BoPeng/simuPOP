/***************************************************************************
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

#ifndef _MATING_H
#define _MATING_H
/**
\file
\brief head file of class Mating and its subclasses
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
  The Mating classes describe various mating scheme --- a required parameter
  of Simulator.

  */
  template<class Pop>
    class Mating
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
      virtual bool isCompatible( Pop& pop)
      {
        return true;
      }

      /// constructor
      Mating(double numOffspring=1.0,
        PyObject* numOffspringFunc=NULL,
        UINT maxNumOffspring = 0,
        UINT mode=MATE_NumOffspring,
        vectorlu newSubPopSize=vectorlu(),
        string newSubPopSizeExpr="",
        PyObject* newSubPopSizeFunc=NULL)
        :m_numOffspring(numOffspring),
        m_numOffspringFunc(NULL),
        m_maxNumOffspring(maxNumOffspring),
        m_mode(mode), m_firstOffspring(true),
        m_subPopSize(newSubPopSize),
        m_subPopSizeExpr(newSubPopSizeExpr,""),
        m_subPopSizeFunc(NULL)
      {
        DBG_FAILIF( !m_subPopSizeExpr.empty() && newSubPopSizeFunc != NULL,
          ValueError, "Please only specify one of newSubPopSizeExpr and newSubPopSizeFunc.");

        if ( numOffspringFunc != NULL)
        {
          if(!PyCallable_Check(numOffspringFunc))
            throw ValueError("Passed variable is not a callable python function.");

          Py_INCREF(numOffspringFunc);
          m_numOffspringFunc = numOffspringFunc;
        }

        if( newSubPopSizeFunc != NULL)
        {
          if(!PyCallable_Check(newSubPopSizeFunc))
            throw ValueError("Passed variable is not a callable python function.");

          Py_INCREF(newSubPopSizeFunc);
          m_subPopSizeFunc = newSubPopSizeFunc;
        }

        DBG_FAILIF(mode == MATE_BinomialDistribution && maxNumOffspring < 2,
          ValueError, "If mode is MATE_BinomialDistribution, maxNumOffspring should be > 1");
      }

      Mating(const Mating& rhs)
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
      virtual ~Mating()
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
        foo( mate = Mating() )
        // in C++ implementation, foo keeps a pointer to Mating()
        // Mating * p = mate;
        foo.p->fun()
      \endcode
      will fail since Mating() is released after the first line
      being executed.

      With the help of clone() const, the C++ implementation can avoid this problem by
      \code
      foo( Mating* mate = &Mating() )
      {
      Mating * p = mate.clone() const;
      }
      \endcode
      */
      virtual Mating* clone() const
      {
        return new Mating(*this);
      }

      /// return name of the mating type
      /// used primarily in logging.
      virtual string __repr__()
      {
        return "<simuPOP::generic mating scheme>";
      }

      /// mate: This is not supposed to be called for base Mating class.
      /**
      \param pop Population
      \param scratch scratch population
      \param ops during mating operators
      \return return false when mating fail.
      */
      virtual bool mate( Pop& pop, Pop& scratch, vector<Operator<Pop>* >& ops)
      {
        throw SystemError("You are not supposed to call base mating scheme.");
      }

      UINT numOffspring(int gen )
      {
        static double numOS = 0.;

        // use the same numOffspings each generation
        if( ! m_firstOffspring && m_mode == MATE_NumOffspring )
          return static_cast<UINT>(numOS);

        if(m_numOffspringFunc == NULL)
          numOS = m_numOffspring;
        else if( m_mode == MATE_NumOffspring ||   // case 1, first time
          // case 2, all the time; case 3,4,5, first time
          m_mode == MATE_NumOffspringEachFamily ||
          ( m_firstOffspring &&
          ( m_mode == MATE_GeometricDistribution  ||
          m_mode == MATE_PoissonDistribution   ||
          m_mode == MATE_BinomialDistribution )))
        {
          PyObject* arglist = Py_BuildValue("(i)", gen);
          PyObject* result = PyEval_CallObject(m_numOffspringFunc, arglist);
          Py_XDECREF(arglist);

          if( result == NULL)
          {
            PyErr_Print();
            throw ValueError("NumOfOffspring function call failed.");
          }

          PyObj_As_Double(result, numOS);
          Py_DECREF(result);
        }

        m_firstOffspring = false;

        if( m_mode == MATE_NumOffspring || m_mode == MATE_NumOffspringEachFamily)
        {
          DBG_FAILIF( numOS < 1, ValueError, "Need at least one offspring.");
          DBG_DO(DBG_MATING, cout << "Number of offspring " <<  static_cast<UINT>(numOS) << endl);
          return static_cast<UINT>(numOS);
        }
        else if( m_mode == MATE_GeometricDistribution)
        {
          DBG_FAILIF( fcmp_lt(numOS, 0) || fcmp_gt(numOS, 1.), ValueError,
            "P for a geometric distribution should be within [0,1], given " + toStr(numOS));
          UINT nos = rng().randGeometric(numOS);
          DBG_DO(DBG_MATING, cout << "Number of offspring " << nos << endl);
          return nos;
        }
        else if( m_mode == MATE_PoissonDistribution)
        {
          DBG_FAILIF( fcmp_lt(numOS, 0) || fcmp_gt(numOS, 1.), ValueError,
            "P for a Poisson distribution should be within [0,1], given " + toStr(numOS));
          UINT nos = rng().randPoisson(numOS)+1;
          DBG_DO(DBG_MATING, cout << "Number of offspring " << nos << endl);
          return nos;
        }
        else if( m_mode == MATE_BinomialDistribution)
        {
          DBG_FAILIF( fcmp_lt(numOS, 0) || fcmp_gt(numOS, 1.), ValueError,
            "P for a Bionomial distribution should be within [0,1], given " + toStr(numOS));
          DBG_FAILIF( m_maxNumOffspring < 1, ValueError,
            "Max number of offspring should be greater than 1. Given "
            + toStr(m_maxNumOffspring));
          UINT nos = rng().randBinomial(m_maxNumOffspring-1, numOS)+1;
          DBG_DO(DBG_MATING, cout << "Number of offspring " << nos << endl);
          return nos;
        }
        else
          throw ValueError("Wrong mating numoffspring mode. Should be one of \n"
            "MATE_NumOffspring, MATE_NumOffspringEachFamily and MATE_GEometricDistribution");
      }

      void resetNumOffspring()
      {
        m_firstOffspring =  true;
      }

    public:
      /// whether or not to generate offspring genotype
      /// this is true when none of the during-mating operator can do this.
      bool formOffGenotype(const vector<Operator<Pop>* >& ops)
      {
        bool res = true;
        for( typename vector<Operator<Pop> *>::const_iterator iop = ops.begin(),
          iopEnd = ops.end(); iop != iopEnd;  ++iop)
        {
          if( (*iop)->formOffGenotype() )
          {
            DBG_DO(DBG_SIMULATOR, cout << "Operator " << (*iop)->__repr__() <<
              " will form offspring" << endl);

            res= false;
            break;
          }
        }
        return res;

      }

      /// dealing with pop/subPop size change, copy of structure etc.
      void prepareScratchPop(Pop& pop, Pop& scratch)
      {
        /// force new subPop size?
        if( m_subPopSize.empty() && m_subPopSizeExpr.empty()
          && m_subPopSizeFunc == NULL )
        {
          // use population structure of pop
          // genotype of scratch pop will not be kept.
          scratch.setSubPopStru(pop.subPopSizes(), true);
        }
        else if(! m_subPopSize.empty())           // set subPoplation size
        {
          // allow  change of pop size of scratch
          scratch.setSubPopStru(this->m_subPopSize, true);

          DBG_FAILIF( scratch.numSubPop() != pop.numSubPop() ,
            ValueError, "number of subPopulaitons must agree. \n Pre: "
            + toStr( pop.numSubPop()) + " now: " + toStr(scratch.numSubPop() ));
        }
        else if(! m_subPopSizeExpr.empty())       // evaluate subPopulation size
        {
          m_subPopSizeExpr.setLocalDict( pop.dict());
          vectorf sizef = m_subPopSizeExpr.valueAsArray();
          vectorlu sz(sizef.size());

          for(size_t i=0, iEnd = sizef.size(); i<iEnd; i++)
            sz[i] = static_cast<ULONG>(sizef[i]);

          DBG_DO(DBG_SIMULATOR, cout << "New subpop size " << sz << endl);

          // allow population size change
          scratch.setSubPopStru(sz, true);
          DBG_FAILIF( scratch.numSubPop() != pop.numSubPop() ,
            ValueError, "number of subPopulaitons must agree.\n Pre: "
            + toStr( pop.numSubPop()) + " now: " + toStr(scratch.numSubPop() ));
        }
        else                                      // use m_subPopSizeFunc
        {
          // get generation number
          int gen = pop.gen();
          //          int gen = mainVars().getVarAsInt("gen");
          // convert current pop size to a tuple
          PyObject* curSize = PyTuple_New( pop.numSubPop());

          DBG_ASSERT(curSize!=NULL, SystemError, "Can not convert current pop size to a list");

          for(size_t i=0; i<pop.numSubPop(); ++i)
            PyTuple_SetItem(curSize, i, PyInt_FromLong( pop.subPopSize(i)));

          PyObject* arglist = Py_BuildValue("(iO)", gen, curSize);
          PyObject* result = PyEval_CallObject(m_subPopSizeFunc, arglist);
          Py_XDECREF(curSize);                    // just to be safe
          Py_XDECREF(arglist);

          if( result == NULL)
          {
            PyErr_Print();
            throw ValueError("Function call failed.");
          }

          vectorf res;
          PyObj_As_Array(result, res);
          Py_DECREF(result);

          vectorlu sz(res.size());

          for(size_t i=0; i< res.size(); i++)
            sz[i] = static_cast<ULONG>(res[i]);

          DBG_DO(DBG_SIMULATOR, cout << "New subpop size " << sz << endl);

          // allow  change of pop size of scratch
          scratch.setSubPopStru(sz, true);
          DBG_FAILIF( scratch.numSubPop() != pop.numSubPop() ,
            ValueError, "number of subPopulaitons must agree.\n Pre: "
            + toStr( pop.numSubPop()) + " now: " + toStr(scratch.numSubPop() ));
        }
      }

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

      /// new subPopulation size. mostly used to 'keep' subPopsize
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
     No mating. parent generation will be the offspring generation
     during mating operators will be applied though.

     */
  template<class Pop>
    class NoMating: public Mating<Pop>
  {
    public:

      /// constructor, no new subPopsize parameter
      NoMating()
        :Mating<Pop>(1, NULL, 0, MATE_NumOffspring,
        vectorlu(),"",NULL)
      {

      }

      /// destructor
      ~NoMating(){}

      /// clone() const. The same as Mating::clone() const.
      /** \sa Mating::clone() const
       */
      virtual Mating<Pop>* clone() const
      {
        return new NoMating(*this);
      }

      /// return name of the mating type
      virtual string __repr__()
      {
        return "<simuPOP::no mating>";
      }

      /**
       \brief do the mating. --- no mating :-)

       All individuals will be passed to during mating operators but
       no one will die (ignore during mating failing signal).
      */
      virtual bool mate( Pop& pop, Pop& scratch, vector<Operator<Pop> *>& ops)
      {
        // first: copy structure. Make sure scratch and pop are the same except
        // for genotype
        // apply during mating operators
        if( ! ops.empty() )
        {
          for(typename Pop::IndIterator it = pop.indBegin(), itEnd = pop.indEnd(); it != itEnd;  ++it)
          {
            for( typename vector<Operator<Pop> *>::iterator iop = ops.begin(), iopEnd = ops.end(); iop != iopEnd;  ++iop)
            {
              (*iop)->applyDuringMating(pop, it, NULL, NULL);
            }                                     // all during-mating operators
          }
        }
        return true;
      }

  };

  /**
    binomial random selection

    No sex. Choose one individual from last generation.
    Subpopulations are dealt separately.

    */
  template<class Pop>
    class BinomialSelection: public Mating<Pop>
  {
    public:

      /// constructor
      BinomialSelection(int numOffspring=1,
        vectorlu newSubPopSize=vectorlu(),
        string newSubPopSizeExpr="",
        PyObject* newSubPopSizeFunc=NULL)
        :Mating<Pop>(numOffspring, NULL, 0,
        MATE_NumOffspring, newSubPopSize, newSubPopSizeExpr,
        newSubPopSizeFunc)
        {}

      /// destructor
      ~BinomialSelection(){}

      /// clone() const. The same as Mating::clone() const.
      /** \sa Mating::clone() const
       */
      virtual Mating<Pop>* clone() const
      {
        return new BinomialSelection(*this);
      }

      /// return name of the mating type
      virtual string __repr__()
      {
        return "<simuPOP::binomial random selection>";
      }

      /**
       \brief do the mating.
       \param pop population
       \param scratch scratch population, will be used in this mating scheme.
       \param ops during mating operators
       \return return false when mating fails.
      */
      virtual bool mate( Pop& pop, Pop& scratch, vector<Operator<Pop> *>& ops)
      {

        this->resetNumOffspring();

        // get generation number
        //int gen = mainVars().getVarAsInt("gen");
        int gen = pop.gen();

        this->prepareScratchPop(pop, scratch);

        /// determine if mate() will generate offspring genotype
        bool formOffGeno = this->formOffGenotype(ops);

        bool select = false;

        // if selection is on
        // false: return 1 if selector does not exist.
        if( pop.getVarAsInt("selector", false) == 1 )
        {
          select = true;

          pop.setIntVar("selector", 0);
#ifndef OPTIMIZED
          // allow migration etc later.
          pop.setIntVar("fixIndOrder", 0);
#endif
        }

        // for each subPopulation
        for(UINT sp=0, spEnd = pop.numSubPop(); sp < spEnd;  ++sp)
        {
          UINT spSize = pop.subPopSize(sp);

          // if selection is on
          // generate accumulative fitness values
          WeightedSampler ws(rng());
          if( select )
          {
            m_fitness.resize(spSize);

            // get fitness values
            double* fitness;
#ifndef OPTIMIZED
            int len = pop.getVarAsDoubleNumArray("fitness", fitness);

            DBG_ASSERT( static_cast<size_t>(len) == pop.popSize(), ValueError,
              "Length of var fitness should equal to popsize");
#else
            pop.getVarAsDoubleNumArray("fitness", fitness);
#endif

            size_t ind;
            m_fitness[0] = fitness[ pop.subPopBegin(sp) ];
            if( spSize > 1 )
              for( ind = 1; ind < spSize; ++ind)
                m_fitness[ind] = m_fitness[ind-1] + fitness[ pop.subPopBegin(sp) + ind ];

            for(ind = 0; ind < spSize; ++ind)
              m_fitness[ind] /= m_fitness[spSize-1];
            ws.set(m_fitness);
          }

          // choose a parent and genreate m_numOffspring offspring
          ULONG spInd = 0;
          while( spInd < scratch.subPopSize(sp))
          {
            typename Pop::IndType * parent;
            // choose a parent
            if( select )
              parent = (&*pop.indBegin(sp) + ws.get());
            else
              parent = (&*pop.indBegin(sp) + rng().randInt(spSize));

            // generate m_numOffspring offspring
            for(UINT numOS=0, numOSEnd = this->numOffspring(gen); numOS < numOSEnd;  numOS++)
            {
              typename Pop::IndIterator it = scratch.indBegin(sp) + spInd++;

              if(  formOffGeno )
                /// use deep copy!!!!!!!
                it->copyFrom(*parent);

              // apply during mating operators
              for( typename vector<Operator<Pop> *>::iterator iop = ops.begin(),
                iopEnd = ops.end(); iop != iopEnd;  ++iop)
              {

                try
                {
                  // During Mating operator might reject this offspring.
                  if(!(*iop)->applyDuringMating(pop, it, parent, NULL))
                  {
                    spInd --;
                    break;
                  }
                }
                catch(...)
                {
                  cout << "Duringmating operator " << (*iop)->__repr__() << " throws an exception." << endl << endl;
                  throw;
                }
              }                                   // all during-mating operators

              // success
              if( spInd == scratch.subPopSize(sp) )
                break;
            }                                     // offsrping for each parent
          }                                       // all offspring
        }                                         // all subPopulation.
        // use scratch population,
        pop.pushAndDiscard(scratch);
        return true;
      }

    private:
      /// accumulative fitness
      vectorf m_fitness;
  };

  /**
  basic sexual random mating.

  Within each subPopulation, choose male and female randomly
  randmly get one copy of chromosome from father/mother.

  require: sexed individual; ploidy == 2

  apply during mating operators and put into the next generation.

  if ignoreParentsSex is set, parents will be chosen regardless of sex.

  Otherwise, male and female will be collected and be chosen randomly.

  If there is no male or female in a subPopulation,
  if m_UseSameSexIfUniSex is true, an warning will be generated and same
  sex mating (?) will be used
  otherwise, RandomMating will return false.

  if there is no during mating operator to copy alleles, a direct copy
  will be used.
  */

  template<class Pop>
    class RandomMating : public Mating<Pop>
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
      RandomMating( double numOffspring=1.,
        PyObject* numOffspringFunc=NULL,
        UINT maxNumOffspring=0,
        UINT mode=MATE_NumOffspring,
        vectorlu newSubPopSize=vectorlu(),
        PyObject* newSubPopSizeFunc=NULL,
        string newSubPopSizeExpr="", bool contWhenUniSex=true)
        :Mating<Pop>(numOffspring, numOffspringFunc, maxNumOffspring, mode,
        newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc),
        m_contWhenUniSex(contWhenUniSex),
        m_sexIndex(0),
        m_maleFitness(0), m_femaleFitness(0),
        m_maleSampler(rng()), m_femaleSampler(rng()),
        m_bt(rng())
      {
      }

      /// destructor
      ~RandomMating(){}

      /// clone() const. Generate a copy of itself and return pointer
      /// this is to make sure the object is persistent and
      /// will not be freed by python.
      virtual Mating<Pop>* clone() const
      {
        return new RandomMating(*this);
      }

      virtual bool isCompatible( Pop& pop)
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

      /// do the mating. parameters see Mating::mate .
      /**
      Within each subPopulation, choose male and female randomly
      randmly get one copy of chromosome from father/mother.

      require: sexed individual; ploidy == 2

      apply during mating operators and put into the next generation.

      Otherwise, male and female will be collected and be chosen randomly.

      - If there is no male or female in a subPopulation,
      - if m_contWhenUniSex is true, an warning will be generated and same
      sex mating (?) will be used
      - otherwise, RandomMating will return false.

      */
      virtual bool mate( Pop& pop, Pop& scratch, vector<Operator<Pop> *>& ops)
      {

        bool hasSexChrom = pop.sexChrom();

        this->resetNumOffspring();

        // get generation number
        int gen = pop.gen();

        this->prepareScratchPop(pop, scratch);

        /// determine if any during-mating operator will generate offspring genotype
        bool formOffGeno = this->formOffGenotype(ops);

        // if need to do recombination here
        size_t btIdx=0;
        if( formOffGeno && m_bt.size() != pop.numChrom()*2*scratch.popSize() )
        {
          vectorf rates(1);
          // for r = 0.5, 'fast bernulli' algorithm does not do
          // any good at all.... :-)
          // may be we should use the old way?
          rates[0] = 0.5;
          // 2 is for mom/dad
          m_bt.setParameter(rates, pop.numChrom()*2*scratch.popSize());
        }
        DBG_DO(DBG_MATING, cout << "doTrial " << endl);
        if( formOffGeno )
          m_bt.doTrial();

        bool select = false;

        /// if selection is on
        if( pop.getVarAsInt("selector", false) == 1 )
        {
          select = true;
          pop.setIntVar("selector", 0);
#ifndef OPTIMIZED
          // allow migration etc later.
          pop.setIntVar("fixIndOrder", 0);
#endif
        }

        // to save some resizing time, use constant size male and female index
        // there might be a big waste of RAM!!!
        ULONG maxSpSize = 0;
        for(size_t sp=0; sp < pop.numSubPop(); ++sp)
          if( maxSpSize < pop.subPopSize(sp) )
            maxSpSize = pop.subPopSize(sp);
        m_sexIndex.resize( maxSpSize );
        UINT numMale, numFemale;
        ULONG* maleIndex=&m_sexIndex[0], *femaleIndex;

        /// random mating happens within each subPopulation
        for(UINT sp=0, spEnd = pop.numSubPop(); sp < spEnd;  ++sp)
        {

          DBG_DO(DBG_MATING, cout << "SP " << sp << endl);
          ULONG spSize = pop.subPopSize(sp);

          if( spSize == 0 ) continue;

          numMale = 0;
          numFemale = maxSpSize;

          for( ULONG it=pop.subPopBegin(sp), itEnd =pop.subPopEnd(sp); it < itEnd;  it++)
          {
            if( pop.individual(it).sex() == Male)
              m_sexIndex[numMale++] = it;
            else
              m_sexIndex[--numFemale] = it;
          }

          // point to the fist female index
          femaleIndex = &m_sexIndex[numFemale];
          numFemale = maxSpSize - numFemale;

          /// now, all individuals of needToFind sex is collected
          if( (numMale == 0 || numFemale ==0 ) && !m_contWhenUniSex )
          {
            throw ValueError("Subpopulation becomes uni-sex. Can not continue. \n"
              "You can use ignoreParentsSex (do not check parents' sex) or \ncontWhenUnixSex "
              "(same sex mating if have to) options to get around this problem.");
          }

          /// if selection is on
          if( select )
          {
            m_maleFitness.resize(numMale);
            m_femaleFitness.resize(numFemale);

            size_t ind;

            // get fitness values
            double* fitness;
#ifndef OPTIMIZED
            int len = pop.getVarAsDoubleNumArray("fitness", fitness);

            DBG_ASSERT( static_cast<size_t>(len) == pop.popSize(),
              ValueError, "Length of var fitness should equal to popsize");
#else
            // does not check fitness length
            pop.getVarAsDoubleNumArray("fitness", fitness);
#endif

            for( ind = 0; ind < numMale; ++ind)
              m_maleFitness[ind] = fitness[ maleIndex[ind] ];
            for( ind = 0; ind < numFemale; ++ind)
              m_femaleFitness[ind] = fitness[ femaleIndex[ind] ];

            /*
            /// normalize and cusum fitness
            if( numMale > 1 )
              for( ind = 1; ind < numMale; ++ind)
                m_maleFitness[ind] += m_maleFitness[ind-1];
            for(ind = 0; ind < numMale; ++ind)
              m_maleFitness[ind] /= m_maleFitness[numMale-1];

            if( numFemale > 1 )
              for( ind = 1; ind < numFemale; ++ind)
                m_femaleFitness[ind] += m_femaleFitness[ind-1];
            for( ind = 0; ind < numFemale; ++ind)
            m_femaleFitness[ind] /= m_femaleFitness[numFemale-1];
            */
            m_maleSampler.set(m_maleFitness);
            m_femaleSampler.set(m_femaleFitness);
          }
          ULONG spInd = 0;
          ULONG spIndEnd = scratch.subPopSize(sp);
          while( spInd < spIndEnd)
          {
            // randomly choose parents
            typename Pop::IndType * dad, *mom;

            if( select)
            {
              if( numMale != 0 )
                dad = &pop.individual( maleIndex[ m_maleSampler.get() ] );
              else
                dad = &pop.individual( femaleIndex[ m_femaleSampler.get() ] );

              if( numFemale != 0 )
                mom = &pop.individual( femaleIndex[ m_femaleSampler.get() ] );
              else
                mom = &pop.individual( maleIndex[ m_maleSampler.get() ]);
            }
            else
            {
              if( numMale != 0 )
                dad = &pop.individual( maleIndex[ rng().randInt(numMale) ]);
              else
                dad = &pop.individual( femaleIndex[ rng().randInt(numFemale) ]);

              if( numFemale != 0 )
                mom = &pop.individual( femaleIndex[ rng().randInt(numFemale) ]);
              else
                mom = &pop.individual( maleIndex[ rng().randInt(numMale) ]);
            }

            // generate m_numOffspring offspring
            for(UINT numOS=0, numOSEnd = this->numOffspring(gen); numOS < numOSEnd;  numOS++)
            {
              typename Pop::IndIterator it = scratch.indBegin(sp) + spInd++;

              //
              // assign sex randomly
              if( ! hasSexChrom)
              {
                int offSex = rng().randInt(2);
                it->setSex( offSex==0?Male:Female);
              }

              if( formOffGeno )                   // use the default no recombination random mating.
              {
                int dadPloidy=0, momPloidy=1;     // initialize to avoid compiler complains
                const BitSet& succ = m_bt.succ(0);

                for(UINT ch=0, chEnd = dad->numChrom(); ch < chEnd;  ++ch)
                {
                  dadPloidy = succ[btIdx++];
                  momPloidy = succ[btIdx++];

                  DBG_ASSERT((dadPloidy==0 || dadPloidy==1) &&
                    ( momPloidy==0 || momPloidy==1), ValueError,
                    "Ploidy must be 0 or 1");

                  copy(dad->genoBegin(dadPloidy, ch),
                    dad->genoEnd(dadPloidy,ch) ,
                    it->genoBegin(0,ch));
                  copy(mom->genoBegin(momPloidy,ch),
                    mom->genoEnd(momPloidy,ch) ,
                    it->genoBegin(1,ch));
                }

                // last chromosome (sex chromosomes)
                if( hasSexChrom )
                {
                  if( dadPloidy == 1)             // Y of XY
                    it->setSex(Male);
                  else
                    it->setSex(Female);
                }
              }

              for( typename vector<Operator<Pop> *>::iterator iop = ops.begin(), iopEnd = ops.end(); iop != iopEnd;  ++iop)
              {
                try
                {
                  // During Mating operator might reject this offspring.
                  if(!(*iop)->applyDuringMating(pop, it, dad, mom))
                  {
                    spInd --;
                    break;
                  }
                }
                catch(...)
                {
                  cout << "Duringmating operator " << (*iop)->__repr__() << " throws an exception." << endl << endl;
                  throw;
                }
              }

              // appy operators
              // success
              if( spInd == spIndEnd )
                break;

            }                                     // each offspring
          }
        }                                         // each subPop

        // use scratch populaiton
        pop.pushAndDiscard(scratch);
        return true;
      }

    private:

      /// if no other sex exist in a subPopulation,
      /// same sex mating will occur if m_contWhenUniSex is set.
      /// otherwise, an exception will be thrown.
      bool m_contWhenUniSex;

      /// internal index to female/males.
      vectorlu m_sexIndex;

      vectorf m_maleFitness, m_femaleFitness;

      // weighted sampler
      WeightedSampler m_maleSampler, m_femaleSampler;

      // for use when there is RandomMating need to do recombination
      BernulliTrials m_bt;
  };

  /**
  random mating
  */
  /*
    class RandomMating: public Mating
    {
      public:
        /// whether or not the mating scheme will change population
        /// size.
        virtual bool willChangePopSize(){ return false; }

        /// constructor
        RandomMating(){};

  /// destructor
  ~RandomMating(){};

  /// return name of the mating type
  string __repr__(){ return "<simuPOP::random mating>";
  }

  bool operator()( Pop& pop, vector<Operator<Pop> *>& ops) {
  // get current population

  // get temporary population

  return true;
  }
  };
  */
  /**
  Polygamy
  */
  /*
  class Polygamy: public Mating
  {
    public:
      /// whether or not the mating scheme will change population
      /// size.
      virtual bool willChangePopSize(){ return false; }

      /// constructor
      Polygamy(){};

  /// destructor
  ~Polygamy(){};

  /// return name of the mating type
  string __repr__(){ return "<simuPOP::polygamy>";
  }

  bool operator()( Pop& pop, vector<Operator<Pop> *>& ops) {
  // get current population

  // get temporary population

  return true;
  }
  };
  */
  /*
  class Monogamy:public Mating
  {
  public:
      /// whether or not the mating scheme will change population
      /// size.
      virtual bool willChangePopSize(){ return false; }

      /// constructor
      Monogamy(){};

  /// destructor
  ~Monogamy(){};

  /// return name of the mating type
  string __repr__(){ return "<simuPOP::monogamy>";
  }

  bool operator()( Pop& pop, vector<Operator<Pop> *>& ops){
  // get current population

  // get temporary population

  return true;
  }
  };

  class Polygyny:public Mating
  {
  public:
  /// whether or not the mating scheme will change population
  /// size.
  virtual bool willChangePopSize(){ return false; }

  /// constructor
  Polygyny(){};

  /// destructor
  ~Polygyny(){};

  /// return name of the mating type
  string __repr__(){ return "<simuPOP::polygyny>";
  }

  bool operator()( Pop& pop, vector<Operator<Pop> *>& ops){
  // get current population

  // get temporary population

  return true;
  }
  };

  class Polyandry: public Mating
  {
  public:
  /// whether or not the mating scheme will change population
  /// size.
  virtual bool willChangePopSize(){ return false; }

  /// constructor
  Polyandry(){};

  /// destructor
  ~Polyandry(){};

  /// return name of the mating type
  string __repr__(){ return "<simuPOP::polyandry>";
  }

  bool operator()( Pop& pop, vector<Operator<Pop> *>& ops){
  // get current population

  // get temporary population

  return true;
  }
  };
  */

  /**
     Hybrid mating scheme.

     */
  template<class Pop>
    class PyMating: public Mating<Pop>
  {
    public:

      /// constructor
      /** This operator takes only one parameter: a mate function.
      During mating, this function will be called with pop as parameter.
      The expected return value should be a list (or carray) of indices
      to parents in the order of dad1,mom1,dad2,mom2,...

      \param mateFunc a python function that takes a population as parameter and
      return a list of parental indices. Currently, only diploid population is
      supported.

      \param keepSubPopStru if true, mating is strictly between subpop and
      subpopulation structure will be kept. Otherwise, mating is considered
      as in the big population regardless of subpopulation strcture. The resulting
      population does not have subpopulation strcture.

      Note that these indices should be absolte indices and mating across
      subpops is not allowed.

      In this way, you can organize arbitrary complex mating
      scheme (but also with considerable work.)
      */
      PyMating(PyObject* mateFunc, bool keepSubPopStru=true)
        :Mating<Pop>(1, NULL, 0, MATE_NumOffspring, vectorlu(),"",NULL),
        m_mateFunc(NULL), m_keepSubPopStru(keepSubPopStru)
      {
        if( ! PyCallable_Check(mateFunc) )
          throw ValueError("mateFunc is not a valid Python function.");

        Py_XINCREF(mateFunc);
        m_mateFunc = mateFunc;
      }

      /// destructor
      ~PyMating()
      {
        if( m_mateFunc != NULL)
          Py_DECREF(m_mateFunc);
      }

      PyMating(const PyMating& rhs):
      Mating<Pop>(rhs),
        m_mateFunc(rhs.m_mateFunc),
        m_keepSubPopStru(rhs.m_keepSubPopStru)
      {
        if( m_mateFunc != NULL)
          Py_INCREF(m_mateFunc);
      }

      /// clone() const. The same as Mating::clone() const.
      /** \sa Mating::clone() const
       */
      virtual Mating<Pop>* clone() const
      {
        return new PyMating(*this);
      }

      /// return name of the mating type
      virtual string __repr__()
      {
        return "<simuPOP::hybrid mating scheme>";
      }

      /**
       \brief do the mating with specified mating function.

       All individuals will be passed to during mating operators but
       no one will die (ignore during mating failing signal).
      */
      virtual bool mate( Pop& pop, Pop& scratch, vector<Operator<Pop> *>& ops)
      {
        throw SystemError("PyMating is not implemented yet.");
        /*
                /// determine if mate() will generate offspring genotype
                bool formOffGeno = this->formOffGenotype(ops);
                // if need to do recombination here
                size_t btIdx=0;
                // for use when there is RandomMating need to do recombination
                BernulliTrials bt(rng());

                // first: call the function, determine the
                // size of next generation.
                PyObject* popObj=pointer2pyObj((void*)(&pop), PopSWIGType);
        Py_INCREF(popObj);
        PyObject* arglist = Py_BuildValue("(O)", popObj );
        PyObject* result = PyEval_CallObject(m_mateFunc, arglist);

        // get result (a list)
        if( ! PySequence_Check(result))
        throw ValueError("Return value of mateFunc should be a sequence.");

        // get population size
        ULONG size = PySequence_Size(result);
        vectorlu spSize = vectorlu(1);
        spSize[0] = size/2;

        // prepare scratch pop
        scratch.setSubPopStru(spSize, true);
        // keep track of subpopulation stru
        spSize[0] = 0;

        if( formOffGeno && bt.size() != pop.numChrom()*2*scratch.popSize() )
        {
        vectorf rates(1);
        // for r = 0.5, 'fast bernulli' algorithm does not do
        rates[0] = 0.5;
        // 2 is for mom/dad
        bt.setParameter(rates, pop.numChrom()*2*scratch.popSize());
        }
        DBG_DO(DBG_MATING, cout << "doTrial " << endl);
        if( formOffGeno )
        bt.doTrial();

        //
        ULONG indIdx = 0;
        while( indIdx < size/2)
        {
        ULONG dadIdx = PyInt_AsLong( PySequence_GetItem(result, indIdx*2) );
        ULONG momIdx = PyInt_AsLong( PySequence_GetItem(result, indIdx*2+1) );

        if( m_keepSubPopStru)
        {
        ULONG spIdx = pop.subPopIndPair(dadIdx).first;
        ULONG spIdx1 = pop.subPopIndPair(momIdx).first;

        DBG_ASSERT( spIdx == spIdx1, ValueError,
        "Two parents should come from the same subpopulation.");

        if( spIdx >= spSize.size())
        spSize.resize(spIdx+1, 0);

        // the offspring is in this subpop.
        spSize[spIdx]++;
        }

        typename Pop::IndType * dad = &*pop.indBegin() + dadIdx ;
        typename Pop::IndType * mom = &*pop.indBegin() + momIdx ;
        typename Pop::IndIterator it = scratch.indBegin() + indIdx;
        //
        // assign sex randomly
        int offSex = rng().randInt(2);
        it->setSex( offSex==0?Male:Female);

        if( formOffGeno )
        {
        int dadPloidy;
        int momPloidy;
        const BitSet& succ = bt.succ(0);

        for(UINT ch=0, chEnd = dad->numChrom(); ch < chEnd;  ++ch)
        {

        //dadPloidy = rng().randInt(2);
        //momPloidy = rng().randInt(2);
        dadPloidy = succ[btIdx++];
        momPloidy = succ[btIdx++];

        DBG_ASSERT((dadPloidy==0 || dadPloidy==1) &&
        ( momPloidy==0 || momPloidy==1), ValueError,
        "Ploidy must be 0 or 1");

        copy(dad->genoBegin(dadPloidy, ch),
        dad->genoEnd(dadPloidy,ch) ,
        it->genoBegin(0,ch));
        copy(mom->genoBegin(momPloidy,ch),
        mom->genoEnd(momPloidy,ch) ,
        it->genoBegin(1,ch));
        }
        }

        // apply during mating operators
        for( typename vector<Operator<Pop> *>::iterator iop = ops.begin(),
        iopEnd = ops.end(); iop != iopEnd;  ++iop)
        {
        try
        {
        (*iop)->applyDuringMating(pop, it, dad, mom);
        }
        catch(...)
        {
        cout << "Duringmating operator " << (*iop)->__repr__() << " throws an exception." << endl << endl;
        throw;
        }
        }                                       // all during-mating operators
        // success
        indIdx++;

        }                                         // offsrping for each parent

        if(m_keepSubPopStru)
        scratch.setSubPopStru( spSize, false);

        Py_DECREF(popObj);
        Py_DECREF(arglist);

        // use scratch populaiton
        pop.pushAndDiscard(scratch);
        */
        return true;
      }

    private:
      PyObject* m_mateFunc;

      bool m_keepSubPopStru;
  };

}
#endif
