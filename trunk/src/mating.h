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

      virtual void submitScratch(Pop& pop, Pop& scratch)
      {
      }

      /// mate: This is not supposed to be called for base Mating class.
      /**
      \param pop Population
      \param scratch scratch population
      \param ops during mating operators
      \return return false when mating fail.
      */
      virtual bool mate( Pop& pop, Pop& scratch, vector<Operator<Pop>* >& ops, bool submit)
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
          PyCallFunc(m_numOffspringFunc, "(i)", gen, numOS, PyObj_As_Double);
        }

        m_firstOffspring = false;

        if( m_mode == MATE_NumOffspring || m_mode == MATE_NumOffspringEachFamily)
        {
          DBG_FAILIF( numOS < 1, ValueError, "Need at least one offspring.");
          return static_cast<UINT>(numOS);
        }
        else if( m_mode == MATE_GeometricDistribution)
        {
          DBG_FAILIF( fcmp_lt(numOS, 0) || fcmp_gt(numOS, 1.), ValueError,
            "P for a geometric distribution should be within [0,1], given " + toStr(numOS));
          UINT nos = rng().randGeometric(numOS);
          return nos;
        }
        else if( m_mode == MATE_PoissonDistribution)
        {
          DBG_FAILIF( fcmp_lt(numOS, 0) || fcmp_gt(numOS, 1.), ValueError,
            "P for a Poisson distribution should be within [0,1], given " + toStr(numOS));
          UINT nos = rng().randPoisson(numOS)+1;
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

          vectorf res;
          PyCallFunc2(m_subPopSizeFunc, "(iO)", gen, curSize,
            res, PyObj_As_Array);
          Py_XDECREF(curSize);                    // just to be safe

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
     No mating. No subpopulation change.
     During mating operator will be applied, but
     the return values are not checked.
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

      virtual void submitScratch(Pop& pop, Pop& scratch)
      {
      }

      /**
       \brief do the mating. --- no mating :-)

       All individuals will be passed to during mating operators but
       no one will die (ignore during mating failing signal).
      */
      virtual bool mate( Pop& pop, Pop& scratch, vector<Operator<Pop> *>& ops, bool submit)
      {
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

    1. numOffspring protocol is honored
    2. population size changes are allowed
    3. selection is possible.

    So this works just like a sexless random mating.
    If ploidy is one, this is chromosomal mating.
  */
  template<class Pop>
    class BinomialSelection: public Mating<Pop>
  {
    public:

      /// constructor
      BinomialSelection(double numOffspring=1.,
        PyObject* numOffspringFunc=NULL,
        UINT maxNumOffspring=0,
        UINT mode=MATE_NumOffspring,
        vectorlu newSubPopSize=vectorlu(),
        string newSubPopSizeExpr="",
        PyObject* newSubPopSizeFunc=NULL)
        :Mating<Pop>(numOffspring,
        numOffspringFunc, maxNumOffspring, mode,
        newSubPopSize, newSubPopSizeExpr,
        newSubPopSizeFunc),
        m_sampler(rng())
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

      virtual void submitScratch(Pop& pop, Pop& scratch)
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
      virtual bool mate( Pop& pop, Pop& scratch, vector<Operator<Pop> *>& ops, bool submit)
      {
        this->resetNumOffspring();
        // scrtach will have the right structure.
        this->prepareScratchPop(pop, scratch);

        DBG_DO(DBG_MATING, m_famSize.clear());

        DBG_ASSERT( pop.numSubPop() == scratch.numSubPop(), SystemError,
          "Number of subpopulation can not be changed.");

        vectorf& fitness = pop.fitness();

        /// determine if mate() will generate offspring genotype
        bool formOffGeno = this->formOffGenotype(ops);

        // for each subPopulation
        for(UINT sp=0; sp < pop.numSubPop(); ++sp)
        {
          UINT spSize = pop.subPopSize(sp);
          if( spSize == 0 ) continue;

          // if selection is on
          if( ! fitness.empty() )
            m_sampler.set( vectorf(fitness.begin()+pop.subPopBegin(sp),
              fitness.begin()+pop.subPopEnd(sp) ) );

          // choose a parent and genreate m_numOffspring offspring
          ULONG spInd = 0;
          ULONG spIndEnd = scratch.subPopSize(sp);
          while( spInd < spIndEnd)
          {
            typename Pop::IndType * parent;
            // choose a parent
            if( !fitness.empty() )
              parent = &pop.individual( m_sampler.get(), sp);
            else
              parent = &pop.individual( rng().randInt(spSize), sp);

            // generate m_numOffspring offspring
            UINT numOS, numOSEnd;
            for(numOS=0, numOSEnd = this->numOffspring(pop.gen() ); numOS < numOSEnd;  numOS++)
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
                    numOS --;
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
              if( spInd == spIndEnd )
              {
                numOS++;
                break;
              }
            }                                     // offsrping for each parent
            // record family size
            DBG_DO(DBG_MATING, m_famSize.push_back( numOS ));
          }                                       // all offspring
        }                                         // all subPopulation.

        if(submit)
          submitScratch(pop, scratch);
        return true;
      }

    private:
      /// accumulative fitness
      WeightedSampler m_sampler;

#ifndef OPTIMIZED
      ///
      vectori m_famSize;
#endif
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
        string newSubPopSizeExpr="",
        bool contWhenUniSex=true)
        :Mating<Pop>(numOffspring,
        numOffspringFunc, maxNumOffspring, mode,
        newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc),
        m_contWhenUniSex(contWhenUniSex),
        m_maleIndex(0), m_femaleIndex(0),
        m_maleFitness(0), m_femaleFitness(0),
        m_maleSampler(rng()), m_femaleSampler(rng())
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

      virtual void submitScratch(Pop& pop, Pop& scratch)
      {
        pop.fitness().clear();
        // use scratch population,
        pop.pushAndDiscard(scratch);
        DBG_DO(DBG_MATING, pop.setIntVectorVar("famSizes", m_famSize));
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
      virtual bool mate( Pop& pop, Pop& scratch, vector<Operator<Pop> *>& ops, bool submit)
      {
        bool hasSexChrom = pop.sexChrom();

        this->resetNumOffspring();
        // scrtach will have the right structure.
        this->prepareScratchPop(pop, scratch);

        DBG_DO(DBG_MATING, m_famSize.clear());

        DBG_ASSERT( pop.numSubPop() == scratch.numSubPop(), SystemError,
          "Number of subpopulation can not be changed.");

        // empty fitness means no selection
        vectorf& fitness = pop.fitness();

        /// determine if any during-mating operator will generate offspring genotype
        bool formOffGeno = this->formOffGenotype(ops);

        UINT numMale, numFemale;

        /// random mating happens within each subPopulation
        for(UINT sp=0; sp < pop.numSubPop(); ++sp)
        {
          ULONG spSize = pop.subPopSize(sp);
          if( spSize == 0 ) continue;

          numMale = 0;
          for( typename Pop::IndIterator it=pop.indBegin(sp), itEnd = pop.indEnd(sp); it < itEnd;  ++it)
            if(it->sex() == Male)
              numMale ++;

          // to gain some performance, allocate memory at first.
          m_maleIndex.resize(numMale);
          m_femaleIndex.resize(spSize-numMale);

          numMale = 0;
          numFemale = 0;

          for( ULONG it=pop.subPopBegin(sp), itEnd =pop.subPopEnd(sp); it < itEnd;  it++)
          {
            if( pop.individual(it).sex() == Male)
              m_maleIndex[numMale++] = it;
            else
              m_femaleIndex[numFemale++] = it;
          }

          /// now, all individuals of needToFind sex is collected
          if( (numMale == 0 || numFemale ==0 ) && !m_contWhenUniSex )
          {
            throw ValueError("Subpopulation becomes uni-sex. Can not continue. \n"
              "You can use ignoreParentsSex (do not check parents' sex) or \ncontWhenUnixSex "
              "(same sex mating if have to) options to get around this problem.");
          }

          DBG_ASSERT( numFemale + numMale == spSize, SystemError,
            "Wrong number of male/female.");

          /// if selection is on
          if( ! fitness.empty() )
          {
            m_maleFitness.resize(numMale);
            m_femaleFitness.resize(numFemale);

            size_t ind;

            DBG_ASSERT( fitness.size() == pop.popSize(),
              ValueError, "Length of var fitness should equal to popsize");

            for( ind = 0; ind < numMale; ++ind)
              m_maleFitness[ind] = fitness[ m_maleIndex[ind] ];
            for( ind = 0; ind < numFemale; ++ind)
              m_femaleFitness[ind] = fitness[ m_femaleIndex[ind] ];

            m_maleSampler.set(m_maleFitness);
            m_femaleSampler.set(m_femaleFitness);
          }

          // generate scratch.subPopSize(sp) individuals.
          ULONG spInd = 0;
          ULONG spIndEnd = scratch.subPopSize(sp);
          while( spInd < spIndEnd)
          {
            // randomly choose parents
            typename Pop::IndType * dad, *mom;
            RNG& rnd = rng();

            if( !fitness.empty() )                // with selection
            {
              // using weidhted sampler.
              if( numMale != 0 )
                dad = &pop.individual( m_maleIndex[ m_maleSampler.get() ] );
              else
                dad = &pop.individual( m_femaleIndex[ m_femaleSampler.get() ] );

              if( numFemale != 0 )
                mom = &pop.individual( m_femaleIndex[ m_femaleSampler.get() ] );
              else
                mom = &pop.individual( m_maleIndex[ m_maleSampler.get() ]);
            }
            else
            {
              // using random sample.
              if( numMale != 0 )
                dad = &pop.individual( m_maleIndex[ rnd.randInt(numMale) ]);
              else
                dad = &pop.individual( m_femaleIndex[ rnd.randInt(numFemale) ]);

              if( numFemale != 0 )
                mom = &pop.individual( m_femaleIndex[ rnd.randInt(numFemale) ]);
              else
                mom = &pop.individual( m_maleIndex[ rnd.randInt(numMale) ]);
            }

            // generate m_numOffspring offspring per mating
            UINT numOS=0, numOSEnd;
            for(numOS=0, numOSEnd = this->numOffspring(pop.gen()); numOS < numOSEnd;  numOS++)
            {
              typename Pop::IndIterator it = scratch.indBegin(sp) + spInd++;

              //
              // assign sex randomly
              if( ! hasSexChrom)
              {
                int offSex = rnd.randInt(2);
                it->setSex( offSex==0?Male:Female);
              }

              if( formOffGeno )                   // use the default no recombination random mating.
              {
                int dadPloidy=0, momPloidy=1;     // initialize to avoid compiler complains

                for(UINT ch=0, chEnd = dad->numChrom(); ch < chEnd;  ++ch)
                {
                  dadPloidy = rnd.randInt(2);
                  momPloidy = rnd.randInt(2);

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

              /// apply all during mating operators
              for( typename vector<Operator<Pop> *>::iterator iop = ops.begin(), iopEnd = ops.end(); iop != iopEnd;  ++iop)
              {
                try
                {
                  // During Mating operator might reject this offspring.
                  if(!(*iop)->applyDuringMating(pop, it, dad, mom))
                  {
                    spInd --;
                    numOS --;
                    break;
                  }
                }
                catch(...)
                {
                  cout << "DuringMating operator " << (*iop)->__repr__() << " throws an exception." << endl << endl;
                  throw;
                }
              }
              // appy operators
              // success
              if( spInd == spIndEnd )
              {
                numOS++;
                break;
              }
            }                                     // each offspring
            // and then break here since spInd == spIndEnd.
            // record family size
            DBG_DO(DBG_MATING, m_famSize.push_back( numOS ));
          }
        }                                         // each subPop

        if(submit)
          submitScratch(pop, scratch);
        return true;
      }

    private:

      /// if no other sex exist in a subPopulation,
      /// same sex mating will occur if m_contWhenUniSex is set.
      /// otherwise, an exception will be thrown.
      bool m_contWhenUniSex;

      /// internal index to female/males.
      vectorlu m_maleIndex, m_femaleIndex;

      vectorf m_maleFitness, m_femaleFitness;

      // weighted sampler
      WeightedSampler m_maleSampler, m_femaleSampler;

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
  //    T      maximum generation number. The process will restart if it can not
  //           reach allele zero after T generations. Default to 1,000,000
  //
  // Of course, you should specify only one of N/NtFunc and one of s/sFunc
  //
  // Tracking the allele frequency of allele a.
  //
  //
  vectorf FreqTrajectoryStoch( double freq, long N=0,
    PyObject* NtFunc=NULL, vectorf s=vectorf(), PyObject* sFunc=NULL,
    ULONG T=1000000)
  {
    // is NtFunc callable?
    if( NtFunc != NULL )
    {
      if( ! PyCallable_Check(NtFunc) )
        throw ValueError("NtFunc is not a valid Python function.");
      else
        // increase the ref, just to be safe
        Py_INCREF(NtFunc);
    }

    // 1, 1+s1, 1+s2
    double s1, s2;
    if( sFunc!= NULL )
    {
      if( ! PyCallable_Check(sFunc) )
        throw ValueError("sFunc is not a valid Python function.");
      else
        // increase the ref, just to be safe
        Py_INCREF(sFunc);
    }
    else if( s.empty() )
      // default to [1,1,1]
    {
      s1 = 0.;
      s2 = 0.;
    }
    else if( s.size() != 3 || s[0] == 0)
    {
      throw ValueError("s should be a vector of length 3. (for AA, Aa and aa)");
    }
    else
    {
      // convert to the form 1, s1, s2
      s1 = s[1] / s[0] - 1.;
      s2 = s[2] / s[0] - 1.;
    }

    // get current population size
    int Ntmp;
    if( NtFunc == NULL)
      Ntmp = N;
    else
    {
      PyCallFunc(NtFunc, "(i)", 0, Ntmp, PyObj_As_Int);
    }

    // all calculated population size
    vectorlu Nt(1, Ntmp);
    // copies of allele a at each genertion.
    vectorlu it(1, static_cast<long>(Ntmp*freq));
    // allele frequency of allele a at each geneation
    vectorf xt(1, freq);
    // store calculated fitness s1, s2, if necessary
    // note that s_T will never be used.
    vectorf s1_cache(1,0), s2_cache(1,0);
    // will be used to store return value of sFunc
    vectorf s_vec;

    // t is current generation number.
    ULONG idx = 0;

    // a,b,c etc for solving the quadratic equation
    double a,b,c,b2_4ac,y1,y2,y;
    while( true )
    {
      // first get N(t-1), if it has not been calculated
      if( idx+1 >= Nt.size() )
      {
        if( NtFunc != NULL)
        {
          PyCallFunc(NtFunc, "(i)", idx+1, Ntmp, PyObj_As_Int);
        }
        Nt.push_back(Ntmp);
      }
      //
      // get fitness
      if( sFunc != NULL )
      {
        if( idx+1 >= s1_cache.size() )
        {
          PyCallFunc(sFunc, "(i)", idx+1, s_vec, PyObj_As_Array);

          DBG_ASSERT(s_vec.size()==3 || s_vec[0] != 0., ValueError,
            "Returned value from sFunc should be a vector of size 3");

          s1 = s_vec[1]/s_vec[0] - 1.;
          s2 = s_vec[2]/s_vec[0] - 1.;
          s1_cache.push_back(s1);
          s2_cache.push_back(s2);
        }
        else
        {
          s1 = s1_cache[idx];
          s2 = s2_cache[idx];
        }
      }
      // given x(t)
      // calculate y=x(t-1)' by solving an equation
      //
      //  x_t = y(1+s2 y+s1 (1-y))/(1+s2 y+2 s1 y(1-y))
      //
      if( s1 == 0. && s2 == 0. )
      {
        // special case when a = 0.
        y = xt[idx];
      }
      else
      {
        a = s2*xt[idx] - 2*s1*xt[idx] - s2 + s1;
        b = 2*s1*xt[idx] - 1 - s1;
        c = xt[idx];
        b2_4ac = b*b-4*a*c;

        if( b2_4ac < 0 )
          throw ValueError("Quadratic function does not yield a valid solution");

        y1 = (-b+sqrt(b2_4ac))/(2*a);
        y2 = (-b-sqrt(b2_4ac))/(2*a);

        // choose one of the solutions
        if( y1 < 0. || y1 > 1.0 )
        {
          if( y2 >= 0. && y2 <= 1.0)
            y = y2;
          else
            throw ValueError("None of the solutions is valid");
        }
        else
        {
          // if y1 is valid
          if( y2 >= 0. && y2 <= 1.0)
            throw ValueError("Both solutions are valid. Further decision is needed. y1: " +
              toStr(y1) + " y2: " + toStr(y2) + " xt: " + toStr(xt[idx]));
          else
            y = y1;
        }
      }

      if( it.size() < idx+2 )
      {
        it.push_back(0);
        xt.push_back(0);
      }
      // y is obtained, is the expected allele frequency for the next generation t-1
      it[idx+1] = rng().randBinomial( Nt[idx+1], y);
      xt[idx+1] = it[idx+1]/static_cast<double>(Nt[idx+1]);

      if( it[idx+1] == 0 )
      {
        // need 0, 1, ....,
        if( it[idx] == 1 )
          break;
        else
        {
          DBG_DO(DBG_MATING, cout << "Reaching 0, but next gen has more than 1 allele a" << endl);
          // restart
          idx = 0;
        }
      }
      else if( it[idx+1] == Nt[idx+1] )
        // when the allele get fixed, restart
      {
        idx = 0;
      }
      // if not done, but t already reaches T
      else if( idx == T )
      {
        cout << "Warning: reaching T gnerations. Try again." << endl;
        idx = 0;
      }
      else
        // go to next generation
        idx++;
    }
    // clean up
    if( NtFunc != NULL)
      Py_DECREF(NtFunc);
    if( sFunc != NULL)
      Py_DECREF(sFunc);

    // return a reversed version of xt, without trailing 0
    // also we may only need part of it since it has been
    // extended by previous attemps
    DBG_ASSERT( it[idx+1] == 0, SystemError,
      "The last generation still has allele a");

    // number of valid generation is idx+1
    vectorf traj(idx+1);
    for(ULONG i=0; i<=idx; ++i)
      traj[i] = xt[idx-i];
    return traj;
  }

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

      vectorf& setTraj(const vectorf freq, size_t idx=0)
      {
        DBG_FAILIF( idx >= m_freqs.size(), IndexError,
          "Index out of range");

        m_freqs[idx] = freq;
      }

      vectorf traj(size_t idx=0)
      {
        DBG_FAILIF( idx >= m_freqs.size(), IndexError,
          "Index out of range");

        return m_freqs[idx];
      }

      vectorf& trajRef(size_t idx=0)
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
  //    T      maximum generation number. The process will restart if it can not
  //           reach allele zero after T generations. Default to 1,000,000
  //
  // Of course, you should specify only one of N/NtFunc and one of s/sFunc
  //
  // Tracking the allele frequency of allele a.
  //
  //
  trajectory FreqTrajectoryMultiStoch( vectorf freq=vectorf(), long N=0,
    PyObject* NtFunc=NULL, vectorf s=vectorf(), PyObject* sFunc=NULL,
    ULONG T=1000000)
  {
    size_t nLoci = freq.size();
    size_t i, j;

    DBG_ASSERT( nLoci > 0, ValueError, "Number of loci should be at least one");

    trajectory result(nLoci);

    // in the cases of independent and constant selection pressure
    // easy case.
    if( sFunc == NULL)
    {
      for( i=0; i<nLoci; ++i)
      {
        result.setTraj(FreqTrajectoryStoch(freq[i], N, NtFunc,
          vectorf(s.begin()+3*i, s.begin()+3*(i+1)),
          NULL, T), i);
      }
      return result;
    }

    // other wise, use a vectorized version of ...
    // is NtFunc callable?
    if( NtFunc != NULL )
    {
      if( ! PyCallable_Check(NtFunc) )
        throw ValueError("NtFunc is not a valid Python function.");
      else
        // increase the ref, just to be safe
        Py_INCREF(NtFunc);
    }

    // sAll will store returned value of sFunc,
    // convert to 1, 1+s1, 1+s2 format ...
    vectorf sAll;
    if( sFunc != NULL )
    {
      if( ! PyCallable_Check(sFunc) )
        throw ValueError("sFunc is not a valid Python function.");
      else
        // increase the ref, just to be safe
        Py_INCREF(sFunc);
    }
    else if( s.empty() )
      // default to [1,1,1]
    {
      sAll.resize(nLoci*3, 0.);
    }
    else if( s.size() != 3*nLoci)
    {
      throw ValueError("s should be a vector of length 3 times number of loci. (for AA, Aa and aa etc)");
    }
    else
    {
      for(i=0; i<nLoci; ++i)
      {
        // convert to the form 1, s1, s2
        sAll[3*i+1] = s[3*i+1] / s[3*i] - 1.;
        sAll[3*i+2] = s[3*i+2] / s[3*i] - 1.;
        sAll[3*i] = 0.;
      }
    }

    // get current population size
    int Ntmp;
    if( NtFunc == NULL)
      Ntmp = N;
    else
    {
      PyCallFunc(NtFunc, "(i)", 0, Ntmp, PyObj_As_Int);
    }

    // all calculated population size
    vectorlu Nt(1, Ntmp);
    // copies of allele a at each genertion.
    vectorlu it(nLoci);
    for(i=0; i<nLoci; ++i)
      it[i] = static_cast<long>(Ntmp*freq[i]);
    // allele frequency of allele a at each geneation
    vectorf xt(nLoci);
    for(i=0; i<nLoci; ++i)
      xt[i] = freq[i];

    ULONG idx = 0;

    // a,b,c etc for solving the quadratic equation
    double s1,s2,x,a,b,c,b2_4ac,y1,y2,y;
    // whether or not each locus is done
    vector<bool> done(nLoci, false);
    //
    while( true )
    {
      // first get N(t-1), if it has not been calculated
      if( idx+1 >= Nt.size() )
      {
        if( NtFunc != NULL)
        {
          PyCallFunc(NtFunc, "(i)", idx+1, Ntmp, PyObj_As_Int);
        }
        Nt.push_back(Ntmp);
      }
      //
      // get fitness, since it will change according to
      // xt, I do not cache the result
      if( sFunc != NULL )
      {
        // compile allele frequency... and pass
        PyObject* freqObj = Double_Vec_As_NumArray( xt.begin()+nLoci*idx, xt.end()+nLoci*(idx+1) );
        PyCallFunc2(sFunc, "(iO)", idx+1, freqObj, sAll, PyObj_As_Array);

        DBG_ASSERT(sAll.size()==3*nLoci, ValueError,
          "Returned value from sFunc should be a vector of size 3");

        for(i=0; i<nLoci; ++i)
        {
          // convert to the form 1, s1, s2
          sAll[3*i+1] = sAll[3*i+1] / sAll[3*i] - 1.;
          sAll[3*i+2] = sAll[3*i+2] / sAll[3*i] - 1.;
          sAll[3*i] = 0.;
        }
      }
      //
      for(i = 0; i<nLoci; ++i)
      {
        // allocate space
        if( it.size() < (idx+2)*nLoci )
        {
          for(j=0;j<nLoci;++j)
          {
            it.push_back(0);
            xt.push_back(0);
          }
        }
        // if done
        if( done[i] )
        {
          it[(idx+1)*nLoci+i] = 0;
          xt[(idx+1)*nLoci+i] = 0.;
          continue;
        }
        // given x(t)
        // calculate y=x(t-1)' by solving an equation
        //
        //  x_t = y(1+s2 y+s1 (1-y))/(1+s2 y+2 s1 y(1-y))
        //
        s1 = sAll[i*3+1];
        s2 = sAll[i*3+2];
        x  = xt[idx*nLoci+i];
        if( s1 == 0. && s2 == 0. )
        {
          // special case when a = 0.
          y = x;
        }
        else
        {
          a = s2*x - 2*s1*x - s2 + s1;
          b = 2*s1*x - 1 - s1;
          c = x;
          b2_4ac = b*b-4*a*c;

          if( b2_4ac < 0 )
            throw ValueError("Quadratic function does not yield a valid solution");

          y1 = (-b+sqrt(b2_4ac))/(2*a);
          y2 = (-b-sqrt(b2_4ac))/(2*a);

          // choose one of the solutions
          if( y1 < 0. || y1 > 1.0 )
          {
            if( y2 >= 0. && y2 <= 1.0)
              y = y2;
            else
              throw ValueError("None of the solutions is valid");
          }
          else
          {
            // if y1 is valid
            if( y2 >= 0. && y2 <= 1.0)
              throw ValueError("Both solutions are valid. Further decision is needed. y1: " +
                toStr(y1) + " y2: " + toStr(y2) + " xt: " + toStr(xt[idx]));
            else
              y = y1;
          }
        }

        // y is obtained, is the expected allele frequency for the next generation t-1
        ULONG i_t = rng().randBinomial( Nt[idx+1], y);
        it[nLoci*(idx+1)+i] = i_t;
        xt[nLoci*(idx+1)+i] = i_t/static_cast<double>(Nt[idx+1]);

        if( i_t == 0 )
        {
          // need 0, 1, ...., good...
          if( it[idx] == 1 )
          {
            done[i] = true;
          }
          else
          {
            DBG_DO(DBG_MATING, cout << "Reaching 0, but next gen has more than 1 allele a" << endl);
            // restart
            idx = 0;
            break;
          }
        }
        else if( i_t == Nt[idx+1] )
          // when the allele get fixed, restart
        {
          idx = 0;
          break;
        }
      }                                           // end of for each locus
      // if all done?
      if( find(done.begin(), done.end(), false) == done.end() )
        break;
      //
      // if not done, but t already reaches T
      if( idx == T )
      {
        cout << "Warning: reaching T gnerations. Try again." << endl;
        idx = 0;
      }
      else
        // go to next generation
        idx ++;
    }
    // clean up
    if( NtFunc != NULL)
      Py_DECREF(NtFunc);
    if( sFunc != NULL)
      Py_DECREF(sFunc);

    // number of valid generation is idx+1
    vectorf traj(idx+1);
    for(i=0; i<nLoci; ++i)
    {
      for(j=0; j<=idx; ++j)
        traj[j] = xt[nLoci*(idx-j)+i];
      result.setTraj(traj, i);
    }
    return result;
  }

  // simulate trajectory
  vectorf FreqTrajectorySelSim(
    double sel,                                   // strength of selection coef  ::8
    long Ne,                                      // effective population size ::9
    double freq,                                  // initial freq ::10
    double dom_h,                                 // strength of dominance ::27
    int selection                                 // selection ::5
    )
  {
    //
    // first step, prepeare u:
    //
    // probability of zero for traj condition
    // originally prep_table
    long   N = Ne;
    double s = sel*4.0;
    double h = dom_h;

    vector<double> u(N+1);
    double ratio = ((1-s/4.0)/(1+s/4.0));
    vector<double> prod(N);
    prod[0] = 1;

    for(long i=1;i<N;i++)
    {
      if(selection==1)
        prod[i] = prod[i-1]*ratio;
      if(selection==2)
        prod[i] = prod[i-1]*((1-(s)*(1-2*((double)i/N)))/(1+(s)*(1-2*((double)i/N))));
      if(selection==3)
        prod[i] = prod[i-1]*( 1 - (s/2.0)*(((double)i/N)+h*(1-2*((double)i/N)) ) )/
          ( 1 + (s/2.0)*(((double)i/N)+h*(1-2*((double)i/N)) ) );
    }
    u[N-1] = prod[N-1];
    for(long i = N-2; i>0; i--)
      u[i] = u[i+1] + prod[i];
    double first= u[1];
    for(long i=1; i<N; i++)
      u[i] /= (1.0 + first);

    //
    // second step
    //
    // prepare log look up table, necessary?
    vector<double> log_lookup(N+1);
    for(int i=1;i<=N;i++)
      log_lookup[i] = -(double)N*log((double)i/(i+1));

    //
    // third step
    //
    // start
    int jump=-1;
    double weight_up;

    int j = int(double(N)*freq);
    if((j == 0)||(j>=N))
      throw " Starting number selected types = 0 or >=N ";

    double * pos_log_ptr = &log_lookup[j];
    double * neg_log_ptr = &log_lookup[N-j-1];
    double * u0_pnt_j   = &u[j];
    double * u0_pnt_jp1 = &u[j+1];

    double lambda_driving = 1+s/4.0;
    double tot_birth = 2.0;
    double prob_up = lambda_driving/tot_birth;

    // accumulative t_time
    double past_traj_time = 0;

    vector<double> traj_time(1, 0);
    vector<double> traj_freq(1, double(j)/N);
    vector<double> traj_integral_pos(1, 0);
    vector<double> traj_integral_neg(1, 0);

    double prev_pos  = 0;
    double prev_neg  = 0;
    double t_time = 0;

    traj_time.push_back(0);
    traj_freq.push_back(0);

    double position_coeff;

    while(true)
    {
      traj_integral_pos.back() = prev_pos + t_time*(*pos_log_ptr);
      traj_integral_neg.back() = prev_neg + t_time*(*neg_log_ptr);

      position_coeff = (double)j * (N-j)/N;

      //Removed Factor of 2
      t_time = -(1/position_coeff)*log(rng().randUniform01() ) /double(N);
      past_traj_time += t_time;
      weight_up = (*u0_pnt_jp1)/(*u0_pnt_j);

      if(selection==2)
        prob_up = (1+(s)*(1-2*((double)j/N))) / 2.0;
      if(selection==3)
        prob_up = ( 1 + (s/2.0)*(((double)j/N) +
          h*(1-2*(double(j)/N)) ) ) / 2.0 ;

      //Jump
      if( rng().randUniform01() < prob_up*weight_up)
      {
        if(jump==1)
        {
          pos_log_ptr++;
          neg_log_ptr--;
        }
        j++;
        u0_pnt_jp1++;
        u0_pnt_j++;
        jump = 1;
      }
      else
      {
        if(jump==-1)
        {
          pos_log_ptr--;
          neg_log_ptr++;
        }
        j--;
        u0_pnt_jp1--;
        u0_pnt_j--;
        jump=-1;
      }

      prev_pos = traj_integral_pos.back();
      prev_neg = traj_integral_neg.back();
      traj_time.back() = past_traj_time;
      traj_freq.back() = (double)j/N;

      if ( j == 0)
        break;
      if ( j == N )
        cout << "Reach allele N, this should be rare" << endl;

      if ( traj_freq.size() % 10000 == 1)
        cout << "Size " << traj_freq.size() << " at freq " << traj_freq.back() << endl;

      traj_freq.push_back(0);
      traj_time.push_back(0);
      traj_integral_pos.push_back(0);
      traj_integral_neg.push_back(0);
    }

    // now we need to translate 4Nt to generations
    // g = 4Nt, g is generation
    double maxTime = traj_time.back();
    vectorf gen_freq( static_cast<size_t>(maxTime * 4 * Ne)+1);
    size_t ngen = gen_freq.size();
    int curTime = traj_time.size()-1;

    for(size_t g=0; g<ngen; ++g)
    {
      double t = (static_cast<size_t>(maxTime*4*Ne) - g)/(4.*Ne);
      // find, from the right, most close to t
      // idea case: x< t <cur
      while( curTime >= 1 && traj_time[curTime-1] > t)
        --curTime;
      gen_freq[g] = traj_freq[curTime-1] +
        (traj_freq[curTime]-traj_freq[curTime-1])*(t-traj_time[curTime-1]);
    }

    return gen_freq;
  }

  /*simulate the sample path of the frequency of disease allele,
  conditional on non-extinction and non-fixation
  array of the mutant allele frequency is backword, start from present-day,
  end at the founding time of the disease
  But, simulation is forward, start from the past when only one copy
  of disease allele, until the
  present-day (if #disease allele <5, poisson dis. conditional
  on #disease >0, else, normal dis.)*/
  vectorf FreqTrajectoryForward(double lowbound, double highbound,
    int disAge, double grate, long N0, double seleCo)
  {
    vector<double> DisSamplePath(disAge+1);

    int ftime = disAge;
    // DisSamplePath is the current allele frequency
    DisSamplePath[0] = 0;

    int trying=0;
    while(DisSamplePath[0] <= lowbound || DisSamplePath[0] >= highbound)
    {
      if((++trying)%1000==0)
        cout<<"Trying path "<<trying<<" times\n";
      double num = 1;
      //initially 1 copy of mutant allele
      double Nt = N0*exp(-ftime*grate);
      double fre = 1/(2*Nt);

      // initial allele frequency
      DisSamplePath[ftime] = fre;

      // backward in array, but forward in time.
      for(int gth = ftime-1; gth>= 0; gth--)
      {
        // population size at generation gth.
        Nt = N0*exp(-grate*gth);
        // number of alleles is less than four, use
        // poisson approximation.
        if(num<4)
        {
          double lamda;
          // expected allele frequency: sp(1-p)
          lamda = (fre + seleCo*fre*(1-fre))*2*Nt;
          num = rng().randPoisson(lamda);
          // restart if num<=0
          if(num<= 0)
          {
            DisSamplePath[0] = 0;
            break;
          }
          fre = num/(2*Nt);
          // restart?
          if(fre>= 1)
          {
            DisSamplePath[0] = 1;
            break;
          }
          DisSamplePath[gth] = fre;
        }
        // diffusion process approximation
        else
        {
          // use uniform approximation???
          double min = seleCo*fre*(1-fre)-sqrt(3*fre*(1-fre)/(2*Nt));
          double max = seleCo*fre*(1-fre)+sqrt(3*fre*(1-fre)/(2*Nt));
          double psv = min + (max-min)*rng().randUniform01();
          fre = fre+psv;
          if(fre<= 0)
          {
            DisSamplePath[0] = 0;
            break;
          }
          if(fre>= 1)
          {
            DisSamplePath[0] = 1;
            break;
          }
          num = fre* 2* Nt;
          DisSamplePath[gth] = fre;
        }
      }
    }                                             // while
    // reverse the result and return
    vectorf gen_freq(disAge+1);
    for(int i=0; i<disAge+1; ++i)
      gen_freq[i] = DisSamplePath[disAge-i];

    return gen_freq;
  }

  /**
    controlled mating
  */
  template<class Pop>
    class ControlledMating : public Mating<Pop>
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
      ControlledMating(
        Mating<Pop>& matingScheme,
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
      ControlledMating(const ControlledMating& rhs)
        : m_loci(rhs.m_loci),
        m_alleles(rhs.m_alleles),
        m_freqFunc(rhs.m_freqFunc),
        m_range(rhs.m_range)
      {
        Py_INCREF(m_freqFunc);
        m_matingScheme = rhs.m_matingScheme->clone();
      }

      /// destructor
      ~ControlledMating()
      {
        if( m_freqFunc != NULL)
          Py_DECREF(m_freqFunc);
        delete m_matingScheme;
      }

      /// clone() const. Generate a copy of itself and return pointer
      /// this is to make sure the object is persistent and
      /// will not be freed by python.
      virtual Mating<Pop>* clone() const
      {
        return new ControlledMating(*this);
      }

      virtual bool isCompatible( Pop& pop)
      {
        return m_matingScheme->isCompatible(pop);
      }

      /// return name of the mating type
      virtual string __repr__()
      {
        return "<simuPOP::controlled mating>";
      }

      vectorlu countAlleles(Pop& pop, const vectori& loci, const vectori& alleles)
      {
        vectorlu alleleNum(loci.size(), 0L);
        for(size_t l=0; l < loci.size(); ++l)
        {
          int loc = loci[l];
          Allele ale = alleles[l];
          // go through all alleles
          for(GappedAlleleIterator a=pop.alleleBegin(loc),
            aEnd=pop.alleleEnd(loc); a != aEnd; ++a)
          {
            if( AlleleUnsigned(*a) == ale )
              alleleNum[l]++;
          }
        }
        return alleleNum;
      }

      virtual bool mate( Pop& pop, Pop& scratch, vector<Operator<Pop> *>& ops, bool submit)
      {
        // first call the function and get the range
        vectorf freqRange;

        PyCallFunc( m_freqFunc, "(i)", pop.gen(),
          freqRange, PyObj_As_Array);

        DBG_ASSERT( freqRange.size() == m_loci.size() || freqRange.size() == 2 * m_loci.size(),
          ValueError, "Length of returned range should have the same or double the number of loci.");

        vectorlu alleleNum;

#ifndef OPTIMIZED
        // calculate allele frequen at these loci
        alleleNum = countAlleles(pop, m_loci, m_alleles);

        DBG_DO(DBG_MATING, cout << "Range of allele frequencies at generation "
          << pop.gen() << " is " << freqRange << endl);
#endif

        // compare integer is easier and more acurate, and we need to make sure
        // there is one allele when the frequency is greater than 0.
        vectorlu alleleRange(2*m_loci.size(), 0);
        if( freqRange.size() == m_loci.size() )
        {
          for(size_t i=0; i<m_loci.size(); ++i)
          {
            if (freqRange[i] < 0. )
              freqRange[i] = 0.;
            if (freqRange[i] > 1. )
              freqRange[i] = 1.;

            alleleRange[2*i] = static_cast<ULONG>(freqRange[i]*pop.popSize()*pop.ploidy());
            if( freqRange[i]>0 && alleleRange[2*i] == 0)
              alleleRange[2*i] = 1;
            // upper bound using range.
            alleleRange[2*i+1] = static_cast<ULONG>((freqRange[i]+m_range)*pop.popSize()*pop.ploidy())+1;
#ifndef OPTIMIZED
            if( alleleNum[i] == 0 && alleleRange[2*i] > 0 )
              throw ValueError("No allele exists so there is no way to reach specified allele frequency.");
#endif
          }
        }
        else
          // returned values are [low1, high1, low2, high2...]
        {
          for(size_t i=0; i<m_loci.size(); ++i)
          {
            if (freqRange[2*i] < 0. )
              freqRange[2*i] = 0.;
            if (freqRange[2*i+1] > 1. )
              freqRange[2*i+1] = 1.;

            alleleRange[2*i] = static_cast<ULONG>(freqRange[2*i]*pop.popSize()*pop.ploidy());

            DBG_FAILIF( freqRange[2*i] > freqRange[2*i+1], ValueError,
              "Incorrect frequency range: " + toStr(freqRange[2*i]) +
              " - " + toStr(freqRange[2*i+1]));

            if( freqRange[2*i]>0 && alleleRange[2*i] == 0)
              alleleRange[2*i] = 1;
            alleleRange[2*i+1] = static_cast<ULONG>(freqRange[2*i+1]*pop.popSize()*pop.ploidy())+1;
#ifndef OPTIMIZED
            if( alleleNum[i] == 0 && alleleRange[2*i] > 0 )
              throw ValueError("No allele exists so there is no way to reach specified allele frequency.");
#endif
          }
        }

        while(true)
        {
          // do the mating
          m_matingScheme->mate(pop, scratch, ops, false);

          // check allele frequency
          alleleNum = countAlleles(scratch, m_loci, m_alleles);

          DBG_DO(DBG_MATING, cout << "Mating finished, new count "
            << alleleNum << " range " << alleleRange << endl);
          //
          bool succ = true;
          for(size_t i=0; i<m_loci.size(); ++i)
          {
            if( alleleNum[i] < alleleRange[2*i] || alleleNum[i] > alleleRange[2*i+1] )
            {
              succ = false;
              break;
            }
          }
          // cout << "succ " << succ << "submit " << submit << endl;
          if( succ && submit)
          {
            // cout << "success" << endl;
            m_matingScheme->submitScratch(pop, scratch);
            break;
          }
        }
        return true;
      }

    private:

      /// mating scheme
      Mating<Pop>* m_matingScheme;

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
                PyObject* popObj=pyPopObj((void*)(&pop));
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
