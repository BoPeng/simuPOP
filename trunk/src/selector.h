/***************************************************************************
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

#ifndef _SELECTOR_H
#define _SELECTOR_H
/**
\file
\brief head file of class selector:public Operator
*/
#include "utility.h"
#include "operator.h"

#include <numeric>
using std::min;

namespace simuPOP
{
  /** \brief selection

  Genetic selection is tricky to simulate. In simuPOP, I employee
  an ability (fitness) to mate approach. Namely, the probability
  that an individual will be chosen for mating is proportional
  to its fitness value. More specifically,

  - PreMating selectors assign fitness values to each individual.

  - Sexless mating (e.g. binomialSelection) : individuals are chosen
  at probabilities that are proportional to their fitness values. More
  specifically, if there are N individuals with fitness values
  \f$f_i, i=1,...,N \f$, individual \f$i\f$ will have probability
  \f$ \frac{f_i}{\sum_{j=1}^N f_j} \f$ to be chosen to be passed
  to the next generation.

  - Random mating with sex (e.g. randommating): males and females are
  separated and each are chosen as described above.

  Please refer to the user's guide for details.
  */

  class selector: public Operator
  {
    public:
      /// constructor. default to be always active.
      selector( int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL)
        :Operator("","",stage, begin, end, step, at, rep, grp)
      {
      }

      /// destructor
      virtual ~selector()
      {
      }

      virtual Operator* clone() const
      {
        return new selector(*this);
      }

      /// calculate/return w11 etc
      virtual double indFitness(individual * ind)
      {
        ///
        throw ValueError("This selector is not supposed to be called directly");
        return 1.;
      }

      /// set fitness to all individual
      bool apply(population& pop)
      {
        vectorf& fitness = pop.fitness();

        fitness.resize(pop.popSize());

        size_t i=0;
        for (population::IndIterator it = pop.indBegin(); it != pop.indEnd(); ++it)
          fitness[i++] = indFitness(&*it) ;

        return true;
      }

      virtual string __repr__()
      {
        return "<simuPOP::selector>" ;
      }
  };

  /** \brief selection according to genotype at one locus

  map selector. Assign fitness value according to a
  given dictionary.
  */

  class mapSelector: public selector
  {
    public:
      /** \brief create a map selector (selection according to genotype at one locus

      \param locus the locus index. The genotype of this locus will be axamed.
      \param loci the loci index. The genotype of this locus will be axamed.
      \param fitness a dictionary of fitness. The genotype must be in the form of 'a-b'
         for single locus, and 'a-b|c-d|e-f' for multi-locus..
      \param phase if true, a/b and b/a will have different fitness value. Default to false.
      \param output and other parameters please refer to help(baseOperator.__init__)
      */
      mapSelector( vectoru loci, const strDict& fitness, bool phase=false,
        int stage=PreMating, int begin=0, int end=-1, int step=1,
        vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL):
      selector(stage, begin, end, step, at, rep, grp),
        m_loci(loci), m_dict(fitness), m_phase(phase)
      {
      };

      virtual ~mapSelector()
      {
      }

      virtual Operator* clone() const
      {
        return new mapSelector(*this);
      }

      /// currently assuming diploid
      virtual double indFitness(individual * ind);

      virtual string __repr__()
      {
        return "<simuPOP::selector::map selector>" ;
      }

    private:
      /// one locus
      vectoru m_loci;

      /// fitness for each genotype
      strDict m_dict;

      ///
      bool m_phase;
  };

  /** \brief selection according to genotype at one locus

  multiple allele selector. This selector group alleles to disease
  and wild type and return fitness to AA,Aa,aa. (A is wildtype).
  */

  class maSelector: public selector
  {
    public:
      /** \brief create a multiple allele selector (selection according to diseased or wildtype
      alleles)
      Note that maSelector only work for diploid population now.

      \param locus the locus index. The genotype of this locus will be axamed.
      \param loci the loci index.
      \param fitness For the single locus case, fitness is an array of fitness of AA,Aa,aa. A is the
         wild type group. In the case of multiple loci, fitness should be in the order of
              BB Bb bb
           AA 1  2  3
           Aa 4  5  6
      aa 7  8  9
      The length for such table is 3^(#loci).
      \param wildtype an array of alleles in the wildtype group. Anything else is disease allele.
      default = [1]
      NOTE that wildtype at all loci are the same.
      \param output and other parameters please refer to help(baseOperator.__init__)
      */
      maSelector( vectoru loci, const vectorf& fitness, const vectora& wildtype,
        int stage=PreMating, int begin=0, int end=-1, int step=1,
        vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL):
      selector(stage, begin, end, step, at, rep, grp),
        m_loci(loci), m_fitness(fitness), m_wildtype(wildtype)
      {
        DBG_ASSERT( m_fitness.size() == static_cast<UINT>(pow(3, loci.size())),
          ValueError, "Please specify fitness for each combination of genotype.");
      };

      virtual ~maSelector()
      {
      }

      virtual Operator* clone() const
      {
        return new maSelector(*this);
      }

      /// currently assuming diploid
      virtual double indFitness(individual * ind);

      virtual string __repr__()
      {
        return "<simuPOP::selector::multiple-alleles selector>" ;
      }

    private:
      /// one locus
      vectoru m_loci;

      /// fitness for each genotype
      vectorf m_fitness;

      ///
      vectora m_wildtype;
  };

  /** \brief selection according to genotype at multiple loci multiplicative model

   multiple loci selector. This selector takes several selectors and
   multiply their fitness values...
   e.g.
     mlmselector( [mapselector(...), maselector(...) ])
   */

  class mlSelector: public selector
  {
    public:

#define SEL_None 0
#define SEL_Multiplicative 1
#define SEL_Additive 2
#define SEL_Heterogeneity 3

      typedef std::vector< Operator * > vectorop;

    public:
      /** \brief multiple loci selector using a multiplicative model.

      \param selectors a list of selectors.
      */
      mlSelector( const vectorop selectors, int mode = SEL_Multiplicative,
        int stage=PreMating, int begin=0, int end=-1, int step=1,
        vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL):
      selector(stage, begin, end, step, at, rep, grp),
        m_selectors(0), m_mode(mode)
      {
        DBG_FAILIF( selectors.empty(), ValueError, "Please specify at least one selector.");
        for(vectorop::const_iterator s = selectors.begin(), sEnd=selectors.end(); s != sEnd; ++s)
        {
          DBG_ASSERT( (*s)->__repr__().substr(10,8) == "selector", ValueError,
            "Expecting a list of fitness calculator. Given " + (*s)->__repr__());
          m_selectors.push_back( (*s)->clone());
        }
      };

      virtual ~mlSelector()
      {
        for(vectorop::iterator s = m_selectors.begin(), sEnd=m_selectors.end(); s != sEnd; ++s)
          delete *s;
      }

      virtual Operator* clone() const
      {
        throw ValueError("Multi-loci selector can not be nested.");
      }

      /// currently assuming diploid
      virtual double indFitness(individual * ind);

      virtual string __repr__()
      {
        return "<simuPOP::selector::multiple-loci selector>" ;
      }

    private:
      /// a list of selectors
      vectorop m_selectors;

      /// mode
      int m_mode;
  };

  /** \brief selection using user supplied function

  Assign fitness value by calling a user supplied function
  */

  class pySelector: public selector
  {
    public:
      /** \brief create a python hybrid selector

      \param loci susceptibility loci. The genotype at these loci will be
      passed to func.
      \param func a Python function that accept genotypes at susceptibility loci
      and return fitness value.
      \param output and other parameters please refer to help(baseOperator.__init__)
      */
      /// provide locus and fitness for 11, 12, 13 (in the form of dictionary)
      pySelector( vectoru loci, PyObject* func,
        int stage=PreMating, int begin=0, int end=-1, int step=1,
        vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL):
      selector(stage, begin, end, step, at, rep, grp),
        m_loci(loci), m_alleles(0), m_len(0), m_numArray(NULL)
      {
        if( !PyCallable_Check(func))
          throw ValueError("Passed variable is not a callable python function.");

        Py_XINCREF(func);
        m_func = func;

        DBG_FAILIF( loci.empty(), ValueError,
          "Please specify susceptibility loci");
      };

      /// destructor
      virtual ~pySelector()
      {
        if( m_func != NULL)
          Py_DECREF(m_func);
      }

      /// CPPONLY
      pySelector(const pySelector& rhs):
      selector(rhs),
        m_loci(rhs.m_loci),
        m_func(rhs.m_func),
        m_alleles(rhs.m_alleles),
        m_len(rhs.m_len),
        m_numArray(NULL)
      {
        if( m_func != NULL)
          Py_INCREF(m_func);
      }

      virtual Operator* clone() const
      {
        return new pySelector(*this);
      }

      /// currently assuming diploid
      virtual double indFitness(individual * ind);

      virtual string __repr__()
      {
        return "<simuPOP::selector::python selector>" ;
      }

    private:

      /// susceptibility loci
      vectoru m_loci;

      /// user supplied python function
      PyObject* m_func;

      /// copy of alleles of each individual a time.
      vectora m_alleles;

      /// length of m_alleles
      int m_len;

      /// the object that passed to func
      PyObject * m_numArray;

  };

  // ///////////////////////// PENETRANCE ///////////////////////////////

  /** \brief penetrance

  Please refer to the user's guide for details.
  */

  class penetrance: public Operator
  {
    public:
      /// constructor. default to be always active.
      /// default to post mating
      penetrance(bool exposePenetrance=false, int stage=DuringMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL)
        :Operator("","",stage, begin, end, step, at, rep, grp),
        m_exposePenetrance(exposePenetrance)
      {
        DBG_FAILIF( exposePenetrance==true && stage==DuringMating, ValueError,
          "Can not expose penetrance values when applied as a during mating operator.");
      }

      /// destructor
      virtual ~penetrance()
      {
      }

      /// this function is very important
      virtual Operator* clone() const
      {
        return new penetrance(*this);
      }

      /// calculate/return penetrance etc
      virtual double penet(individual * ind)
      {
        throw ValueError("This penetrance calculator is not supposed to be called directly");
        return 1.;
      }

      /// set pentrance to all individuals and record penetrance if requested.
      virtual bool apply(population& pop)
      {
        double p;
        vectorf pVec;
        size_t i=0;

        if( m_exposePenetrance )
          pVec.resize(pop.popSize());

        for (population::IndIterator it = pop.indBegin(); it != pop.indEnd(); ++it)
        {
          p = penet(&*it);
          if( m_exposePenetrance)
            pVec[i++] = p;

          if( rng().randUniform01() < p )
            it->setAffected(true);
          else
            it->setAffected(false);
        }
        if( m_exposePenetrance)
          pop.setDoubleVectorVar("penetrance", pVec);

        return true;
      }

      /// set penetrance to all individual
      virtual bool applyDuringMating(population& pop, population::IndIterator offspring,
        individual* dad=NULL, individual* mom=NULL)
      {
        if( rng().randUniform01() < penet(&*offspring) )
          offspring->setAffected(true);
        else
          offspring->setAffected(false);
        return true;
      }

      virtual string __repr__()
      {
        return "<simuPOP::penetrance>" ;
      }

    private:
      bool m_exposePenetrance;
  };

  /** \brief penetrance according to genotype at one locus

  map selector. Assign penetrance value according to a
  given dictionary.
  */

  class mapPenetrance: public penetrance
  {
    public:
      /** \brief create a map penetrance function (penetrance according to genotype at one locus

      \param locus the locus index. The genotype of this locus will be axamed.
      \param loci the loci index. The genotype of this locus will be axamed.
      \param penetrance a dictionary of penetrance. The genotype must be in the form of 'a-b' for single
         locus.
      \param phase if true, a/b and b/a will have different penetrance value. Default to false.
      \param output and other parameters please refer to help(baseOperator.__init__)
      */
      mapPenetrance( vectoru loci, const strDict& penet, bool phase=false,
        bool exposePenetrance=false, int stage=DuringMating, int begin=0, int end=-1, int step=1,
        vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL):
      penetrance(exposePenetrance, stage, begin, end, step, at, rep, grp),
        m_loci(loci), m_dict(penet), m_phase(phase)
      {
      };

      virtual ~mapPenetrance()
      {
      }

      virtual Operator* clone() const
      {
        return new mapPenetrance(*this);
      }

      /// currently assuming diploid
      virtual double penet(individual * ind)
      {
        string key;

        for(vectoru::iterator loc=m_loci.begin(); loc!=m_loci.end(); ++loc)
        {

          /// get genotype of ind
          Allele a = ind->allele(*loc, 0);
          Allele b = ind->allele(*loc, 1);

          if( loc != m_loci.begin() )
            key += '|';

          if( ! m_phase && a > b )                // ab=ba
            key +=  toStr(static_cast<int>(b)) + "-" + toStr(static_cast<int>(a));
          else
            key +=  toStr(static_cast<int>(a)) + "-" + toStr(static_cast<int>(b));
        }

        strDict::iterator pos = m_dict.find(key);

        DBG_ASSERT( pos != m_dict.end(), ValueError,
          "No penetrance value for genotype " + key);

        return( pos->second);
      }

      virtual string __repr__()
      {
        return "<simuPOP::penetrance::map penetrance>" ;
      }

    private:
      /// one locus
      vectoru m_loci;

      /// penetrance for each genotype
      strDict m_dict;

      ///
      bool m_phase;
  };

  /** \brief penetrance according to genotype at one locus

  multiple allele selector. This selector group alleles to disease
  and wild type and return penetrance to AA,Aa,aa. (A is wildtype).
  */

  class maPenetrance: public penetrance
  {
    public:
      /** \brief create a multiple allele selector (penetrance according to diseased or wildtype
      alleles)

      \param locus the locus index. The genotype of this locus will be axamed.
      \param loci the loci index.
      \param penetrance an array of penetrance of AA,Aa,aa. A is the
         wild type group. In the case of multiple loci, fitness should be in the order of
              BB Bb bb
           AA 1  2  3
           Aa 4  5  6
           aa 7  8  9
      \param wildtype an array of alleles in the wildtype group. Anything else is disease allele.,
      default = [1]
      \param output and other parameters please refer to help(baseOperator.__init__)
      */
      maPenetrance( vectoru loci, const vectorf& penet, const vectora& wildtype,
        bool exposePenetrance=false,
        int stage=DuringMating, int begin=0, int end=-1, int step=1,
        vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL):
      penetrance(exposePenetrance, stage, begin, end, step, at, rep, grp),
        m_loci(loci), m_penetrance(penet), m_wildtype(wildtype)
      {
        DBG_ASSERT( m_penetrance.size() ==  static_cast<UINT>(pow(3, loci.size())),
          ValueError, "Please specify penetrance for each combination of genotype.");
      };

      virtual ~maPenetrance()
      {
      }

      virtual Operator* clone() const
      {
        return new maPenetrance(*this);
      }

      /// currently assuming diploid
      virtual double penet(individual * ind)
      {
        UINT index = 0;
        for(vectoru::iterator loc=m_loci.begin(); loc!=m_loci.end(); ++loc)
        {

          /// get genotype of ind
          Allele a = ind->allele(*loc, 0);
          Allele b = ind->allele(*loc, 1);

          int numWildtype=0;

          // count number of wildtype
          if( find(m_wildtype.begin(), m_wildtype.end(), a) != m_wildtype.end() )
            numWildtype ++;

          if( find(m_wildtype.begin(), m_wildtype.end(), b) != m_wildtype.end() )
            numWildtype ++;
          index = index*3 + 2-numWildtype;
        }

        return m_penetrance[index];
      }

      virtual string __repr__()
      {
        return "<simuPOP::penetrance::multiple-alleles penetrance>" ;
      }

    private:
      /// one locus
      vectoru m_loci;

      /// penetrance for each genotype
      vectorf m_penetrance;

      ///
      vectora m_wildtype;
  };

  /** \brief penetrance according to genotype at multiple loci multiplicative model

   multiple loci selector. This selector takes several selectors and
   multiply their penetrance values...
   e.g.
     mlmpenetrance( [mappenetrance(...), mapenetrance(...) ])
   */

  class mlPenetrance: public penetrance
  {
    public:

#define PEN_Multiplicative 1
#define PEN_Additive 2
#define PEN_Heterogeneity 3

      typedef std::vector< Operator * > vectorop;

    public:
      /** \brief multiple loci selector using a multiplicative model.

      \param selectors a list of selectors.
      \param mode one of PEN_Multiplicative, PEN_Additive, PEN_Heterogeneity
      */
      mlPenetrance( const vectorop peneOps, int mode = PEN_Multiplicative,
        bool exposePenetrance=false, int stage=DuringMating, int begin=0, int end=-1, int step=1,
        vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL):
      penetrance(exposePenetrance, stage, begin, end, step, at, rep, grp),
        m_peneOps(0), m_mode(mode)
      {
        DBG_FAILIF( peneOps.empty(), ValueError, "Please specify at least one penetrance operator.");
        for(vectorop::const_iterator s = peneOps.begin(), sEnd=peneOps.end(); s != sEnd; ++s)
        {
          DBG_ASSERT( (*s)->__repr__().substr(10,10) == "penetrance", ValueError,
            "Expecting a list of penetrance calculator");

          m_peneOps.push_back( (*s)->clone() );
        }

      };

      virtual ~mlPenetrance()
      {
        for(vectorop::iterator s = m_peneOps.begin(), sEnd=m_peneOps.end(); s != sEnd; ++s)
          delete *s;
      }

      virtual Operator* clone() const
      {
        throw ValueError("Multi-loci selector can not be nested.");
      }

      /// currently assuming diploid
      virtual double penet(individual * ind)
      {
        if(m_mode == PEN_Multiplicative )
        {
          // x1 x2 x3 ...
          double pen = 1;
          for(vectorop::iterator s = m_peneOps.begin(), sEnd=m_peneOps.end();
            s != sEnd; ++s)
          pen *= static_cast<penetrance *>(*s)->penet(ind);
          return pen;
        }
        else if(m_mode == PEN_Additive )
        {
          // x1 + x2 + x3
          double pen = 0;
          for(vectorop::iterator s = m_peneOps.begin(), sEnd=m_peneOps.end();
            s != sEnd; ++s)
          pen +=  static_cast<penetrance *>(*s)->penet(ind);
          return pen>1?1:pen;
        }
        else if(m_mode == PEN_Heterogeneity )
        {
          // 1-(1-x1)(1-x2)
          double pen = 1;
          for(vectorop::iterator s = m_peneOps.begin(), sEnd=m_peneOps.end();
            s != sEnd; ++s)
          pen *= 1 - static_cast<penetrance *>(*s)->penet(ind);
          return 1 - pen;
        }

        return 0.0;
      }

      virtual string __repr__()
      {
        return "<simuPOP::penetrance::multiple-loci penetrance>" ;
      }

    private:
      /// a list of peneOps
      vectorop m_peneOps;

      /// mode
      int m_mode;
  };

  /** \brief penetrance using user supplied function

  Assign penetrance value by calling a user supplied function
  */

  class pyPenetrance: public penetrance
  {
    public:
      /** \brief create a python hybrid selector

      \param loci susceptibility loci. The genotype at these loci will be
      passed to func.
      \param func a Python function that accept genotypes at susceptibility loci
      and return penetrance value.
      \param output and other parameters please refer to help(baseOperator.__init__)
      */
      /// provide locus and penetrance for 11, 12, 13 (in the form of dictionary)
      pyPenetrance( vectoru loci, PyObject* func, bool exposePenetrance=false,
        int stage=DuringMating, int begin=0, int end=-1, int step=1,
        vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL):
      penetrance(exposePenetrance, stage, begin, end, step, at, rep, grp),
        m_loci(loci), m_alleles(0), m_len(0), m_numArray(NULL)
      {
        if( !PyCallable_Check(func))
          throw ValueError("Passed variable is not a callable python function.");

        Py_XINCREF(func);
        m_func = func;

        DBG_FAILIF( loci.empty(), ValueError,
          "Please specify susceptibility loci");
      };

      /// destructor
      virtual ~pyPenetrance()
      {
        if( m_func != NULL)
          Py_DECREF(m_func);
        if( m_numArray != NULL)
          Py_DECREF(m_numArray);
      }

      /// CPPONLY
      pyPenetrance(const pyPenetrance& rhs):
      penetrance(rhs),
        m_loci(rhs.m_loci),
        m_func(rhs.m_func),
        m_alleles(rhs.m_alleles),
        m_len(rhs.m_len),
        m_numArray(NULL)
      {
        if( m_func != NULL)
          Py_INCREF(m_func);
      }

      virtual Operator* clone() const
      {
        return new pyPenetrance(*this);
      }

      /// currently assuming diploid
      virtual double penet(individual * ind)
      {
        int len = m_loci.size() * ind->ploidy();
        if( m_len != len )
        {
          m_len = len;
          m_alleles.resize(m_len);
          if(m_numArray != NULL)
            Py_DECREF(m_numArray);
          m_numArray = Allele_Vec_As_NumArray( m_alleles.begin(), m_alleles.end() );
        }

        UINT pEnd = ind->ploidy();
        for(size_t i=0, iEnd=m_loci.size(), j=0; i < iEnd; ++i)
          for(UINT p=0; p < pEnd; ++p)
            m_alleles[j++] = ind->allele(m_loci[i], p);

        double resDouble;
        PyCallFunc(m_func, "(O)", m_numArray, resDouble, PyObj_As_Double);

        // make sure the returned value is legitimate.
        DBG_ASSERT( fcmp_ge( resDouble, 0.) && fcmp_le( resDouble, 1.),
          ValueError, "Returned fitness " + toStr(resDouble) + " is out of range [0,1]" );

        return resDouble;
      }

      virtual string __repr__()
      {
        return "<simuPOP::penetrance::python penetrance>" ;
      }

    private:

      /// susceptibility loci
      vectoru m_loci;

      /// user supplied python function
      PyObject* m_func;

      /// copy of alleles of each individual a time.
      vectora m_alleles;

      /// length of m_alleles
      int m_len;

      /// the object that passed to func
      PyObject * m_numArray;

  };

  // ///////////////////////// Quantitative trait ////////////////////////

  /** \brief quantitative trait

  Genetic quantitative trait is tricky to simulate. In simuPOP, I employee
  an ability (fitness) to mate approach. Namely, the probability
  that an individual will be chosen for mating is proportional
  to its fitness value. More specifically,

  - PreMating selectors assign fitness values to each individual.

  - Sexless mating (e.g. binomialSelection) : individuals are chosen
  at probabilities that are proportional to their fitness values. More
  specifically, if there are N individuals with fitness values
  \f$f_i, i=1,...,N \f$, individual \f$i\f$ will have probability
  \f$ \frac{f_i}{\sum_{j=1}^N f_j} \f$ to be chosen to be passed
  to the next generation.

  - Random mating with sex (e.g. randommating): males and females are
  separated and each are chosen as described above.

  Please refer to the user's guide for details.
  */

  class quanTrait: public Operator
  {
    public:
      /// constructor. default to be always active.
      quanTrait( int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL)
        :Operator("","",stage, begin, end, step, at, rep, grp),
        m_qtrait(0), m_len(0)
      {
      }

      /// destructor
      virtual ~quanTrait()
      {
      }

      virtual Operator* clone() const
      {
        return new quanTrait(*this);
      }

      /// calculate/return quantitative trait etc
      virtual double qtrait(individual * ind)
      {
        ///
        throw ValueError("This quantitative trait calculator is not supposed to be called directly");
        return 1.;
      }

      /// set qtrait to all individual
      bool apply(population& pop)
      {
        if(static_cast<size_t>(m_len) != pop.popSize() )
        {
          m_len = pop.popSize();
          m_qtrait.resize(m_len);
          // fill(m_qtrait.begin(), m_qtrait.end(), 0.0);
        }

        size_t i=0;
        for (population::IndIterator it = pop.indBegin(); it != pop.indEnd(); ++it)
          m_qtrait[i++] = qtrait(&*it) ;

        pop.setDoubleVectorVar("qtrait", m_qtrait);
        return true;
      }

      virtual string __repr__()
      {
        return "<simuPOP::qtrait::quantitative trait>" ;
      }

    private:

      /// qtrait values
      vectorf m_qtrait;

      /// length of qtrait values
      int m_len;

  };

  /** \brief quantitative trait according to genotype at one locus

  map selector. Assign qtrait value according to a
  given dictionary.
  */

  class mapQuanTrait: public quanTrait
  {
    public:
      /** \brief create a map selector (quantitative trait according to genotype at one locus

      \param locus the locus index. The genotype of this locus will be axamed.
      \param loci the loci.
      \param qtrait a dictionary of qtrait. The genotype must be in the form of 'a-b'. This is the mean
      of quantitative trait. The actual trait value will be N(mean, sigma^2)
      For multiple loci, the form is 'a-b|c-d|e-f' etc.
      \param sigma standard deviation of the environmental facotr N(0,sigma^2).
      \param phase if true, a/b and b/a will have different qtrait value. Default to false.
      \param output and other parameters please refer to help(baseOperator.__init__)
      */
      mapQuanTrait( vectoru loci, const strDict& qtrait, double sigma=0, bool phase=false,
        int stage=PostMating, int begin=0, int end=-1, int step=1,
        vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL):
      quanTrait(stage, begin, end, step, at, rep, grp),
        m_loci(loci), m_dict(qtrait), m_sigma(sigma), m_phase(phase)
      {
      };

      virtual ~mapQuanTrait()
      {
      }

      virtual Operator* clone() const
      {
        return new mapQuanTrait(*this);
      }

      /// currently assuming diploid
      virtual double qtrait(individual * ind)
      {
        string key;

        for(vectoru::iterator loc=m_loci.begin(); loc!=m_loci.end(); ++loc)
        {
          /// get genotype of ind
          Allele a = ind->allele(*loc, 0);
          Allele b = ind->allele(*loc, 1);

          if( loc != m_loci.begin() )
            key += '|';

          if( ! m_phase && a > b )                // ab=ba
            key +=  toStr(static_cast<int>(b)) + "-" + toStr(static_cast<int>(a));
          else
            key +=  toStr(static_cast<int>(a)) + "-" + toStr(static_cast<int>(b));
        }
        strDict::iterator pos = m_dict.find(key);

        DBG_ASSERT( pos != m_dict.end(), ValueError,
          "No qtrait value for genotype " + key);

        return( rng().randNormal(pos->second, m_sigma) );
      }

      virtual string __repr__()
      {
        return "<simuPOP::qtrait::map quantitative trait>" ;
      }

    private:
      /// one locus
      vectoru m_loci;

      /// qtrait for each genotype
      strDict m_dict;

      ///
      double m_sigma;

      ///
      bool m_phase;
  };

  /** \brief quantitative trait according to genotype at one locus

  multiple allele selector. This selector group alleles to disease
  and wild type and return qtrait to AA,Aa,aa. (A is wildtype).
  */

  class maQuanTrait: public quanTrait
  {
    public:
      /** \brief create a multiple allele selector (quantitative trait according to diseased or wildtype
      alleles)

      \param locus the locus index. The genotype of this locus will be axamed.
      \param qtrait an array of qtrait of AA,Aa,aa. A is the wild type group.
      \param sigma an array of standard deviation for each of the trait genotype (AA, Aa, aa)
      \param wildtype an array of alleles in the wildtype group. Anything else is disease allele.
         default = [1]
      \param output and other parameters please refer to help(baseOperator.__init__)
      */
      maQuanTrait( vectoru loci, const vectorf& qtrait, const vectora& wildtype,
        const vectorf& sigma = vectorf(),
        int stage=PostMating, int begin=0, int end=-1, int step=1,
        vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL):
      quanTrait(stage, begin, end, step, at, rep, grp),
        m_loci(loci), m_qtrait(qtrait), m_sigma(sigma), m_wildtype(wildtype)
      {
        if( m_sigma.empty())
          m_sigma.resize(3,0.);

        DBG_ASSERT( m_qtrait.size() == static_cast<UINT>(pow(3, loci.size()))
          && m_sigma.size() == m_qtrait.size(),
          ValueError, "Please specify qtrait for every combination of genotype.");
      };

      /// destructor
      virtual ~maQuanTrait()
      {
      }

      virtual Operator* clone() const
      {
        return new maQuanTrait(*this);
      }

      /// currently assuming diploid
      virtual double qtrait(individual * ind)
      {
        UINT index = 0;
        for(vectoru::iterator loc=m_loci.begin(); loc!=m_loci.end(); ++loc)
        {
          /// get genotype of ind
          Allele a = ind->allele(*loc, 0);
          Allele b = ind->allele(*loc, 1);

          int numWildtype=0;

          // count number of wildtype
          if( find(m_wildtype.begin(), m_wildtype.end(), a) != m_wildtype.end() )
            numWildtype ++;

          if( find(m_wildtype.begin(), m_wildtype.end(), b) != m_wildtype.end() )
            numWildtype ++;
          index = index*3 + 2-numWildtype;
        }

        return rng().randNormal(m_qtrait[index], m_sigma[index]);
      }

      virtual string __repr__()
      {
        return "<simuPOP::qtrait::multiple-alleles qtrait>" ;
      }

    private:
      /// one locus
      vectoru m_loci;

      /// qtrait for each genotype
      vectorf m_qtrait;

      ///
      vectorf m_sigma;

      ///
      vectora m_wildtype;
  };

  /** \brief quantitative trait according to genotype at multiple loci multiplicative model

   multiple loci selector. This selector takes several selectors and
   multiply their qtrait values...
   e.g.
     mlmquanTrait( [mapquanTrait(...), maquanTrait(...) ])
   */

  class mlQuanTrait: public quanTrait
  {
    public:

#define QT_Multiplicative 1
#define QT_Additive 2

      /// vector of operator pointers.
      typedef std::vector< Operator * > vectorop;

    public:
      /** \brief multiple loci selector using a multiplicative model.

      \param qtraits a list of qtraits.
      */
      mlQuanTrait( const vectorop qtraits, int mode = QT_Multiplicative, double sigma=0,
        int stage=PostMating, int begin=0, int end=-1, int step=1,
        vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL):
      quanTrait(stage, begin, end, step, at, rep, grp),
        m_qtraits(0), m_sigma(sigma), m_mode(mode)
      {
        DBG_FAILIF( qtraits.empty(), ValueError, "Please specify at least one selector.");
        for(vectorop::const_iterator s = qtraits.begin(), sEnd=qtraits.end(); s != sEnd; ++s)
        {
          DBG_ASSERT( (*s)->__repr__().substr(10,6) == "qtrait", ValueError,
            "Expecting a vector of quantitative trait calculator");
          m_qtraits.push_back( (*s)->clone() );
        }

      };

      virtual ~mlQuanTrait()
      {
        for(vectorop::iterator s = m_qtraits.begin(), sEnd=m_qtraits.end(); s != sEnd; ++s)
          delete *s;
      }

      virtual Operator* clone() const
      {
        throw ValueError("Multi-loci selector can not be nested.");
      }

      /// currently assuming diploid
      virtual double qtrait(individual * ind)
      {
        if(m_mode == QT_Multiplicative )
        {
          double fit = 1;
          for(vectorop::iterator s = m_qtraits.begin(), sEnd=m_qtraits.end();
            s != sEnd; ++s)
          fit *= static_cast<quanTrait *>(*s)->qtrait(ind);
          return rng().randNormal(fit, m_sigma);
        }
        else if(m_mode == QT_Additive )
        {
          double fit = 0;
          for(vectorop::iterator s = m_qtraits.begin(), sEnd=m_qtraits.end();
            s != sEnd; ++s)
          fit +=  static_cast<quanTrait *>(*s)->qtrait(ind);
          return rng().randNormal(fit, m_sigma);
        }
        return 0.;
      }

      virtual string __repr__()
      {
        return "<simuPOP::qtrait::multiple-loci qtrait>" ;
      }

    private:
      /// a list of qtraits
      vectorop m_qtraits;

      ///
      double m_sigma;

      /// mode
      int m_mode;
  };

  /** \brief quantitative trait using user supplied function

  Assign qtrait value by calling a user supplied function
  */

  class pyQuanTrait: public quanTrait
  {
    public:
      /** \brief create a python hybrid selector

      \param loci susceptibility loci. The genotype at these loci will be
      passed to func.
      \param func a Python function that accept genotypes at susceptibility loci
      and return qtrait value.
      \param output and other parameters please refer to help(baseOperator.__init__)
      */
      /// provide locus and qtrait for 11, 12, 13 (in the form of dictionary)
      pyQuanTrait( vectoru loci, PyObject* func,
        int stage=PostMating, int begin=0, int end=-1, int step=1,
        vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL):
      quanTrait(stage, begin, end, step, at, rep, grp),
        m_loci(loci), m_alleles(0), m_len(0), m_numArray(NULL)
      {
        if( !PyCallable_Check(func))
          throw ValueError("Passed variable is not a callable python function.");

        Py_XINCREF(func);
        m_func = func;

        DBG_FAILIF( loci.empty(), ValueError,
          "Please specify susceptibility loci");
      };

      virtual ~pyQuanTrait()
      {
        if( m_func != NULL)
          Py_DECREF(m_func);
      }

      /// CPPONLY
      pyQuanTrait(const pyQuanTrait& rhs):
      quanTrait(rhs),
        m_loci(rhs.m_loci),
        m_func(rhs.m_func),
        m_alleles(rhs.m_alleles),
        m_len(rhs.m_len),
        m_numArray(NULL)
      {
        if( m_func != NULL)
          Py_INCREF(m_func);
      }

      virtual Operator* clone() const
      {
        return new pyQuanTrait(*this);
      }

      /// currently assuming diploid
      virtual double qtrait(individual * ind)
      {
        if( m_len == 0)
        {
          m_len = m_loci.size() * ind->ploidy();
          m_alleles.resize( m_len);
          m_numArray = Allele_Vec_As_NumArray( m_alleles.begin(), m_alleles.end() );
        }

        DBG_FAILIF( static_cast<size_t>(m_len) != ind->ploidy() * m_loci.size(),
          SystemError,
          "Length of m_len is wrong. Have you changed pop type?" );

        UINT pEnd = ind->ploidy();
        for(size_t i=0, iEnd=m_loci.size(), j=0; i < iEnd; ++i)
          for(UINT p=0; p < pEnd; ++p)
            m_alleles[j++] = ind->allele(m_loci[i], p);

        double resDouble;
        PyCallFunc(m_func, "(O)", m_numArray, resDouble, PyObj_As_Double);
        return resDouble;
      }

      virtual string __repr__()
      {
        return "<simuPOP::qtrait::python qtrait>" ;
      }

    private:

      /// susceptibility loci
      vectoru m_loci;

      /// user supplied python function
      PyObject* m_func;

      /// copy of alleles of each individual a time.
      vectora m_alleles;

      /// length of m_alleles
      int m_len;

      /// the object that passed to func
      PyObject * m_numArray;

  };

  // ///////////////////// SUBSET ///////////////////////////////////////////
  // // ascertainment ........................................

  /// thrink population accroding to some outside value

  class pySubset: public Operator
  {

    public:
      /// create a directmigrator
      /**
      \param keep a carray of the length of population.
      its values will be assigned to info.
      \stage and other parameters please see help(baseOperator.__init__)
      */
      pySubset( PyObject* keep=NULL,
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL)
        : Operator( "", "", stage, begin, end, step, at, rep, grp)
      {
        DBG_ASSERT( PyObj_Is_IntNumArray(keep), ValueError,
          "Passed vector is not a Python/Numeric int array");
        Py_INCREF(keep);
        m_keep = keep;
      }

      /// destructor
      virtual ~pySubset()
      {
        if( m_keep != NULL)
          Py_DECREF(m_keep);
      }

      /// CPPONLY
      pySubset(const pySubset& rhs):
      Operator(rhs),
        m_keep(rhs.m_keep)
      {
        if( m_keep != NULL)
          Py_INCREF(m_keep);
      }

      /// this function is very important
      virtual Operator* clone() const
      {
        return new pySubset(*this);
      }

      virtual bool apply(population& pop)
      {

        DBG_ASSERT( NumArray_Size(m_keep) >= static_cast<int>(pop.popSize()) ,
          ValueError, "Given subpopid array has a length of "
          + toStr( NumArray_Size(m_keep)) + " which is less than population size "
          + toStr(pop.popSize()));

        long * id = reinterpret_cast<long*>(NumArray_Data(m_keep));

        for(size_t i=0, iEnd=pop.popSize(); i<iEnd; ++i)
          pop.ind(i).setInfo( id[i] );

        pop.setSubPopByIndInfo();
        return true;
      }

      virtual string __repr__()
      {
        return "<simuPOP::pySubset>" ;
      }

    private:
      PyObject* m_keep;

  };

  /// sample from population and save samples

  /// sample operator will generate a new subpopulation in pop namespace.

  class sample: public Operator
  {

    public:
      /// create a sample
      /**
      \param name variable name of the sampled population (will be put in pop local namespace)
      \param nameExpr expression version of name. If both name and nameExpr is empty, do not store pop.
      \param times how many times to run the sample process? This is usually one, but we
         may want to take several random samples.
      \param saveAs filename to save the population.
      \param saveAsExpr expression for save filename
      \param format to save sample(s)
      \param stage and other parameters please see help(baseOperator.__init__)
      */
      sample( const string& name="sample", const string& nameExpr="", UINT times=1,
        const string& saveAs="", const string& saveAsExpr="",   const string& format="auto",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL)
        : Operator( "", "", stage, begin, end, step, at, rep, grp),
        m_name(name), m_nameExpr(nameExpr,""), m_times(times), m_saveAs(saveAs),
        m_saveAsExpr(saveAsExpr), m_format(format)
      {
        if(name=="" && nameExpr=="" && saveAs=="" && saveAsExpr=="")
          throw ValueError("Please specify name (or nameExpr) or saveAs (or saveAsExpr) to save sample.");
      }

      /// destructor
      virtual ~sample(){};

      /// this function is very important
      virtual Operator* clone() const
      {
        return new sample(*this);
      }

      virtual bool prepareSample(population& )
      {
        return true;
      }

      virtual population& drawsample(population& pop)
      {
        // keep pop untouched.
        // use ind.info() directly
        throw SystemError("This function is not supposed to be called directly");
        return pop;
      }

      /// return the samples
      PyObject* samples(population& pop);

      virtual bool apply(population& pop);

      virtual string __repr__()
      {
        return "<simuPOP::sample>" ;
      }

    private:

      /// name to save sample, default to 'sample'
      string m_name;

      /// pop name
      Expression m_nameExpr;

      /// sample times
      UINT m_times;

      /// filename to save sample
      string  m_saveAs;

      /// saveas expression
      Expression m_saveAsExpr;

      /// format to save samples
      string m_format;
  };

  /// thrink population accroding to some outside value

  class randomSample: public sample
  {

    public:
      /// draw random sample, regardless of affected status
      /**
      \param size size of sample. It can be either a number, representing
        the overall sample size, regardless of population strucutre;
        or an array, representing number of samples drawn from each subpopulation.
      \param stage and other parameters please see help(baseOperator.__init__)
      \param name variable name of the sampled population (will be put in pop local namespace)
      \param nameExpr expression version of name. If both name and nameExpr is empty, do not store pop.
      \param times how many times to run the sample process? This is usually one, but we
         may want to take several random samples.
      \param saveAs filename to save the population.
      \param saveAsExpr expression for save filename
      \param format to save sample(s)
      \param stage and other parameters please see help(baseOperator.__init__)

      \note ancestral populations will not be copied to the samples
      */
      randomSample( vectorlu size=vectorlu(),
        const string& name="sample", const string& nameExpr="", UINT times=1,
        const string& saveAs="", const string& saveAsExpr="",   const string& format="auto",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL)
        : sample(name, nameExpr, times, saveAs, saveAsExpr, format,
        stage, begin, end, step, at, rep, grp),
        m_size(size)
      {
      }

      /// destructor
      virtual ~randomSample(){};

      /// this function is very important
      virtual Operator* clone() const
      {
        return new randomSample(*this);
      }

      /// value checking
      virtual bool prepareSample(population& pop);

      virtual population& drawsample(population& pop);

      virtual string __repr__()
      {
        return "<simuPOP::random sample>" ;
      }

    private:
      /// sample size
      vectorlu m_size;
  };

  /// thrink population accroding to some outside value

  class caseControlSample: public sample
  {

    public:
      /// draw cases and controls
      /**
      \param cases number of cases, or an array of number of cases from
        each subpopulation.
      \param controls number of controls, or an array of number of controls
        from each subpopulation.
      \param name variable name of the sampled population (will be put in pop local namespace)
      \param nameExpr expression version of name. If both name and nameExpr is empty, do not store pop.
      \param times how many times to run the sample process? This is usually one, but we
         may want to take several random samples.
      \param saveAs filename to save the population.
      \param saveAsExpr expression for save filename
      \param format to save sample(s)
      \param stage and other parameters please see help(baseOperator.__init__)
      */
      caseControlSample( const vectori& cases=vectori(), const vectori& controls = vectori(),
        bool spSample=false, const string& name="sample", const string& nameExpr="", UINT times=1,
        const string& saveAs="", const string& saveAsExpr="",   const string& format="auto",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL)
        : sample(name, nameExpr, times, saveAs, saveAsExpr, format,
        stage, begin, end, step, at, rep, grp),
        m_numCases(cases), m_numControls(controls), m_spSample(spSample),
        m_caseIdx(0), m_controlIdx(0)
      {
      }

      /// destructor
      virtual ~caseControlSample(){};

      /// this function is very important
      virtual Operator* clone() const
      {
        return new caseControlSample(*this);
      }

      virtual bool prepareSample(population& pop );

      virtual population& drawsample(population& pop);

      virtual string __repr__()
      {
        return "<simuPOP::case control sample>" ;
      }

    private:

      /// number of cases, use vectori instead of vectorlu because
      /// this will be post to setIntVectorVar
      vectori m_numCases, m_numControls;

      /// whether or not sample from each subpop
      bool m_spSample;

      vector< vectori > m_caseIdx, m_controlIdx;
  };

  /// thrink population accroding to some outside value

  class affectedSibpairSample: public sample
  {

    public:
      /// draw cases and controls
      /**
      \param size number of affected sibpairs to be sampled.
        Can be a number or an array. If a number is given, it is
        the total number of sibpairs, ignoring population structure.
        Otherwise, given number of sibpairs are sampled from
        subpopulations. If size is unspecified, this operator
        will return all affected sibpairs.
      \param countOnly set variables about number of affected sibpairs,
      do not actually draw the sample
      \param name variable name of the sampled population (will be put in pop local namespace)
      \param nameExpr expression version of name. If both name and nameExpr is empty, do not store pop.
      \param times how many times to run the sample process? This is usually one, but we
      may want to take several random samples.
      \param saveAs filename to save the population.
      \param saveAsExpr expression for save filename
      \param format to save sample(s)
      \param stage and other parameters please see help(baseOperator.__init__)
      */
      affectedSibpairSample( vectoru size = vectoru(),
        bool chooseUnaffected=false,
        bool countOnly=false,
        const string& name="sample", const string& nameExpr="", UINT times=1,
        const string& saveAs="", const string& saveAsExpr="",
        const string& format="auto",
        int stage=PostMating, int begin=0, int end=-1,
        int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL)
        : sample(name, nameExpr, times, saveAs, saveAsExpr, format,
        stage, begin, end, step, at, rep, grp),
        m_size(size), m_affectedness(!chooseUnaffected), m_countOnly(countOnly),
        m_sibpairs(0)
      {
      }

      /// destructor
      virtual ~affectedSibpairSample(){};

      /// this function is very important
      virtual Operator* clone() const
      {
        return new affectedSibpairSample(*this);
      }

      virtual bool prepareSample(population& pop);

      virtual population& drawsample(population& pop);

      virtual string __repr__()
      {
        return "<simuPOP::affected sibpair sample>" ;
      }

    private:
      /// sample size
      vectoru m_size;

      bool m_affectedness;

      // do not draw sample
      bool m_countOnly;

      /// sibs for all subpopulations
      vector< vector< std::pair<ULONG, ULONG> > > m_sibpairs;

  };

  /// thrink population accroding to some outside value

  class pySample: public sample
  {

    public:
      /// create a python sampler
      /**
      \param keep a carray of the length of population.
      its values will be assigned to info.
      \stage and other parameters please see help(baseOperator.__init__)
      */
      pySample( PyObject * keep, bool keepAncestralPops,
        const string& name="sample", const string& nameExpr="", UINT times=1,
        const string& saveAs="", const string& saveAsExpr="",   const string& format="auto",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL)
        : sample(name, nameExpr, times, saveAs, saveAsExpr, format,
        stage, begin, end, step, at, rep, grp),
        m_keepAncestralPops(keepAncestralPops)
      {
        DBG_ASSERT( PyObj_Is_IntNumArray(keep), ValueError,
          "Passed vector is not a Python/Numeric int array");
        Py_INCREF(keep);
        m_keep = keep;
      }

      /// destructor
      virtual ~pySample()
      {
        if( m_keep != NULL)
          Py_DECREF(m_keep);
      }

      /// CPPONLY
      pySample(const pySample& rhs):
      sample(rhs),
        m_keep(rhs.m_keep)
      {
        if( m_keep != NULL)
          Py_INCREF(m_keep);
      }

      /// this function is very important
      virtual Operator* clone() const
      {
        return new pySample(*this);
      }

      virtual population& drawsample(population& pop)
      {
        DBG_ASSERT( NumArray_Size(m_keep) >= static_cast<int>(pop.popSize()) ,
          ValueError, "Given subpopid array has a length of "
          + toStr( NumArray_Size(m_keep)) + " which is less than population size "
          + toStr(pop.popSize()));

        long * id = reinterpret_cast<long*>(NumArray_Data(m_keep));

        for(size_t i=0, iEnd=pop.popSize(); i<iEnd; ++i)
          pop.ind(i).setInfo( id[i] );

        return pop.newPopByIndInfo(m_keepAncestralPops);
      }

      virtual string __repr__()
      {
        return "<simuPOP::pySubset>" ;
      }

    private:
      PyObject* m_keep;

      bool m_keepAncestralPops;
  };

}
#endif
