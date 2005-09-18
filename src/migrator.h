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

#ifndef _MIGRATOR_H
#define _MIGRATOR_H
/**
\file
\brief head file of class Migrator:public Operator
*/
#include "operator.h"
#include <list>
using std::list;

#include <iostream>
using std::cout;
using std::endl;

#include <algorithm>
using std::sort;
using std::random_shuffle;

#include "simupop_cfg.h"
#include <math.h>
// using std::sqrt;

namespace simuPOP
{
  /**
  @author Bo Peng
  */
  template<class Pop>
    class Migrator: public Operator<Pop>
  {

    public:

#define MigrByProbability  1
#define MigrByProportion   2
#define MigrByCounts       3

    public:
      /// create a migrator
      /**
      \param rate migration rate, proportion or count. Determined by parameter mode.
         rate should be a m by n matrix. If a number is given,
         the migration rate will be r*ones(m,n).
      \param mode one of MigrByProbability (default), MigrByProportion or MigrByCounts
      \param fromSubPop an array of 'from' subpops, default to all subpopulations.
        If a single subpop is specified, [] can be ignored. I.e., [a] is equvalent to a.
      \param toSubPop an array of 'to' subpops, default to all subpopulations. If
        a single subpop is specified, [] can be ignored.
      \param stage is default to PreMating. For details about other parameters,
      please refer to help(baseOperator.__init__)

      rate is a matrix with dimensions determined by fromSubPop and toSubPop. By default,
      rate is a matrix with element (i,j) being the migration rate, probability or count
      from subpop i to subpop j. If fromSubPop and/or toSubPop are given, migration only
      happen between these subpopulations. An extreme case is 'point migration'

      rate=[[r]], fromSubPop=a, toSubPop=b

      which migrate from subpop a to b with given rate r.
      */
      Migrator( const matrix& rate, int mode = MigrByProbability,
        vectoru fromSubPop=vectoru(), vectoru toSubPop=vectoru(),
        int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL, string sep="\t")
        : Operator<Pop>( "", "", stage, begin, end, step, at, rep, grp, sep),
        m_rate(0), m_mode(0), m_from(fromSubPop), m_to(toSubPop)
      {

        DBG_FAILIF( rate.size() > 1 && !m_from.empty() && m_from.size() != rate.size(),
          ValueError, "Length of param from must match rows of rate matrix.");

        DBG_FAILIF( rate.size() > 1 && !m_to.empty() && m_to.size() != rate[0].size(),
          ValueError, "Length of param to must match columns of rate matrix.");

        DBG_FAILIF( rate.size() == 1 && (m_from.size() != 1 || m_to.size() != 1),
          ValueError, "If rate=r (assumed to be [[r]], fromSubPOp and toSubPop should not be empty()");

        setRates(rate, mode);
      };

      /// destructor
      virtual ~Migrator(){};

      virtual Operator<Pop>* clone() const
      {
        return new Migrator<Pop>(*this);
      }

      /// return rate
      matrix rate()
      {
        return m_rate;
      }

      /// set migration rate
      /** format 0-0 0-1 0-2, 1-0 1-1 1-2, 2-0, 2-1, 2-2.
      for mode 1 or 2, 00,11,22 will be set automatically.
      regardless of input.
      */
      void setRates(const matrix& rate, int mode)
      {
        if( rate.empty() )
        {
          m_from.clear();
          m_to.clear();
          return;
        }

        UINT szFrom = rate.size();
        UINT szTo = rate[0].size();

        if( m_from.empty() )
          for(UINT i=0; i<szFrom; ++i)
            m_from.push_back(i);

        if( m_to.empty() )
          for(UINT i=0; i<szTo; ++i)
            m_to.push_back(i);

        m_mode = mode;

        if( m_mode != MigrByProbability && m_mode !=MigrByProportion
          && m_mode != MigrByCounts)
          throw ValueError("Migration mode can only be MigrByProbability), "
            " MigrByProportion or MigrByCounts");

        // check parameters
        for(UINT i=0; i < szFrom; ++i)
        {
          if( rate[i].size() != szTo )
            throw ValueError("Expecting a matrix of migration rate.");

          for(size_t j=0; j < szTo; ++j )
            if( fcmp_lt( rate[i][j], 0.)
            || (m_mode != MigrByCounts &&  fcmp_gt(rate[i][j], 1.)))
              throw ValueError("Migration rate should be in the range of [0,1]");
        }

        m_rate = rate;

        /// set r[i][i]--- may need to extend rate (to add i->i)
        if( m_mode == MigrByProbability || m_mode == MigrByProportion)
        {
          for(UINT i=0; i< szFrom; i++)           // from
          {
            // look for from=to cell.
            UINT spFrom = m_from[i];
            double sum = accumulate(m_rate[i].begin(), m_rate[i].end(), 0.0);
            //
            vectoru::iterator spTo = find(m_to.begin(), m_to.end(), spFrom);
            if( spTo == m_to.end() )              // if no to, only check if sum <= 1
            {
              if( fcmp_gt(sum, 1.0) )
                throw ValueError("Sum of migrate rate from one subPop should <= 1");
              // adding i->i item
              m_rate[i].push_back(1.0 - sum);
            }
            else                                  // if not, set r[i][i]
            {
              double & self = m_rate[i][ spTo - m_to.begin() ];
              sum -= self;
              if( fcmp_gt(sum, 1.0) )
                throw ValueError("Sum of migrate rate from one subPop should <= 1");
              // reset to-my-self probability/proportion
              self = 1.0 - sum;
            }
          }
        }

        if( m_mode == MigrByProbability )         // probability, we need cumulative probability
        {
          for(size_t i=0; i < szFrom; i++)        // from
          {
            size_t jEnd = m_rate[i].size();
            for(size_t j=1; j < jEnd; j++)
              // rate might have been extended to save i->i info
              // the last one should be 1
              m_rate[i][j] += m_rate[i][j-1];

            DBG_FAILIF( fcmp_ne(m_rate[i][ jEnd-1] , 1.0) , ValueError,
              "Accumulative rate should have 1 at last.");
          }
        }
      }

      virtual bool apply(Pop& pop)
      {
        // set info of individual
        pop.setIndInfoWithSubPopID();

        typename Pop::IndIterator ind, indEd;

        vectorlu toIndices(0);

        for(UINT from=0, fromEnd=m_from.size(); from < fromEnd; ++from)
        {
          UINT spFrom = m_from[from];
          // rateSize might be toSize + 1, the last one is from->from
          UINT toSize = m_to.size(), toIndex;
          WeightedSampler ws(rng(), m_rate[from]);

          // m_from out of range.... ignore.
          if( spFrom >= pop.numSubPop() )
            continue;

          if(m_mode == MigrByProbability )        // migrate by probability
          {
            // for each individual, migrate according to migration probability
            for(ind=pop.indBegin( spFrom), indEd = pop.indEnd(spFrom);
              ind != indEd;  ++ind)
            {
              //toIndex = rng().randIntByFreq( rateSize, &m_rate[from][0] ) ;
              toIndex = ws.get();

              DBG_ASSERT( toIndex < m_rate[from].size(), ValueError,
                "Return index out of range.");

              // rateSize==toSize (no i->i addition)
              //   toIndex < toSize
              // rateSize = toSize + 1, ignore i->1 (last one)
              //  toIndex < toSize
              if( toIndex < toSize && m_to[toIndex] != spFrom )
                ind->setInfo( m_to[toIndex] );
            }
            continue;
          }

          /// 2nd, or 3rd method
          // first find out how many people will move to other subPop
          // then randomly assign individuals to move
          vectorlu toNum(toSize);
          ULONG spSize = pop.subPopSize(spFrom);
          if( m_mode == MigrByProportion )
          {
            for(UINT i=0; i<toSize; ++i)
              toNum[i] = static_cast<ULONG>(spSize * m_rate[from][i]);
          }
          else                                    // by count
          {
            for(UINT i=0; i<toSize; ++i)
              toNum[i] = static_cast<ULONG>(m_rate[from][i]);
          }
          /// create a vector and assign indices, then random shuffle
          /// and assign info
          toIndices.resize( spSize);
          UINT k=0;
          for(UINT i=0; i<toSize && k<spSize; ++i)
            for(UINT j=0; j< toNum[i] && k<spSize; ++j)
              toIndices[k++] = m_to[i];

          while( k < spSize)
            toIndices[k++] = spFrom;

          random_shuffle(toIndices.begin(), toIndices.end());
          ind = pop.indBegin(spFrom);
          // set info
          for(UINT i=0; i<spSize; ++i)
            (ind + i )->setInfo( toIndices[i] );
        }                                         /// for all subPop.

        // do migration.
        // true: rearrange individuals
        pop.setSubPopByIndInfo();

        return true;
      }

      virtual string __repr__()
      {
        return "<simuPOP::migrator>" ;
      }

    private:

      /// migration rate. its meaning is controled by m_mode
      matrix m_rate;

      /// asProbability (1), asProportion (2), or asCounts.
      int m_mode;

      /// from->to subPop index.
      /// default to 0 - rows of rate - 1, 0 - columns of rate - 1
      vectoru m_from, m_to;
  };

  /** \brief migrate using given info vector

  You can use directMigrator to accomplish any migration: that is to say
  you directly specify subpopulation numbers for each individual and
  this operator will do the rest.

  */
  template<class Pop>
    class PyMigrator: public Operator<Pop>
  {

    public:
      /// create a directMigrator
      /**
      This operator accept a one-dimensional Numeric Python int array. (created by Numeric.array ).
      The contend of the array will be considered as subpopulation id.

      \param subPopID a 1-d Numeric Python int array. Must has length greater or equal to
      population size.
      \param stage is default to PreMating, please refer to help(baseOperator.__init__)
      for details about other parameters.
      */
      PyMigrator( PyObject* subPopID=NULL,
        int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL, string sep="\t")
        : Operator<Pop>( "", "", stage, begin, end, step, at, rep, grp, sep)
      {
        DBG_ASSERT( PyObj_Is_IntNumArray(subPopID), ValueError,
          "Passed vector is not a Python/Numeric int array");
        Py_INCREF(subPopID);
        m_subPopID = subPopID;
      }

      /// destructor
      virtual ~PyMigrator()
      {
        if( m_subPopID != NULL)
          Py_DECREF(m_subPopID);
      }

      /// CPPONLY
      PyMigrator(const PyMigrator& rhs):Operator<Pop>(rhs), m_subPopID(rhs.m_subPopID)
      {
        if( m_subPopID != NULL)
          Py_INCREF(m_subPopID);
      }

      virtual Operator<Pop>* clone() const
      {
        return new PyMigrator<Pop>(*this);
      }

      virtual bool apply(Pop& pop)
      {

        DBG_ASSERT( NumArray_Size(m_subPopID) >= static_cast<int>(pop.popSize()) ,
          ValueError, "Given subpopid array has a length of "
          + toStr( NumArray_Size(m_subPopID)) + " which is less than population size "
          + toStr(pop.popSize()));

        long * id = reinterpret_cast<long*>(NumArray_Data(m_subPopID));

        for(size_t i=0, iEnd=pop.popSize(); i<iEnd; ++i)
          pop.individual(i).setInfo( id[i] );

        // do migration.
        // true: rearrange individuals
        pop.setSubPopByIndInfo();

        return true;
      }

      virtual string __repr__()
      {
        return "<simuPOP::python migrator>" ;
      }

    private:
      PyObject* m_subPopID;

  };

  /** \brief
   split subpopulation
   */
  template<class Pop>
    class SplitSubPop: public Operator<Pop>
  {

    public:
      /// split a subpopulation (or whole population as subpop 0)
      /**
      \param which which subpop to split (if there is no subpop structure, 0 is the only subpop)
      \param subPop new subpop size, should add up to the size of subpop to be splitted
      \param proportions proportions of new subpop. (use one of subPop or proportions). Should add up to one.
      \param subPopID optional. new subpop IDs. If given, should have the same length
         as subPop or proportions. SInce subpop with negative id will be removed.
         You can remove part of a subpop by setting a new negative id.
      */
      SplitSubPop( UINT which=0,  vectorlu sizes=vectorlu(), vectorf proportions=vectorf(),
        vectori subPopID=vectori(),
        int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL, string sep="\t")
        : Operator<Pop>( "", "", stage, begin, end, step, at, rep, grp, sep),
        m_which(which), m_subPopSizes(sizes), m_proportions(proportions), m_subPopID(subPopID)
      {
        DBG_FAILIF( sizes.empty() && proportions.empty(), ValueError,
          "Please specify one of subPop and proportions.");
        DBG_FAILIF( !sizes.empty() && !proportions.empty(), ValueError,
          "Please specify only one of subPop and proportions.");
      }

      /// destructor
      virtual ~SplitSubPop()
      {
      }

      virtual Operator<Pop>* clone() const
      {
        return new SplitSubPop<Pop>(*this);
      }

      virtual bool apply(Pop& pop)
      {
        if( !m_subPopSizes.empty())
          pop.splitSubPop(m_which, m_subPopSizes, m_subPopID);
        else
          pop.splitSubPopByProportion(m_which, m_proportions, m_subPopID);
        return true;
      }

      virtual string __repr__()
      {
        return "<simuPOP::split population>" ;
      }

    private:
      /// which subpop to split
      UINT m_which;

      /// new subpopulation size
      vectorlu m_subPopSizes;

      /// new subpopulation proportions.
      vectorf m_proportions;

      /// subpopulation id, optional
      vectori m_subPopID;

  };

  /** \brief
    merge subpopulations
    */
  template<class Pop>
    class MergeSubPops: public Operator<Pop>
  {

    public:
      /// merge subpopulations
      /**
      \param subPops subpops to be merged, default to all subpops.
      */
      MergeSubPops( vectoru subPops=vectoru(), bool removeEmptySubPops=false,
        int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL, string sep="\t")
        : Operator<Pop>( "", "", stage, begin, end, step, at, rep, grp, sep),
        m_subPops(subPops), m_removeEmptySubPops(removeEmptySubPops)
      {
      }

      /// destructor
      virtual ~MergeSubPops()
      {
      }

      virtual Operator<Pop>* clone() const
      {
        return new MergeSubPops<Pop>(*this);
      }

      virtual bool apply(Pop& pop)
      {
        pop.mergeSubPops(m_subPops, m_removeEmptySubPops);
        return true;
      }

      virtual string __repr__()
      {
        return "<simuPOP::merge subpopulations>" ;
      }

    private:
      ///
      vectoru m_subPops;

      ///
      bool m_removeEmptySubPops;
  };

}
#endif
