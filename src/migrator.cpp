/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu                                                        *
 *                                                                         *
 *   $LastChangedDate: 2006-02-21 15:27:25 -0600 (Tue, 21 Feb 2006) $
 *   $Rev: 191 $
 *
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

#include "simupop_cfg.h"
#include "utility.h"

#include "migrator.h"

namespace simuPOP
{
  void migrator::setRates(const matrix& rate, int mode)
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
      for(UINT i=0; i< szFrom; i++)               // from
      {
        // look for from=to cell.
        UINT spFrom = m_from[i];
        double sum = accumulate(m_rate[i].begin(), m_rate[i].end(), 0.0);
        //
        vectoru::iterator spTo = find(m_to.begin(), m_to.end(), spFrom);
        if( spTo == m_to.end() )                  // if no to, only check if sum <= 1
        {
          if( fcmp_gt(sum, 1.0) )
            throw ValueError("Sum of migrate rate from one subPop should <= 1");
          // adding i->i item
          m_rate[i].push_back(1.0 - sum);
        }
        else                                      // if not, set r[i][i]
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
  }

  bool migrator::apply(population& pop)
  {
    // set info of individual
    pop.setIndInfoWithSubPopID();

    population::IndIterator ind, indEd;

    vectorlu toIndices(0);

    Weightedsampler ws(rng());

    for(UINT from=0, fromEnd=m_from.size(); from < fromEnd; ++from)
    {
      UINT spFrom = m_from[from];
      // rateSize might be toSize + 1, the last one is from->from
      UINT toSize = m_to.size(), toIndex;

      // m_from out of range.... ignore.
      if( spFrom >= pop.numSubPop() )
        continue;

      if(m_mode == MigrByProbability )            // migrate by probability
      {
        ws.set(m_rate[from]);

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
      else                                        // by count
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
    }                                             /// for all subPop.

    // do migration.
    // true: rearrange individuals
    pop.setSubPopByIndInfo();

    return true;
  }

  bool pyMigrator::apply(population& pop)
  {
    if(PyObj_Is_IntNumArray(m_subPopID) )
    {
      DBG_ASSERT( NumArray_Size(m_subPopID) >= static_cast<int>(pop.popSize()) ,
        ValueError, "Given subpopid array has a length of "
        + toStr( NumArray_Size(m_subPopID)) + " which is less than population size "
        + toStr(pop.popSize()));

      long * id = reinterpret_cast<long*>(NumArray_Data(m_subPopID));

      for(size_t i=0, iEnd=pop.popSize(); i<iEnd; ++i)
        pop.ind(i).setInfo( id[i] );
    }
    else
    {
      DBG_ASSERT( PySequence_Size(m_subPopID) >=  static_cast<int>(pop.popSize()) ,
        ValueError, "Given subpopid array has a length of "
        + toStr( PySequence_Size(m_subPopID)) + " which is less than population size "
        + toStr(pop.popSize()));

      int id;
      for(size_t i=0, iEnd=pop.popSize(); i<iEnd; ++i)
      {
        PyObj_As_Int(PySequence_GetItem(m_subPopID, i), id);
        pop.ind(i).setInfo(id);
      }
    }
    // do migration.
    // true: rearrange individuals
    pop.setSubPopByIndInfo();
    return true;
  }

}
