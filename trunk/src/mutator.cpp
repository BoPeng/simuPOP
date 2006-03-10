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
#include "mutator.h"

namespace simuPOP
{

  void mutator::initialize(population& pop)
  {
#ifndef BINARYALLELE
    if( m_maxAllele == 0 )
      m_maxAllele = pop.maxAllele();
    else if ( m_maxAllele > 0 && m_maxAllele > pop.maxAllele() )
      throw ValueError("maxAllele exceeds population max allele.");
#endif

    DBG_DO(DBG_MUTATOR, cout << "initialize mutator" << endl);

    // deal with applicable loci
    if(m_atLoci.empty() )
    {
      // all loci
      m_atLoci.resize(pop.totNumLoci() );
      for(UINT i=0, iEnd = pop.totNumLoci(); i < iEnd;  ++i)
        m_atLoci[i] = i;
    }

    /// all use the same rate
    if( m_rate.size() < m_atLoci.size() )
    {
      m_rate.resize( m_atLoci.size());
      fill(m_rate.begin()+1, m_rate.end(), m_rate[0]);
    }

    m_bt.setParameter(m_rate, pop.ploidy() * pop.popSize());

#ifndef OPTIMIZED
    for(size_t i=0; i<m_rate.size(); ++i)
      if( fcmp_lt( m_rate[i], 0.) || fcmp_gt( m_rate[i], 1.) )
        throw ValueError("Migration rate should be between [0,1], given " + toStr(m_rate[i]));
#endif
    m_mutCount.resize(pop.totNumLoci(), 0);
    m_initialized = true;
  }

  bool mutator::apply(population& pop)
  {
    if( !m_initialized || m_bt.size() != pop.ploidy() * pop.popSize())
    {
      initialize(pop);
      DBG_DO(DBG_MUTATOR, cout << "Mutate at loci " << m_atLoci <<
        " at rate " << m_rate << endl);
    }

    DBG_DO(DBG_MUTATOR, cout <<"Mutate replicate " << pop.rep() << endl);

    m_bt.doTrial();

    // mutate each mutable locus
    for( size_t i=0, iEnd=m_atLoci.size(); i < iEnd; ++i)
    {
      const BitSet& succ = m_bt.succ(i);

      int locus = m_atLoci[i];

      BitSet::size_type pos = succ.find_first();
      if( pos != BitSet::npos)
      {
        do
        {
#ifndef OPTIMIZED
          AlleleRef ptr = *(pop.alleleBegin( locus ) + pos).ptr();
          DBG_DO(DBG_MUTATOR, cout << "Mutate locus " << locus
            << " of individual " << (pos/pop.ploidy()) << " from " << int(ptr) );
          mutate(ptr);
#else
          mutate( *(pop.alleleBegin( locus ) + pos).ptr() );
#endif
          m_mutCount[ locus ]++;
        }while( (pos = succ.find_next(pos)) != BitSet::npos );
      }                                           // succ.any
    }                                             // each applicable loci

    return true;
  }

  /// mutate to a state other than current state with equal probability
  void kamMutator::mutate(AlleleRef allele)
  {
#ifdef BINARYALLELE
    allele = !allele;
#else
    Allele new_allele = rng().randInt(this->maxAllele());
    if(new_allele >= allele)
      allele = new_allele+1;
    else
      allele = new_allele;
#endif
  }

  void gsmMutator::mutate(AlleleRef allele)
  {
    int step;

    if( m_func == NULL)                           // use a geometric distribution.
      step = rng().randGeometric(m_p);
    else
    {
      PyObject* arglist = Py_BuildValue("()");
      PyObject* result = PyEval_CallObject(m_func, arglist);
      Py_DECREF(arglist);
      if( result == NULL)
      {
        PyErr_Print();
        throw ValueError("Function call failed.");
      }

      PyObj_As_Int(result, step);
      Py_DECREF(result);
    }

    DBG_DO(DBG_MUTATOR, cout << "step is " << step << endl);

    if( rng().randUniform01() < m_incProb)
    {
      if( static_cast<UINT>(allele + step) < this->maxAllele() )
        AlleleAdd(allele, step);
      else
        allele = this->maxAllele();
    }
    else
    {
      if( allele > step)
        AlleleMinus(allele, step);
      else
        allele = 0;
    }
  }

}
