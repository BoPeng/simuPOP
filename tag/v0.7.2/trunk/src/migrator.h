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

#ifndef _MIGRATOR_H
#define _MIGRATOR_H
/**
\file
\brief head file of class migrator:public Operator
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

	class migrator: public Operator
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
			migrator( const matrix& rate, int mode = MigrByProbability,
				vectoru fromSubPop=vectoru(), vectoru toSubPop=vectoru(),
				int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL)
				: Operator( "", "", stage, begin, end, step, at, rep, grp),
				m_rate(0), m_mode(0), m_from(fromSubPop), m_to(toSubPop)
			{

				DBG_FAILIF( !m_from.empty() && m_from.size() != rate.size(),
					ValueError, "Length of param fromSubPop must match rows of rate matrix.");

				DBG_FAILIF( !m_to.empty() && m_to.size() != rate[0].size(),
					ValueError, "Length of param toSubPop must match columns of rate matrix.");

				setRates(rate, mode);
			};

			/// destructor
			virtual ~migrator(){};

			virtual Operator* clone() const
			{
				return new migrator(*this);
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
			void setRates(const matrix& rate, int mode);

			virtual bool apply(population& pop);

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

	You can use directmigrator to accomplish any migration: that is to say
	you directly specify subpopulation numbers for each individual and
	this operator will do the rest.

	*/

	class pyMigrator: public Operator
	{

		public:
			/// create a directmigrator
			/**
			This operator accept a one-dimensional Numeric Python int array. (created by Numeric.array ).
			The contend of the array will be considered as subpopulation id.

			\param subPopID a 1-d array (list, typle, carray). Must has length greater or equal to
			population size.
			\param stage is default to PreMating, please refer to help(baseOperator.__init__)
			for details about other parameters.
			*/
			pyMigrator( PyObject* subPopID=NULL,
				int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL)
				: Operator( "", "", stage, begin, end, step, at, rep, grp)
			{
				// carray of python list/typle
				DBG_ASSERT( PyObj_Is_IntNumArray(subPopID) ||
					PySequence_Check(subPopID), ValueError,
					"Passed vector is not a sequence (Python list, tuple or carray)");
				Py_INCREF(subPopID);
				m_subPopID = subPopID;
			}

			/// destructor
			virtual ~pyMigrator()
			{
				if( m_subPopID != NULL)
					Py_DECREF(m_subPopID);
			}

			/// CPPONLY
			pyMigrator(const pyMigrator& rhs):Operator(rhs), m_subPopID(rhs.m_subPopID)
			{
				if( m_subPopID != NULL)
					Py_INCREF(m_subPopID);
			}

			virtual Operator* clone() const
			{
				return new pyMigrator(*this);
			}

			virtual bool apply(population& pop);

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

	class splitSubPop: public Operator
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
			splitSubPop( UINT which=0,  vectorlu sizes=vectorlu(), vectorf proportions=vectorf(),
				vectoru subPopID=vectoru(),
				bool randomize=true,
				int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL)
				: Operator( "", "", stage, begin, end, step, at, rep, grp),
				m_which(which), m_subPopSizes(sizes), m_proportions(proportions),
				m_subPopID(subPopID), m_randomize(randomize)
			{
				DBG_FAILIF( sizes.empty() && proportions.empty(), ValueError,
					"Please specify one of subPop and proportions.");
				DBG_FAILIF( !sizes.empty() && !proportions.empty(), ValueError,
					"Please specify only one of subPop and proportions.");
			}

			/// destructor
			virtual ~splitSubPop()
			{
			}

			virtual Operator* clone() const
			{
				return new splitSubPop(*this);
			}

			virtual bool apply(population& pop);

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
			vectoru m_subPopID;

			/// random split
			/// randomize population before split.
			/// this is because some mating schemes generate
			/// individuals non-randomly, for example,
			/// put affected individuals at the beginning.
			bool m_randomize;
	};

	/** \brief
	  merge subpopulations
	  */

	class mergeSubPops: public Operator
	{

		public:
			/// merge subpopulations
			/**
			\param subPops subpops to be merged, default to all subpops.
			*/
			mergeSubPops( vectoru subPops=vectoru(), bool removeEmptySubPops=false,
				int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL)
				: Operator( "", "", stage, begin, end, step, at, rep, grp),
				m_subPops(subPops), m_removeEmptySubPops(removeEmptySubPops)
			{
			}

			/// destructor
			virtual ~mergeSubPops()
			{
			}

			virtual Operator* clone() const
			{
				return new mergeSubPops(*this);
			}

			virtual bool apply(population& pop)
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
