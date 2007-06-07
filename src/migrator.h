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

#include "simuPOP_cfg.h"
#include <math.h>
// using std::sqrt;

namespace simuPOP
{

    /// migrate individuals from a (sub)population to another (sub)population
	class migrator: public Operator
	{

		public:

#define MigrByProbability  1
#define MigrByProportion   2
#define MigrByCounts       3

		public:
			/// create a migrator
			/**
            Migrator is the only way to mix genotypes of several subpopulations 
            because mating is strictly within subpopulations in simuPOP. Migrators
            are quite flexible in simuPOP in the sense that
            \li Migration can happen from and to a subset of subpopulations.
            \li Migration can be done by probability, proportion or by counts. In
                the case of probability, if the migration rate from subpopulation
                \c a to \c b is \c r, then everyone in subpopulation \c a will have this
                probability to migrate to \c b. In the case of proportion, exactly
                <tt>r*size_of_subPop_a</tt> individuals (chosen by random) will migrate
                to subpopulation \c b. In the last case, a given number of individuals will
                migrate.
            \li New subpopulation can be generated through migration. You simply
                need to migrate to a new subpopulation number.
                
			\param rate migration rate, can be a proportion or counted number. Determined by
                parameter \c mode. \c rate should be an m by n matrix. If a number is given,
                the migration rate will be <tt>r*ones(m,n)</tt>???.
			\param mode one of \c MigrByProbability (default), \c MigrByProportion or \c MigrByCounts
			\param fromSubPop an array of 'from' subpopulations. Default to all. If a single \c subpop
                is specified, <tt>[]</tt> can be ignored. I.e., <tt>[a]<tt> is equvalent to \c a.
			\param toSubPop an array of 'to' subpopulations. Default to all subpopulations. If a single
                \c subpop is specified, <tt>[]</tt> can be ignored.
			\param stage default to \c PreMating

            \note
                \li The overall population size will not be changed. (Mating schemes can
                    do that). If you would like to keep the subpopulation size after migration, you
                    can use the \c newSubPopSize or \c newSubPopSizeExpr parameter of a mating scheme.
                \li \c rate is a matrix with dimensions determined by \c fromSubPop and \c toSubPop.
                    By default, \c rate is a matrix with element \c r(i,j), where \c r(i, j) is the
                    migration rate, probability or count from subpopulation \c i to \c j. If \c fromSubPop
                    and/or \c toSubPop are given, migration will only happen between these subpopulations.
                    An extreme case is 'point migration', <tt>rate=[[r]], fromSubPop=a, toSubPop=b</tt>
                    which migrate from subpopulation \c a to \c b with given rate \c r.???
			*/
			migrator( const matrix& rate, int mode = MigrByProbability,
				vectoru fromSubPop=vectoru(), vectoru toSubPop=vectoru(),
				int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				: Operator( "", "", stage, begin, end, step, at, rep, grp, infoFields),
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

            /// deep copy of a migrator
			virtual Operator* clone() const
			{
				return new migrator(*this);
			}

			/// return migration rate
			matrix rate()
			{
				return m_rate;
			}

			/// set migration rate
			/**
			Format should be <tt>0-0 0-1 0-2, 1-0 1-1 1-2, 2-0, 2-1, 2-2</tt>.
			For mode \c MigrByProbability or \c MigrByProportion, <tt>0-0,1-1,2-2</tt> will be set
			automatically regardless of input.
			*/
			void setRates(const matrix& rate, int mode);

            /// apply the \c migrator
			virtual bool apply(population& pop);

            /// used by Python print function to print out the general information of the \c migrator
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

    /// migrate using given information vector
	/**
	You can use directmigrator??? to accomplish any migration: that is to say
	you directly specify subpopulation numbers for each individual and
	this operator will do the rest.
	*/
	class pyMigrator: public Operator
	{

		public:
			/// create a directmigrator
			/**
            For even more complicated migration schemes, you can DIY it using a
            \c pyMigrator. This operator is not strictly hybrid since it does not
            call a python function. However, it takes a carray of subPopulation IDs
            for each individual. \c pyMigrator then complete migration according to
            its content.
        
			This operator accepts a one-dimensional Numeric Python \c int array
            (created by \c Numeric.array ). ???

			\param subPopID a 1-d array (list, typle, carray). Length must be
                greater than or equal to the population size.
			\param stage default to \c PreMating
			*/
			pyMigrator( PyObject* subPopID=NULL,
				int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				: Operator( "", "", stage, begin, end, step, at, rep, grp, infoFields)
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

            /// deep copy of a \c pyMigrator
			virtual Operator* clone() const
			{
				return new pyMigrator(*this);
			}

            /// apply a \c pyMigrator
			virtual bool apply(population& pop);

            /// used by Python print function to print out the general information of the \c pyMigrator
			virtual string __repr__()
			{
				return "<simuPOP::python migrator>" ;
			}

		private:
			PyObject* m_subPopID;

	};

	/// split a subpopulation
	class splitSubPop: public Operator
	{

		public:
			/// split a subpopulation or the whole population as subpopulation \c 0
			/**
			\param which which subpopulation to split. If there is no subpopulation structure,
                use 0 as the first (and only) subpopulation.
			\param sizes new subpopulation sizes. The sizes should be added up to the original
                subpopulation (subpopulation \c which) size.
			\param proportions proportions of new subpopulations. Should be added up to \c 1.
                Optionally, you can use one of \c subPopID or \c proportions to split. ???
			\param subPopID new subpopulation IDs. Otherwise, the operator will automatically
                set new subpopulation IDs to new subpopulations. If given, should have the same length
			   as subPop or proportions.??? Since subpop with negative id will be removed.
			   You can remove part of a subpop by setting a new negative id.???
			*/
			splitSubPop( UINT which=0,  vectorlu sizes=vectorlu(), vectorf proportions=vectorf(),
				vectoru subPopID=vectoru(),
				bool randomize=true,
				int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				: Operator( "", "", stage, begin, end, step, at, rep, grp, infoFields),
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

            /// deep copy of a \c splitSubPop operator
			virtual Operator* clone() const
			{
				return new splitSubPop(*this);
			}

            /// apply a \c splitSubPop operator
			virtual bool apply(population& pop);

            /// used by Python print function to print out the general information of the \c splitSubPop operator
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

	///  merge subpopulations
	class mergeSubPops: public Operator
	{

		public:
			/// merge subpopulations
			/**
            This operator merges subPopulations \c subPops (the only parameter) to a
            single subpopulation. If \c subPops is ignored, all subpopulations will be merged.
			\param subPops subpopulatinos to be merged. Default to all.
			*/
			mergeSubPops( vectoru subPops=vectoru(), bool removeEmptySubPops=false,
				int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				: Operator( "", "", stage, begin, end, step, at, rep, grp, infoFields),
				m_subPops(subPops), m_removeEmptySubPops(removeEmptySubPops)
			{
			}

			/// destructor
			virtual ~mergeSubPops()
			{
			}

            /// deep copy of a \c mergeSubPops operator
			virtual Operator* clone() const
			{
				return new mergeSubPops(*this);
			}

            /// apply a \c mergeSubPops operator
			virtual bool apply(population& pop)
			{
				pop.mergeSubPops(m_subPops, m_removeEmptySubPops);
				return true;
			}

            /// used by Python print function to print out the general information of the \c mergeSubPops operator
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
