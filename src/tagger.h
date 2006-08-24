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

#ifndef _TAGGER_H
#define _TAGGER_H
/**
\file
\brief head file of class tagger: public Operator
*/
#include "operator.h"

const string TAG_InheritFields[2] = {"paternal_tag", "maternal_tag"};
const string TAG_ParentsFields[2] = {"father_idx", "mother_idx"};

namespace simuPOP
{
	/**
	\brief  tagger is a during mating operator that
	tag individual with various information. Potential usages are
	1. record parenting information to track pedigree.
	2. tag a individual/allele and monitor its spread in the population
	   etc.
	3...

	@author Bo Peng
	*/


	class tagger: public Operator
	{

		public:
			/// constructor. default to be always active but no output.
			tagger( int begin=0, int end=-1, int step=1, vectorl at=vectorl(), 
				int rep=REP_ALL, int grp=GRP_ALL, 
				// this is not nice, but is the only way I know how to initialize this array.
				const vectorstr& infoFields=vectorstr()):
			Operator("", "", DuringMating, begin, end, step, at, rep, grp, infoFields)
			{
			};

			/// destructor
			virtual ~tagger(){};

			virtual Operator* clone() const
			{
				return new tagger(*this);
			}
	};

	/// inherite tag from parents.
	/// If both parents have tags, use fathers.

	class inheritTagger: public tagger
	{
		public:
#define TAG_Paternal   0
#define TAG_Maternal   1
#define TAG_Both       2

		public:
			/// constructor. default to be always active.
			inheritTagger(int mode=TAG_Paternal, int begin=0, int end=-1, int step=1, 
				vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr(TAG_InheritFields, TAG_InheritFields+2)):
			tagger( begin, end, step, at, rep, grp, infoFields), m_mode(mode)
			{
			};

			virtual ~inheritTagger()
			{
			}

			virtual string __repr__()
			{
				return "<simuPOP::inherittagger>" ;
			}

			virtual bool applyDuringMating(population& pop, population::IndIterator offspring,
				individual* dad=NULL, individual* mom=NULL)
			{
				UINT id1=0, id2=0;
				if (m_mode == TAG_Paternal)
					id1 = pop.infoIdx(infoField(0));
				else if (m_mode == TAG_Maternal)
					id1 = pop.infoIdx(infoField(1));
				else
				{
					id1 = pop.infoIdx(infoField(0));
					id2 = pop.infoIdx(infoField(1));
				}
				
				if( m_mode == TAG_Paternal)
				{
					if (dad == NULL)
						offspring->setInfo(0, id1);
					else
						offspring->setInfo(dad->info(id1), id1);
				}
				else if( m_mode == TAG_Maternal)
				{
					if (mom == NULL)
						offspring->setInfo(0, id1);
					else
						offspring->setInfo(mom->info(id1), id1);
				}
				else
				{
					if( dad == NULL )
						offspring->setInfo(0, id1);
					else
						offspring->setInfo(dad->info(id1), id1);
					if( mom == NULL)
						offspring->setInfo(0, id2);
					else
						offspring->setInfo(mom->info(id2), id2);
				}
				return true;
			}

		
			virtual Operator* clone() const
			{
				return new inheritTagger(*this);
			}

		private:
			/// mode can be
			/// TAG_Paternal: get dad's info
			/// TAG_Maternal: get mon's info
			/// TAG_BOTH:     get parents' first field
			int m_mode;
	};

	/// inherite tag from parents.
	/// If both parents have tags, use fathers.
	///

	class parentsTagger: public tagger
	{
		public:
			/// constructor. default to be always active.
			/// string can be any string (m_Delimiter will be ignored for this class.)
			///  %r will be replicate number %g will be generation number.
			parentsTagger( int begin=0, int end=-1, int step=1, vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL, 
				const vectorstr& infoFields=vectorstr(TAG_ParentsFields, TAG_ParentsFields+2)):
			tagger( begin, end, step, at, rep, grp, infoFields)
			{
			};

			virtual ~parentsTagger()
			{
			}

			virtual Operator* clone() const
			{
				return new parentsTagger(*this);
			}

			virtual string __repr__()
			{
				return "<simuPOP::parentstagger>" ;
			}

			virtual bool applyDuringMating(population& pop, population::IndIterator offspring,
				individual* dad=NULL, individual* mom=NULL)
			{
				UINT id1 = pop.infoIdx(infoField(0));
				UINT id2 = pop.infoIdx(infoField(1));
				
				if(dad == NULL)
					offspring->setInfo(0, id1);
				else
					offspring->setInfo(dad - &*pop.indBegin(), id1);
				
				if( mom == NULL)
					offspring->setInfo(0, id2);
				else
					offspring->setInfo(mom - &*pop.indBegin(), id2);

				return true;
			}
	};
}
#endif
