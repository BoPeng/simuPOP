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

#include "operator.h"

namespace simuPOP
{

	bool Operator::isActive(UINT rep, UINT numRep, long gen, long end, int grp, bool repOnly )
	{
		// rep does not match
		if( ( m_rep >= 0 && static_cast<UINT>(m_rep) != rep ) ||
			( m_rep == REP_LAST && rep != numRep -1 ) )
			return false;

		// group does not match
		if( m_grp >= 0 && m_grp != grp )
			return false;

		// only check for rep and grp values.
		if( repOnly )
			return true;

		// if gen < 0, we are testing only group info
		if( gen < 0 )
			return true;

		// if all active? (begin=0, end=-1)
		if( ISSETFLAG(m_flags, m_flagAtAllGen))
			return true;

		DBG_FAILIF( end > 0 && gen > end, IndexError,
			"Current generation can not be bigger than ending generation.");

		if( ISSETFLAG( m_flags, m_flagOnlyAtBegin) )
		{
			if(gen == 0) return true;
			else return false;
		}

		if( ISSETFLAG(m_flags, m_flagOnlyAtEnd) && end > 0 )
		{
			if( gen == end ) return true;
			else return false;
		}

		// at Gen has higher priority.
		if( ! m_atGen.empty() )
		{
			// chech atGen.
			for(size_t i=0, iEnd = m_atGen.size(); i < iEnd;  ++i)
			{
				int atGen = m_atGen[i];

				if( atGen > 0 )
				{
					if( gen == atGen )  return true;
					else continue;
				}								  // atGen <=0
				else
				{
					if( end < 0 )				  // can not determine.
						continue;
					// now end >= 0 atGen <=0
					if( end + atGen + 1 == gen )
						return true;
					else
						continue;
				}
			}
			// Do not check other parameters
			return false;
		}

		/// finally, check start, end, every
		if( end < 0 )							  // if we do not know ending generation.
		{
			// can not determine starting gen.
			if( m_beginGen < 0 || m_beginGen > gen )
				return false;

			// can only check positive begin + every
			if( ( ((gen - m_beginGen) % m_stepGen) == 0) && (m_endGen < 0 || m_endGen >= gen) )
				return true;
			else
				return false;
		}										  // know ending generation
		else
		{
			int realStartGen = m_beginGen >= 0 ? m_beginGen : m_beginGen + end + 1;
			int realEndGen = m_endGen >= 0 ? m_endGen : m_endGen + end + 1;

			if(realStartGen > realEndGen)
				return false;

			return(( gen >= realStartGen && gen <= realEndGen && (gen - realStartGen)%m_stepGen == 0 ) );
		}
		return false;
	}

	void Operator::setFlags()
	{
		RESETFLAG(m_flags, m_flagAtAllGen);
		RESETFLAG(m_flags, m_flagOnlyAtBegin);
		RESETFLAG(m_flags, m_flagOnlyAtEnd);

		/// atGen has higher priority: if it is not empty, use it.
		if( m_atGen.empty() )
		{
			if( m_beginGen == 0 && m_endGen == 0 )
				SETFLAG(m_flags, m_flagOnlyAtBegin);

			if( m_endGen == -1 && m_beginGen == -1 )
				SETFLAG(m_flags, m_flagOnlyAtEnd);

			if( m_beginGen == 0 && m_endGen == -1 && m_stepGen == 1 )
				SETFLAG(m_flags, m_flagAtAllGen);

		}
		else if( m_atGen.size() == 1)
		{
			if( m_stepGen < 1 )
				throw IndexError("active generation interval should be at least 1.");

			if( m_atGen[0] == 0 )  SETFLAG(m_flags, m_flagOnlyAtBegin);
			if( m_atGen[0] == -1 ) SETFLAG(m_flags, m_flagOnlyAtEnd);
		}
	}

}
