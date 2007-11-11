/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu                                                        *
*                                                                         *
*   $LastChangedDate: 2006-02-21 15:27:25 -0600 (Tue, 21 Feb 2006)        *
*   $Rev: 191$                                                            *
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

#include "simuPOP_cfg.h"
#include "virtualSubPop.h"
#include "population.h"

namespace simuPOP {

void vspSplitter::resetSubPop(population & pop, SubPopID subPop)
{
	if (!m_activated)
		return;

	RawIndIterator it = pop.rawIndBegin(subPop);
	RawIndIterator it_end = pop.rawIndEnd(subPop);
	for (; it != it_end; ++it)
		it->setVisible(true);
}


ULONG nullSplitter::size(const population & pop, virtualSubPopID subPop) const
{
	DBG_ASSERT(subPop.isVirtual(), ValueError, "Given virtual subpop id is not virtual");
	DBG_ASSERT(subPop.vid() == 0, IndexError, "There is only one virtual subpopulation for this splitter");
	return pop.subPopSize(subPop.id());
}


ULONG duplicateSplitter::size(const population & pop, virtualSubPopID subPop) const
{
	DBG_ASSERT(subPop.isVirtual(), ValueError, "Given virtual subpop id is not virtual");
	DBG_ASSERT(subPop.vid() < m_num, IndexError, "Virtual subpopulation id out of range");
	return pop.subPopSize(subPop.id());
}


ULONG sexSplitter::size(const population & pop, virtualSubPopID subPop) const
{
	ConstRawIndIterator it = pop.rawIndBegin(subPop.id());
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop.id());
	ULONG count = 0;
	Sex s = subPop.vid() == 0 ? Male : Female;

	for (; it != it_end; ++it)
		if (it->sex() == s)
			++count;
	return count;
}


void sexSplitter::activate(population & pop, virtualSubPopID subPop)
{
	Sex s = subPop.vid() == 0 ? Male : Female;

	if (m_activated && m_sex == s)
		return;

	RawIndIterator it = pop.rawIndBegin(subPop.id());
	RawIndIterator it_end = pop.rawIndEnd(subPop.id());
	for (; it != it_end; ++it) {
		if (it->sex() != s)
			// only turn off opposite sex
			it->setVisible(false);
		else if (m_activated)
			// switch? set all
			it->setVisible(true);
	}
	m_sex = s;
	m_activated = true;
}


void sexSplitter::reset(population & pop, SubPopID subPop)
{
	resetSubPop(pop, subPop);
}


ULONG affectionSplitter::size(const population & pop, virtualSubPopID subPop) const
{
	ConstRawIndIterator it = pop.rawIndBegin(subPop.id());
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop.id());
	ULONG count = 0;
	// 0 is unaffected
	bool aff = subPop.vid() == 0 ? false : true;

	for (; it != it_end; ++it)
		if (it->affected() == aff)
			++count;
	return count;
}


void affectionSplitter::activate(population & pop, virtualSubPopID subPop)
{
	bool aff = subPop.vid() == 0 ? false : true;

	if (m_activated && m_affection == aff)
		return;

	RawIndIterator it = pop.rawIndBegin(subPop.id());
	RawIndIterator it_end = pop.rawIndEnd(subPop.id());
	for (; it != it_end; ++it) {
		if (it->affected() != aff)
			// only turn off opposite affection status
			it->setVisible(false);
		else if (m_activated)
			// switch? set all
			it->setVisible(true);
	}
	m_affection = aff;
	m_activated = true;
}


void affectionSplitter::reset(population & pop, SubPopID subPop)
{
	resetSubPop(pop, subPop);
}


}
