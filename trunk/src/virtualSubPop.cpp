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


infoSplitter::infoSplitter(string info, vectorinfo const & values,
                           vectorf const & cutoff)
	: vspSplitter(),
	m_info(info), m_values(values), m_cutoff(cutoff)
{
	DBG_FAILIF(m_values.empty() && m_cutoff.empty(),
	    ValueError, "Please specify either a list of values, or a set of cutoff values");
	DBG_FAILIF(!m_values.empty() && !m_cutoff.empty(),
	    ValueError, "Please specify only a list of values, or a set of cutoff values");
	// cutoff value has to be ordered
	if (m_values.empty())
		for (size_t i = 1; i < m_cutoff.size(); ++i) {
			DBG_ASSERT(m_cutoff[i - 1] < m_cutoff[i],
			    ValueError,
			    "Cutoff values have to be in increasing order");
		}
}


ULONG infoSplitter::size(const population & pop, virtualSubPopID subPop) const
{
	UINT idx = pop.infoIdx(m_info);

	ULONG count = 0;

	ConstRawIndIterator it = pop.rawIndBegin(subPop.id());
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop.id());

	if (m_values.empty()) {
		DBG_FAILIF(subPop.vid() > m_cutoff.size(), IndexError,
		    "Virtual Subpoplation index out of range of 0 ~ "
		    + toStr(m_cutoff.size()));

		// using cutoff, below
		if (subPop.vid() == 0) {
			for (; it != it_end; ++it)
				if (it->info(idx) < m_cutoff[0])
					count++;
			return count;
		} else if (subPop.vid() == m_cutoff.size()) {
			double v = m_cutoff.back();
			for (; it != it_end; ++it)
				if (it->info(idx) >= v)
					count++;
			return count;
		} else {         // in between
			double v1 = m_cutoff[subPop.vid() - 1];
			double v2 = m_cutoff[subPop.vid()];
			for (; it != it_end; ++it) {
				double v = it->info(idx);
				if (v >= v1 && v < v2)
					count++;
			}
			return count;
		}
	} else {
		DBG_FAILIF(subPop.vid() >= m_values.size(), IndexError,
		    "Virtual Subpoplation index out of range of 0 ~ "
		    + toStr(m_values.size() - 1));
		double v = m_values[subPop.vid()];
		for (; it != it_end; ++it)
			if (fcmp_eq(it->info(idx), v))
				count++;
		return count;
	}
	// should never reach here.
	return 0;
}


UINT infoSplitter::numVirtualSubPop()
{
	if (m_values.empty())
		return m_cutoff.size() + 1;
	else
		return m_values.size();
}


void infoSplitter::activate(population & pop, virtualSubPopID subPop)
{
	if (m_activated && m_vsp == subPop.vid())
		return;

	UINT idx = pop.infoIdx(m_info);

	RawIndIterator it = pop.rawIndBegin(subPop.id());
	RawIndIterator it_end = pop.rawIndEnd(subPop.id());
	if (m_values.empty()) {
		DBG_FAILIF(subPop.id() > m_cutoff.size(), IndexError,
		    "Virtual Subpoplation index out of range of 0 ~ "
		    + toStr(m_cutoff.size()));

		// using cutoff, below
		if (subPop.vid() == 0) {
			for (; it != it_end; ++it)
				if (it->info(idx) >= m_cutoff[0])
					it->setVisible(false);
				else if (m_activated)
					it->setVisible(true);
		} else if (subPop.vid() == m_cutoff.size()) {
			double v = m_cutoff.back();
			for (; it != it_end; ++it)
				if (it->info(idx) < v)
					it->setVisible(false);
				else if (m_activated)
					it->setVisible(true);
		} else {         // in between
			double v1 = m_cutoff[subPop.vid() - 1];
			double v2 = m_cutoff[subPop.vid()];
			double v = it->info(idx);
			for (; it != it_end; ++it)
				if (v < v1 || v >= v2)
					it->setVisible(false);
				else if (m_activated)
					it->setVisible(true);
		}
	} else {
		DBG_FAILIF(subPop.id() >= m_values.size(), IndexError,
		    "Virtual Subpoplation index out of range of 0 ~ "
		    + toStr(m_values.size() - 1));
		for (; it != it_end; ++it) {
			double v = m_values[subPop.vid()];
			if (fcmp_ne(it->info(idx), v))
				it->setVisible(false);
			else if (m_activated)
				it->setVisible(true);
		}
	}
	m_activated = true;
	m_vsp = subPop.vid();
}


void infoSplitter::reset(population & pop, SubPopID subPop)
{
	resetSubPop(pop, subPop);
}


string infoSplitter::name(SubPopID sp)
{
	if (m_values.empty()) {
		DBG_FAILIF(sp > m_cutoff.size(), IndexError,
		    "Virtual Subpoplation index out of range of 0 ~ "
		    + toStr(m_cutoff.size()));
		if (sp == 0)
			return "Info < " + toStr(m_cutoff[0]);
		else if (sp == m_cutoff.size())
			return "Info >= " + toStr(m_cutoff[sp - 1]);
		else
			return toStr(m_cutoff[sp - 1]) + " <= Info < "
			       + toStr(m_cutoff[sp]);
	} else {
		DBG_FAILIF(sp >= m_values.size(), IndexError,
		    "Virtual Subpoplation index out of range of 0 ~ "
		    + toStr(m_values.size() - 1));
		return "Info " + toStr(m_values[sp]);
	}


}


} // namespce simuPOP

