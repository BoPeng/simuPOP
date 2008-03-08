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
	m_activated = false;
}


ULONG vspSplitter::countVisibleInds(const population & pop, SubPopID subPop) const
{
	if (!activated())
		return pop.subPopSize(subPop);
	ULONG count = 0;
	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);
	for (; it != it_end; ++it)
		if (it->visible())
			++count;
	return count;
}


combinedSplitter::combinedSplitter(const vectorvsp & splitters)
	: vspSplitter(), m_numVSP(0), m_splitter(), m_vsp(), m_curSplitter(0)
{
	for (size_t i = 0; i < splitters.size(); ++i) {
		m_splitters.push_back(splitters[i]->clone());
		for (size_t j = 0; j < splitters[i]->numVirtualSubPop(); ++j) {
			m_splitter.push_back(i);
			m_vsp.push_back(j);
			m_numVSP ++;
		}
	}
}
	
combinedSplitter::~combinedSplitter()
{
	for (size_t i = 0; i < m_splitters.size(); ++i)
		delete m_splitters[i];
}


vspSplitter * combinedSplitter::clone() const
{
	return new combinedSplitter(m_splitters);
}	


ULONG combinedSplitter::size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	return m_splitters[m_splitter[virtualSubPop]]->size(pop, subPop,
		m_vsp[virtualSubPop]);
}


void combinedSplitter::activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		activateType type)
{
	m_curSplitter = m_splitter[virtualSubPop];
	m_splitters[m_curSplitter]->activate(pop, subPop,
		m_vsp[virtualSubPop], type);
}


void combinedSplitter::deactivate(population & pop, SubPopID sp)
{
	m_splitters[m_curSplitter]->deactivate(pop, sp);
}


string combinedSplitter::name(SubPopID sp)
{
	return m_splitters[m_splitter[sp]]->name(m_vsp[sp]);

}


ULONG sexSplitter::size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	if (virtualSubPop == InvalidSubPopID)
		return countVisibleInds(pop, subPop);
	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);
	ULONG count = 0;
	Sex s = virtualSubPop == 0 ? Male : Female;

	for (; it != it_end; ++it)
		if (it->sex() == s)
			++count;
	return count;
}


void sexSplitter::activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
                           activateType type)
{
	Sex s = virtualSubPop == 0 ? Male : Female;

	RawIndIterator it = pop.rawIndBegin(subPop);
	RawIndIterator it_end = pop.rawIndEnd(subPop);

	for (; it != it_end; ++it)
		if (type == Visible)
			it->setVisible(it->sex() == s);
		else
			it->setIteratable(it->sex() == s);
	if (type == Visible)
		m_activated = true;
}


void sexSplitter::deactivate(population & pop, SubPopID subPop)
{
	resetSubPop(pop, subPop);
}


ULONG affectionSplitter::size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	if (virtualSubPop == InvalidSubPopID)
		return countVisibleInds(pop, subPop);
	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);
	ULONG count = 0;
	// 0 is unaffected
	bool aff = virtualSubPop == 0 ? false : true;

	for (; it != it_end; ++it)
		if (it->affected() == aff)
			++count;
	return count;
}


void affectionSplitter::activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
                                 activateType type)
{
	bool aff = virtualSubPop == 0 ? false : true;

	RawIndIterator it = pop.rawIndBegin(subPop);
	RawIndIterator it_end = pop.rawIndEnd(subPop);

	for (; it != it_end; ++it)
		if (type == Visible)
			it->setVisible(it->affected() == aff);
		else
			it->setIteratable(it->affected() == aff);
	if (type == Visible)
		m_activated = true;
}


void affectionSplitter::deactivate(population & pop, SubPopID subPop)
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


ULONG infoSplitter::size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	if (virtualSubPop == InvalidSubPopID)
		return countVisibleInds(pop, subPop);
	UINT idx = pop.infoIdx(m_info);

	ULONG count = 0;

	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);

	if (m_values.empty()) {
		DBG_FAILIF(static_cast<UINT>(virtualSubPop) > m_cutoff.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_cutoff.size()));

		// using cutoff, below
		if (virtualSubPop == 0) {
			for (; it != it_end; ++it)
				if (it->info(idx) < m_cutoff[0])
					count++;
			return count;
		} else if (static_cast<UINT>(virtualSubPop) == m_cutoff.size()) {
			double v = m_cutoff.back();
			for (; it != it_end; ++it)
				if (it->info(idx) >= v)
					count++;
			return count;
		} else {         // in between
			double v1 = m_cutoff[virtualSubPop - 1];
			double v2 = m_cutoff[virtualSubPop];
			for (; it != it_end; ++it) {
				double v = it->info(idx);
				if (v >= v1 && v < v2)
					count++;
			}
			return count;
		}
	} else {
		DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_values.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_values.size() - 1));
		double v = m_values[virtualSubPop];
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


void infoSplitter::activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
                            activateType type)
{
	UINT idx = pop.infoIdx(m_info);

	RawIndIterator it = pop.rawIndBegin(subPop);
	RawIndIterator it_end = pop.rawIndEnd(subPop);

	if (m_values.empty()) {
		DBG_FAILIF(static_cast<UINT>(subPop) > m_cutoff.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_cutoff.size()));

		// using cutoff, below
		if (virtualSubPop == 0) {
			for (; it != it_end; ++it)
				if (type == Visible)
					it->setVisible(it->info(idx) < m_cutoff[0]);
				else
					it->setIteratable(it->info(idx) < m_cutoff[0]);
		} else if (static_cast<UINT>(virtualSubPop) == m_cutoff.size()) {
			double v = m_cutoff.back();
			for (; it != it_end; ++it)
				if (type == Visible)
					it->setVisible(it->info(idx) >= v);
				else
					it->setIteratable(it->info(idx) >= v);
		} else {         // in between
			double v1 = m_cutoff[virtualSubPop - 1];
			double v2 = m_cutoff[virtualSubPop];
			for (; it != it_end; ++it) {
				double v = it->info(idx);
				if (type == Visible)
					it->setVisible(v >= v1 && v < v2);
				else
					it->setIteratable(v >= v1 && v < v2);
			}
		}
	} else {
		DBG_FAILIF(static_cast<UINT>(subPop) >= m_values.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_values.size() - 1));
		double v = m_values[virtualSubPop];
		for (; it != it_end; ++it)
			if (type == Visible)
				it->setVisible(fcmp_eq(it->info(idx), v));
			else
				it->setIteratable(fcmp_eq(it->info(idx), v));
	}
	if (type == Visible)
		m_activated = true;
}


void infoSplitter::deactivate(population & pop, SubPopID subPop)
{
	resetSubPop(pop, subPop);
}


string infoSplitter::name(SubPopID sp)
{
	if (m_values.empty()) {
		DBG_FAILIF(static_cast<UINT>(sp) > m_cutoff.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_cutoff.size()));
		if (sp == 0)
			return m_info + " < " + toStr(m_cutoff[0]);
		else if (static_cast<UINT>(sp) == m_cutoff.size())
			return m_info + " >= " + toStr(m_cutoff[sp - 1]);
		else
			return toStr(m_cutoff[sp - 1]) + " <= " + m_info + " < "
			       + toStr(m_cutoff[sp]);
	} else {
		DBG_FAILIF(static_cast<UINT>(sp) >= m_values.size(), IndexError,
			"Virtual Subpoplation index out of range of 0 ~ "
			+ toStr(m_values.size() - 1));
		return m_info + " = " + toStr(m_values[sp]);
	}
}


proportionSplitter::proportionSplitter(vectorf const & proportions)
	: vspSplitter(), m_proportions(proportions)
{
	DBG_ASSERT(fcmp_eq(std::accumulate(proportions.begin(),
				proportions.end(), 0.), 1.), ValueError,
		"Passed proportions should add up to one");
}


ULONG proportionSplitter::size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	if (virtualSubPop == InvalidSubPopID)
		return countVisibleInds(pop, subPop);
	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_proportions.size(), IndexError,
		"Virtual subpopulation index out of range");
	//
	if (static_cast<UINT>(virtualSubPop) < m_proportions.size() - 1)
		return static_cast<ULONG>(pop.subPopSize(subPop) * m_proportions[virtualSubPop]);
	// to avoid floating point problem, the last subpop is specially treated
	ULONG size = pop.subPopSize(subPop);
	ULONG spSize = size;
	// virtualSubPop == m_proportions.size() - 1
	for (int i = 0; i < virtualSubPop; ++i)
		spSize -= static_cast<ULONG>(size * m_proportions[i]);
	return spSize;
}


UINT proportionSplitter::numVirtualSubPop()
{
	return m_proportions.size();
}


void proportionSplitter::activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
                                  activateType type)
{
	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_proportions.size(), IndexError,
		"Virtual subpopulation index out of range");

	ULONG size = pop.subPopSize(subPop);
	ULONG spSize = static_cast<ULONG>(size * m_proportions[0]);
	ULONG spCount = 0;
	SubPopID sp = 0;
	SubPopID visibleSP = virtualSubPop;

	RawIndIterator it = pop.rawIndBegin(subPop);
	RawIndIterator it_end = pop.rawIndEnd(subPop);
	for (; it != it_end; ++it, ++spCount) {
		if (spCount == spSize) {
			++sp;
			if (static_cast<UINT>(sp) + 1 < m_proportions.size())
				spSize = static_cast<ULONG>(size * m_proportions[sp]);
			else
				// in the last virtual sp.
				// use a bigger size to make sure the last few
				// individuals are counted.
				spSize = size;
			spCount = 0;
		}
		if (type == Visible)
			it->setVisible(sp == visibleSP);
		else
			it->setIteratable(sp == visibleSP);
	}
	if (type == Visible)
		m_activated = true;
}


void proportionSplitter::deactivate(population & pop, SubPopID subPop)
{
	resetSubPop(pop, subPop);
}


string proportionSplitter::name(SubPopID subPop)
{
	DBG_FAILIF(static_cast<UINT>(subPop) >= m_proportions.size(), IndexError,
		"Virtual subpopulation index out of range");
	return "Prop " + toStr(m_proportions[subPop]);
}


rangeSplitter::rangeSplitter(const intMatrix & ranges)
	: vspSplitter(), m_ranges(ranges)
{
	for (size_t i = 0; i < m_ranges.size(); ++i) {
		DBG_FAILIF(m_ranges[i].size() != 2
			|| m_ranges[i][0] > m_ranges[i][1],
			ValueError, "Wrong range [" +
			toStr(m_ranges[i][0]) + ", " +
			toStr(m_ranges[i][1]) + ")");
	}
}


ULONG rangeSplitter::size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	if (virtualSubPop == InvalidSubPopID)
		return countVisibleInds(pop, subPop);
	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_ranges.size(), IndexError,
		"Virtual subpopulation index out of range");
	if (static_cast<UINT>(m_ranges[virtualSubPop][0]) > pop.subPopSize(subPop))
		return 0;
	if (static_cast<UINT>(m_ranges[virtualSubPop][1]) > pop.subPopSize(subPop))
		return pop.subPopSize(subPop) - m_ranges[virtualSubPop][0];
	return m_ranges[virtualSubPop][1] - m_ranges[virtualSubPop][0];
}


UINT rangeSplitter::numVirtualSubPop()
{
	return m_ranges.size();
}


void rangeSplitter::activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
                             activateType type)
{
	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_ranges.size(), IndexError,
		"Virtual subpopulation index out of range");

	ULONG low = m_ranges[virtualSubPop][0];
	ULONG high = m_ranges[virtualSubPop][1];
	ULONG idx = 0;

	RawIndIterator it = pop.rawIndBegin(subPop);
	RawIndIterator it_end = pop.rawIndEnd(subPop);
	for (; it != it_end; ++it, ++idx)
		if (type == Visible)
			it->setVisible(idx >= low && idx < high);
		else
			it->setIteratable(idx >= low && idx < high);
	if (type == Visible)
		m_activated = true;
}


void rangeSplitter::deactivate(population & pop, SubPopID subPop)
{
	resetSubPop(pop, subPop);
}


string rangeSplitter::name(SubPopID subPop)
{
	DBG_FAILIF(static_cast<UINT>(subPop) >= m_ranges.size(), IndexError,
		"Virtual subpopulation index out of range");
	return "Range [" + toStr(m_ranges[subPop][0]) + ", " +
	       toStr(m_ranges[subPop][1]) + ")";
}


genotypeSplitter::genotypeSplitter(const vectori & loci,
                                   const intMatrix & alleles,
                                   bool phase)
	: vspSplitter(), m_loci(loci), m_alleles(alleles),
	m_phase(phase)
{
}


ULONG genotypeSplitter::size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const
{
	if (virtualSubPop == InvalidSubPopID)
		return countVisibleInds(pop, subPop);
	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_alleles.size(), IndexError,
		"Virtual subpopulation index out of range");
	const vectori & alleles = m_alleles[virtualSubPop];
	ULONG count = 0;
	ConstRawIndIterator it = pop.rawIndBegin(subPop);
	ConstRawIndIterator it_end = pop.rawIndEnd(subPop);
	for (; it != it_end; ++it)
		if (match(& * it, alleles))
			++count;
	return count;
}


UINT genotypeSplitter::numVirtualSubPop()
{
	return m_alleles.size();
}


void genotypeSplitter::activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
                                activateType type)
{
	DBG_FAILIF(static_cast<UINT>(virtualSubPop) >= m_alleles.size(), IndexError,
		"Virtual subpopulation index out of genotype");

	const vectori & alleles = m_alleles[virtualSubPop];
	RawIndIterator it = pop.rawIndBegin(subPop);
	RawIndIterator it_end = pop.rawIndEnd(subPop);
	for (; it != it_end; ++it)
		if (type == Visible)
			it->setVisible(match(& * it, alleles));
		else
			it->setIteratable(match(& * it, alleles));
	if (type == Visible)
		m_activated = true;
}


void genotypeSplitter::deactivate(population & pop, SubPopID subPop)
{
	resetSubPop(pop, subPop);
}


//
// 1: 0 1
// 1: 0 1 1 1
// Two possibilities (depending on ploidy)
// 1, 2: 0 1 1 0 1 1 1 1
string genotypeSplitter::name(SubPopID subPop)
{
	DBG_FAILIF(static_cast<UINT>(subPop) >= m_alleles.size(), IndexError,
		"Virtual subpopulation index out of genotype");
	string label = "Genotype ";
	for (size_t i = 0; i < m_loci.size(); ++i) {
		if (i != 0)
			label += ", ";
		label += toStr(m_loci[i]);
	}
	label += ":";
	for (size_t i = 0; i < m_alleles[subPop].size(); ++i)
		label += " " + toStr(m_alleles[subPop][i]);
	return label;
}


bool genotypeSplitter::match(const individual * it, const vectori & alleles) const
{
	unsigned types = alleles.size() / it->ploidy() / m_loci.size();
	DBG_FAILIF(alleles.size() != types * it->ploidy() * m_loci.size(),
		ValueError, "Given genotype does not match population ploidy.");

	int ploidy = it->ploidy();
	if (types == 1)
		return matchSingle(it, alleles);
	for (unsigned t = 0; t < types; ++t) {
		vectori partial(alleles.begin() + t * ploidy * m_loci.size(),
			alleles.begin() + (t + 1) * ploidy * m_loci.size());
		if (matchSingle(it, partial))
			return true;
	}
	return false;
}

bool genotypeSplitter::matchSingle(const individual * it, const vectori & alleles) const
{
	int ploidy = it->ploidy();
	if (m_phase || ploidy == 1) {
		UINT idx = 0;
		vectori::const_iterator loc = m_loci.begin();
		vectori::const_iterator loc_end = m_loci.end();
		for (; loc != loc_end; ++loc)
			for (int p = 0; p < ploidy; ++p)
				if (it->allele(*loc, p) != alleles[idx++])
					return false;
		return true;
	} else if (ploidy == 2) {
		UINT idx = 0;
		vectori::const_iterator loc = m_loci.begin();
		vectori::const_iterator loc_end = m_loci.end();
		for (; loc != loc_end; ++loc) {
			Allele a1 = it->allele(*loc, 0);
			Allele a2 = it->allele(*loc, 1);
			if ((a1 == alleles[idx] && a2 == alleles[idx + 1]) ||
			    (a1 == alleles[idx + 1] && a2 == alleles[idx]))
				idx += 2;
			else
				return false;
		}
		return true;
	} else
		DBG_FAILIF(true, SystemError, "Ploidy>=3 is not implemented");
	return false;
}


}     // namespce simuPOP


