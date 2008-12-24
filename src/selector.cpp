/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu                                                        *
*                                                                         *
*   $LastChangedDate$
*   $Rev$
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

#include "selector.h"

namespace simuPOP {
bool selector::apply(population & pop)
{
	UINT fit_id = pop.infoIdx(this->infoField(0));

	if (m_subPops.empty()) {
		IndInfoIterator fitness = pop.infoBegin(fit_id);
		// fitness may change with generation so pass generation information
		for (IndIterator it = pop.indBegin(); it.valid(); ++it)
			*fitness++ = indFitness(& * it, pop.gen()) ;
		pop.turnOnSelection();
	} else {
		for (vectoru::iterator sp = m_subPops.begin(); sp != m_subPops.end(); ++sp) {
			IndInfoIterator fitness = pop.infoBegin(fit_id, *sp);
			DBG_FAILIF(*sp > pop.numSubPop(), IndexError,
				"Wrong subpopulation index" + toStr(*sp) + " (number of subpops is " +
				toStr(pop.numSubPop()) + ")");
			for (IndIterator it = pop.indBegin(*sp); it.valid(); ++it)
				*fitness++ = indFitness(& * it, pop.gen());
			pop.turnOnSelection(*sp);
		}
	}

	return true;
}


double mapSelector::indFitness(individual * ind, ULONG gen)
{
	string key;
	size_t ply = ind->ploidy();

	vector<int> alleles(ply);

	for (vectoru::iterator loc = m_loci.begin(); loc != m_loci.end(); ++loc) {
		if (loc != m_loci.begin() )
			key += '|';
		for (int p = 0; p < ply; ++p)
			alleles[p] = ind->allele(*loc, p);
		// if no phase, sort alleles...
		if (!m_phase && ply > 1) {
			if (ply == 2) {
				if (alleles[0] > alleles[1]) { // swap
					int tmp = alleles[0];
					alleles[0] = alleles[1];
					alleles[1] = tmp;
				}
			} else
				std::sort(alleles.begin(), alleles.end());
		}
		// get key
		if (ply == 1)
			key += toStr(alleles[0]);
		else if (ply == 2)
			key += toStr(alleles[0]) + "-" + toStr(alleles[1]);
		else {
			key += toStr(alleles[0]);
			for (size_t i = 1; i < ply; ++i)
				key += "-" + toStr(alleles[i]);
		}
	}

	strDict::iterator pos = m_dict.find(key);

	DBG_ASSERT(pos != m_dict.end(), ValueError,
		"No fitness value for genotype " + key);

	return pos->second;
}


// currently assuming diploid
double maSelector::indFitness(individual * ind, ULONG gen)
{
	UINT index = 0;
	bool singleST = m_wildtype.size() == 1;

	for (vectoru::iterator loc = m_loci.begin(); loc != m_loci.end(); ++loc) {
		// get genotype of ind
		Allele a = ind->allele(*loc, 0);
		Allele b = ind->allele(*loc, 1);

		int numWildtype = 0;

		// count number of wildtype
		// this improve the performance a little bit
		if (singleST) {
			numWildtype = (a == m_wildtype[0]) + (b == m_wildtype[0]);
		} else {
			if (find(m_wildtype.begin(), m_wildtype.end(), a) != m_wildtype.end() )
				numWildtype++;

			if (find(m_wildtype.begin(), m_wildtype.end(), b) != m_wildtype.end() )
				numWildtype++;
		}

		index = index * 3 + 2 - numWildtype;
	}
	return m_fitness[index];
}


double mlSelector::indFitness(individual * ind, ULONG gen)
{
	if (m_mode == SEL_Multiplicative) {
		double fit = 1;
		for (vectorop::iterator s = m_selectors.begin(), sEnd = m_selectors.end();
		     s != sEnd; ++s)
			fit *= static_cast<selector * >(*s)->indFitness(ind, gen);
		return fit;
	} else if (m_mode == SEL_Additive) {
		double fit = 1;
		for (vectorop::iterator s = m_selectors.begin(), sEnd = m_selectors.end();
		     s != sEnd; ++s)
			fit -= 1 - static_cast<selector * >(*s)->indFitness(ind, gen);
		return fit < 0 ? 0. : fit;
	} else if (m_mode == SEL_Heterogeneity) {
		double fit = 1;
		for (vectorop::iterator s = m_selectors.begin(), sEnd = m_selectors.end();
		     s != sEnd; ++s)
			fit *= 1 - static_cast<selector * >(*s)->indFitness(ind, gen);
		return 1 - fit;
	}
	// this is the case for SEL_none.
	return 1.0;
}


double pySelector::indFitness(individual * ind, ULONG gen)
{
	if (m_len == 0) {
		m_len = m_loci.size() * ind->ploidy();
		m_alleles.resize(m_len);
		m_numArray = Allele_Vec_As_NumArray(m_alleles.begin(), m_alleles.end() );
	}

	DBG_FAILIF(static_cast<size_t>(m_len) != ind->ploidy() * m_loci.size(),
		SystemError,
		"Length of m_len is wrong. Have you changed pop type?");

	UINT pEnd = ind->ploidy();
	for (size_t i = 0, iEnd = m_loci.size(), j = 0; i < iEnd; ++i)
		for (UINT p = 0; p < pEnd; ++p)
			m_alleles[j++] = ind->allele(m_loci[i], p);

	if (infoSize() > 1) {
		if (m_info.size() + 1 != infoSize()) {
			m_info.resize(infoSize() - 1);
			m_infoArray = Double_Vec_As_NumArray(m_info.begin(), m_info.end());
		}
		// assign information fields from individusl
		for (size_t i = 1; i < infoSize(); ++i)
			m_info[i - 1] = ind->info(infoField(i));
	}


	if (infoSize() <= 1)
		return m_func.call("(Oi)", m_numArray, gen, PyObj_As_Double);
	else
		return m_func.call("(OiO)", m_numArray, gen, m_infoArray, PyObj_As_Double);
	return 0.;
}


}
