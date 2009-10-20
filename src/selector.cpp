/**
 *  $File: selector.cpp $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "selector.h"

namespace simuPOP {
bool selector::apply(population & pop)
{
	UINT fit_id = pop.infoIdx(this->infoField(0));

	subPopList subPops = applicableSubPops();
	// the usual whole population, easy case.
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);

	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for (; sp != spEnd; ++sp) {
		if (sp->isVirtual())
			pop.activateVirtualSubPop(*sp);
		IndIterator ind = pop.indIterator(sp->subPop());
		for (; ind.valid(); ++ind)
			ind->setInfo(indFitness(& * ind, pop.gen()), fit_id);
		pop.turnOnSelection(sp->subPop());
	}

	return true;
}


double mapSelector::indFitness(individual * ind, ULONG gen)
{
	size_t ply = ind->ploidy();

	vectori alleles(ply * m_loci.size());
	size_t idx = 0;

	vectoru::iterator loc = m_loci.begin();
       	vectoru::iterator locEnd = m_loci.end();
	for (; loc != locEnd; ++loc) {
		for (size_t p = 0; p < ply; ++p, ++idx)
			alleles[idx] = ind->allele(*loc, p);
	}

	tupleDict::iterator pos = m_dict.find(alleles);
	
	if (pos != m_dict.end())
		return pos->second;

	if (ply > 1) {
		// try to look up the key without phase
		tupleDict::iterator it = m_dict.begin();
		tupleDict::iterator itEnd = m_dict.end();
		for (; it != itEnd; ++it) {
			bool ok = true;
			const tupleDict::key_type & key = it->first;
			for (size_t i = 0; i < m_loci.size(); ++i) {
				if (ply == 2) {
					if ((alleles[2 * i ] != key[0] || alleles[2 * i + 1] != key[1]) &&
						(alleles[2 * i ] != key[1] || alleles[2 * i + 1] != key[0])) {
						ok = false;
						break;
					}
				} else {
					std::sort(alleles.begin() + 2 * i, alleles.begin() + 2 * (i + 1));
					tupleDict::key_type sorted_key = it->first;
					std::sort(sorted_key.begin(), sorted_key.end());
					for (size_t j = 0; j < key.size(); ++j) {
						if (alleles[2 * i + j] != sorted_key[j]) {
							ok = false;
							break;
						}
					}
				}
					
			}
			if (ok)
				return it->second;
		}
	}
	// no match
	string allele_string = "(";
	for (size_t i = 0; i < alleles.size(); ++i) {
		if (i != 0)
			allele_string += ", ";
		allele_string += toStr(alleles[i]);
	}
	allele_string += ")";
	DBG_ASSERT(false, ValueError, "No fitness value for genotype " + allele_string);
	// this line should not be reached.
	return 0;
}


// currently assuming diploid
double maSelector::indFitness(individual * ind, ULONG gen)
{
	UINT index = 0;
	bool singleST = m_wildtype.size() == 1;

	for (vectoru::iterator loc = m_loci.begin(); loc != m_loci.end(); ++loc) {
		// get genotype of ind
		Allele a = ToAllele(ind->allele(*loc, 0));
		Allele b = ToAllele(ind->allele(*loc, 1));

		int numWildtype = 0;

		// count number of wildtype
		// this improve the performance a little bit
		if (singleST) {
			numWildtype = (AlleleUnsigned(a) == m_wildtype[0]) + (AlleleUnsigned(b) == m_wildtype[0]);
		} else {
			if (find(m_wildtype.begin(), m_wildtype.end(), AlleleUnsigned(a)) != m_wildtype.end() )
				numWildtype++;

			if (find(m_wildtype.begin(), m_wildtype.end(), AlleleUnsigned(b)) != m_wildtype.end() )
				numWildtype++;
		}

		index = index * 3 + 2 - numWildtype;
	}
	return m_fitness[index];
}


double mlSelector::indFitness(individual * ind, ULONG gen)
{
	if (m_mode == Multiplicative) {
		double fit = 1;
		for (opList::iterator s = m_selectors.begin(), sEnd = m_selectors.end();
		     s != sEnd; ++s)
			fit *= static_cast<selector * >(*s)->indFitness(ind, gen);
		return fit;
	} else if (m_mode == Additive) {
		double fit = 1;
		for (opList::iterator s = m_selectors.begin(), sEnd = m_selectors.end();
		     s != sEnd; ++s)
			fit -= 1 - static_cast<selector * >(*s)->indFitness(ind, gen);
		return fit < 0 ? 0. : fit;
	} else if (m_mode == Heterogeneity) {
		double fit = 1;
		for (opList::iterator s = m_selectors.begin(), sEnd = m_selectors.end();
		     s != sEnd; ++s)
			fit *= 1 - static_cast<selector * >(*s)->indFitness(ind, gen);
		return 1 - fit;
	}
	// this is the case for none.
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
			m_alleles[j++] = ToAllele(ind->allele(m_loci[i], p));

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
		return m_func(PyObj_As_Double, "(Oi)", m_numArray, gen);
	else
		return m_func(PyObj_As_Double, "(OiO)", m_numArray, gen, m_infoArray);
	return 0.;
}


}
