/**
 *  $File: penetrance.cpp $
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

#include "penetrance.h"

namespace simuPOP {
// set pentrance to all individuals and record penetrance if requested.
bool basePenetrance::apply(population & pop)
{
	bool savePene = infoSize() > 0;
	UINT infoIdx = 0;

	if (savePene)
		infoIdx = pop.infoIdx(infoField(0));

	UINT ansGen = 0;
	UINT oldGen = pop.curAncestralGen();
	if (m_ancGen == -1)
		ansGen = pop.ancestralGens();
	else if (m_ancGen > 0) {
		if (static_cast<UINT>(m_ancGen) > pop.ancestralGens())
			ansGen = pop.ancestralGens();
		else
			ansGen = m_ancGen;
	}
	for (UINT i = 0; i <= ansGen; ++i) {
		pop.useAncestralGen(i);

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
			for (; ind.valid(); ++ind) {
				double p = penet(& * ind, pop.gen());

				if (savePene)
					ind->setInfo(p, infoIdx);

				if (GetRNG().randUniform() < p)
					ind->setAffected(true);
				else
					ind->setAffected(false);
			}
			if (sp->isVirtual())
				pop.deactivateVirtualSubPop(sp->subPop());
		}
	}
	pop.useAncestralGen(oldGen);

	return true;
}


bool basePenetrance::applyDuringMating(population & pop, RawIndIterator offspring,
                                       individual * dad, individual * mom)
{
	double p = penet(& * offspring, pop.gen());

	if (infoSize() > 0)
		offspring->setInfo(p, infoField(0));
	if (GetRNG().randUniform() < p)
		offspring->setAffected(true);
	else
		offspring->setAffected(false);
	return true;
}


// this function is the same as mapSelector.
double mapPenetrance::penet(individual * ind, ULONG gen)
{
	vectoru chromTypes;

	for (size_t i = 0; i < m_loci.size(); ++i)
		chromTypes.push_back(ind->chromType(ind->chromLocusPair(m_loci[i]).first));

	size_t ply = ind->ploidy();
	if (ind->isHaplodiploid() && ind->sex() == Male)
		ply = 1;

	vectori alleles;
	alleles.reserve(ply * m_loci.size());

	for (size_t idx = 0; idx < m_loci.size(); ++idx) {
		for (size_t p = 0; p < ply; ++p) {
			if (chromTypes[idx] == ChromosomeY && ind->sex() == Female)
				continue;
			if (((chromTypes[idx] == ChromosomeX && p == 1) ||
			     (chromTypes[idx] == ChromosomeY && p == 0)) && ind->sex() == Male)
				continue;
			alleles.push_back(ind->allele(m_loci[idx], p));
		}
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
			UINT begin_idx = 0;
			UINT end_idx = 0;
			for (size_t i = 0; i < m_loci.size(); ++i) {
				if (chromTypes[i] == ChromosomeY) {
					if (ind->sex() == Female)
						continue;
					else
						++end_idx;
				} else if (chromTypes[i] == ChromosomeX && ind->sex() == Male)
					++end_idx;
				else
					end_idx += ply;
				if (key.size() != end_idx - begin_idx) {
					ok = false;
					break;
				}
				if (ply == 2) {
					if ((alleles[begin_idx] != key[0] || alleles[end_idx - 1] != key[1]) &&
					    (alleles[begin_idx] != key[1] || alleles[end_idx - 1] != key[0])) {
						ok = false;
						break;
					}
				} else {
					std::sort(alleles.begin() + begin_idx, alleles.begin() + end_idx);
					tupleDict::key_type sorted_key = it->first;
					std::sort(sorted_key.begin(), sorted_key.end());
					for (size_t j = 0; j < sorted_key.size(); ++j) {
						if (alleles[ply * i + j] != sorted_key[j]) {
							ok = false;
							break;
						}
					}
				}
				begin_idx = end_idx;
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
	DBG_ASSERT(false, ValueError, "No penetrance value for genotype " + allele_string);
	// this line should not be reached.
	return 0;
}


string maPenetrance::describe(bool format)
{
	return "<simuPOP.maPenetrance> multiple-alleles penetrance" ;
}


// this function is the same as maSelector.
double maPenetrance::penet(individual * ind, ULONG gen)
{
	UINT index = 0;
	bool singleST = m_wildtype.size() == 1;

	DBG_FAILIF((ind->ploidy() == 2 && m_penetrance.size() != static_cast<UINT>(pow(3., static_cast<double>(m_loci.size())))) ||
		(ind->ploidy() == 1 && m_penetrance.size() != static_cast<UINT>(pow(2., static_cast<double>(m_loci.size())))),
		ValueError, "Please specify penetrance for each combination of genotype.");

	for (vectoru::iterator loc = m_loci.begin(); loc != m_loci.end(); ++loc) {
		if (ind->ploidy() == 1) {
			Allele a = ToAllele(ind->allele(*loc));
			if (singleST)
				index = index * 2 + (AlleleUnsigned(a) != m_wildtype[0]);
			else
				index = index * 2 + (find(m_wildtype.begin(), m_wildtype.end(), AlleleUnsigned(a)) == m_wildtype.end());
		} else if (ind->ploidy() == 2) {
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
		} else {
			DBG_FAILIF(true, ValueError, "The maSelector only supports haploid and diploid populations.");
		}
	}
	return m_penetrance[index];
}


double mlPenetrance::penet(individual * ind, ULONG gen)
{
	if (m_mode == Multiplicative) {
		// x1 x2 x3 ...
		double pen = 1;
		for (vectorop::iterator s = m_peneOps.begin(), sEnd = m_peneOps.end();
		     s != sEnd; ++s)
			pen *= static_cast<basePenetrance *>(*s)->penet(ind, gen);
		return pen;
	} else if (m_mode == Additive) {
		// x1 + x2 + x3
		double pen = 0;
		for (vectorop::iterator s = m_peneOps.begin(), sEnd = m_peneOps.end();
		     s != sEnd; ++s)
			pen += static_cast<basePenetrance *>(*s)->penet(ind, gen);
		return pen > 1 ? 1 : pen;
	} else if (m_mode == Heterogeneity) {
		// 1-(1-x1)(1-x2)
		double pen = 1;
		for (vectorop::iterator s = m_peneOps.begin(), sEnd = m_peneOps.end();
		     s != sEnd; ++s)
			pen *= 1 - static_cast<basePenetrance *>(*s)->penet(ind, gen);
		return pen > 1 ? 0 : 1 - pen;
	}

	return 0.0;
}


// the same as pySelector
double pyPenetrance::penet(individual * ind, ULONG gen)
{
	vectoru chromTypes;

	for (size_t i = 0; i < m_loci.size(); ++i)
		chromTypes.push_back(ind->chromType(ind->chromLocusPair(m_loci[i]).first));

	size_t ply = ind->ploidy();
	if (ind->isHaplodiploid() && ind->sex() == Male)
		ply = 1;

	vectori alleles;
	alleles.reserve(ply * m_loci.size());

	for (size_t idx = 0; idx < m_loci.size(); ++idx) {
		for (size_t p = 0; p < ply; ++p) {
			if (chromTypes[idx] == ChromosomeY && ind->sex() == Female)
				continue;
			if (((chromTypes[idx] == ChromosomeX && p == 1) ||
			     (chromTypes[idx] == ChromosomeY && p == 0)) && ind->sex() == Male)
				continue;
			alleles.push_back(ind->allele(m_loci[idx], p));
		}
	}

	vectorf info(m_paramFields.size());
	for (size_t i = 0; i < m_paramFields.size(); ++i)
		info[i] = ind->info(m_paramFields[i]);

	//
	if (m_genotype == NULL || static_cast<UINT>(PySequence_Size(m_genotype)) != alleles.size()) {
		Py_XDECREF(m_genotype);
		m_genotype = PyTuple_New(alleles.size());
	}
	if (!info.empty() && (m_info == NULL || static_cast<UINT>(PySequence_Size(m_info)) != info.size())) {
		Py_XDECREF(m_info);
		m_info = PyTuple_New(info.size());
	}
	// set value
	for (size_t i = 0; i < alleles.size(); ++i)
		PyTuple_SetItem(m_genotype, i, PyInt_FromLong(alleles[i]));
	for (size_t i = 0; i < info.size(); ++i)
		PyTuple_SetItem(m_info, i, PyFloat_FromDouble(info[i]));

	double penetrance;
	if (info.empty())
		penetrance = m_func(PyObj_As_Double, "(Oi)", m_genotype, gen);
	else
		penetrance = m_func(PyObj_As_Double, "(OOi)", m_genotype, m_info, gen);

	DBG_DO(DBG_SELECTOR, cerr << "Genotype " << alleles << " info " << info << " penetrance " << penetrance << endl);
	return penetrance;

}


}


