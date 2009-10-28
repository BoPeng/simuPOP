/**
 *  $File: qtrait.cpp $
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

#include "qtrait.h"

namespace simuPOP {

bool quanTrait::apply(population & pop)
{
	vectoru infoIdx(infoSize());

	for (size_t i = 0; i < infoSize(); ++i)
		infoIdx[i] = pop.infoIdx(infoField(i));

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
	vectorf traits(infoSize());
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
				qtrait(& * ind, pop.gen(), traits);
				for (size_t i = 0; i < infoSize(); ++i)
					ind->setInfo(traits[i], infoIdx[i]);
			}

			if (sp->isVirtual())
				pop.deactivateVirtualSubPop(sp->subPop());
		}
	}
	pop.useAncestralGen(oldGen);
	return true;
}


void pyQuanTrait::qtrait(individual * ind, ULONG gen, vectorf & traits)
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

	PyObject * res = NULL;
	if (info.empty())
		res = m_func("(Oi)", m_genotype, gen);
	else
		res = m_func("(OOi)", m_genotype, m_info, gen);

	if (PyNumber_Check(res)) {
		DBG_ASSERT(infoSize() == 1, RuntimeError,
			"A number is returned from a user-defined function but a sequence is expected.");
		PyObj_As_Double(res, traits[0]);
		Py_DECREF(res);
	} else if (PySequence_Check(res)) {
		DBG_ASSERT(PySequence_Size(res) == traits.size(), RuntimeError,
			"Length of returned sequence does not match number of trait fields");
		PyObj_As_Array(res, traits);
		Py_DECREF(res);
	} else {
		DBG_FAILIF(true, RuntimeError, "Invalid return value");
	}
	return;
}


}
