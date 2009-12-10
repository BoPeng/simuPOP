/**
 *  $File: initializer.cpp $
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

#include "initializer.h"

namespace simuPOP {

string initSex::describe(bool format)
{
	string desc = "<simuPOP.initSex> initialize sex ";

	if (!m_sex.empty())
		desc += "using a list";
	else if (m_maleProp < 0) {
		desc += "randomly";
		if (m_maleFreq != 0.5)
			desc += " with probability " + toStr(m_maleFreq) + " being a male";
	} else
		desc += "randomly with " + toStr(m_maleProp * 100) + " percent of males";
	return desc;
}


bool initSex::apply(population & pop)
{
	subPopList subPops = applicableSubPops();

	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);

	size_t idx = 0;
	subPopList::iterator sp = subPops.begin();
	subPopList::iterator sp_end = subPops.end();
	for (; sp != sp_end; ++sp) {
		weightedSampler ws(GetRNG());
		if (m_maleProp >= 0) {
			vectorf prop(2, m_maleProp);
			prop[1] = 1 - prop[0];
			ws.set(prop, pop.subPopSize(*sp));
		}
		pop.activateVirtualSubPop(*sp);
		IndIterator ind = pop.indIterator(sp->subPop());
		size_t sexSz = m_sex.size();
		if (!m_sex.empty())
			for (; ind.valid(); ++ind, ++idx)
				ind->setSex(m_sex[idx % sexSz] == 1 ? MALE : FEMALE);
		else if (m_maleProp < 0)
			for (; ind.valid(); ++ind)
				ind->setSex(GetRNG().randUniform() < m_maleFreq ? MALE : FEMALE);
		else
			for (; ind.valid(); ++ind)
				ind->setSex(ws.get() == 0 ? MALE : FEMALE);
		pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


string initInfo::describe(bool format)
{
	string desc = "<simuPOP.initInfo> initialize information field";

	if (infoSize() > 1)
		desc += "s";
	for (size_t i = 0; i < infoSize(); ++i) {
		desc += i == 0 ? " " : ", ";
		desc += infoField(i);
	}
	if (m_values.empty()) {
		PyObject * name = PyObject_GetAttrString(m_values.func().func(), "__name__");
		DBG_FAILIF(name == NULL, RuntimeError, "Passwd object does not have attribute __name__");
		desc += " using a Python function " + string(PyString_AsString(name));
	} else
		desc += " using a list of values";
	return desc;
}


bool initInfo::apply(population & pop)
{
	vectoru infoIdx(infoSize());

	if (infoIdx.empty())
		return true;

	for (size_t i = 0; i < infoIdx.size(); ++i)
		infoIdx[i] = pop.infoIdx(infoField(i));

	subPopList subPops = applicableSubPops();

	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);

	subPopList::iterator sp = subPops.begin();
	subPopList::iterator sp_end = subPops.end();
	size_t idx = 0;
	const vectorf & values = m_values.elems();

	for (; sp != sp_end; ++sp) {
		pop.activateVirtualSubPop(*sp);
		IndIterator ind = pop.indIterator(sp->subPop());
		size_t numValues = m_values.size();
		for (; ind.valid(); ++ind, ++idx) {
			for (size_t i = 0; i < infoIdx.size(); ++i) {
				if (values.empty())
					ind->setInfo(m_values.func() (PyObj_As_Double, "()"), infoIdx[i]);
				else
					ind->setInfo(values[idx % numValues], infoIdx[i]);
			}
		}
		pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


initByFreq::initByFreq(const matrix & alleleFreq, const uintList & loci,
	const uintList & ploidy, bool identicalInds,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops,
	const stringList & infoFields)
	: baseOperator("", begin, end, step, at, reps, subPops, infoFields),
	m_alleleFreq(alleleFreq), m_identicalInds(identicalInds),
	m_loci(loci), m_ploidy(ploidy)
{

	DBG_FAILIF(m_alleleFreq.empty(), IndexError, "Allele frequency list is empty.");

	for (size_t i = 0; i < m_alleleFreq.size(); ++i) {
		for (size_t j = 0; j < m_alleleFreq[i].size(); ++j)
			if (fcmp_lt(m_alleleFreq[i][j], 0) || fcmp_gt(m_alleleFreq[i][j], 1))
				throw ValueError("Allele frequency should be between 0 and 1");
		if (fcmp_ne(accumulate(m_alleleFreq[i].begin(), m_alleleFreq[i].end(), 0.), 1.0))
			throw ValueError("Allele frequencies should add up to one.");
	}
}


string initByFreq::describe(bool format)
{
	return "<simuPOP.initByFreq> initialize individual genotype using a list of allele frequencies";
}


bool initByFreq::apply(population & pop)
{
	subPopList subPops = applicableSubPops();

	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);

	DBG_FAILIF(m_alleleFreq.size() > 1 && m_alleleFreq.size() != subPops.size(),
		ValueError, "Ranges and values should have the same length");

	vectoru loci = m_loci.elems();
	if (m_loci.allAvail())
		for (size_t i = 0 ; i < pop.totNumLoci(); ++i)
			loci.push_back(i);

	vectoru ploidy = m_ploidy.elems();
	if (m_ploidy.allAvail())
		for (size_t i = 0 ; i < pop.ploidy(); ++i)
			ploidy.push_back(i);

	pop.sortIndividuals();

	subPopList::iterator sp = subPops.begin();
	subPopList::iterator sp_end = subPops.end();
	for (size_t idx = 0; sp != sp_end; ++sp, ++idx) {
		//
		vectorf & alleleFreq = m_alleleFreq.size() == 1 ? m_alleleFreq[0] : m_alleleFreq[idx];

		// weightedSampler ws(GetRNG(), incFreq);
		weightedSampler ws(GetRNG(), alleleFreq);

		DBG_ASSERT(fcmp_eq(std::accumulate(alleleFreq.begin(), alleleFreq.end(), 0.), 1),
			SystemError, "Allele frequecies should add up to one.");

		// will go through virtual subpopulation if sp is virtual
		pop.activateVirtualSubPop(*sp);
		IndIterator left = pop.indIterator(sp->subPop());
		IndIterator it = pop.indIterator(sp->subPop());
		for (; it.valid(); ++it) {
			if (!m_identicalInds || it == left) {
				for (vectoru::iterator loc = loci.begin(); loc != loci.end(); ++loc)
					for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
						it->setAllele(ToAllele(ws.get()), *loc, *p);
			} else {
				// identical individuals
				for (vectoru::iterator loc = loci.begin(); loc != loci.end(); ++loc)
					for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
						it->setAllele(ToAllele(left->allele(*loc, *p)), *loc, *p);
			}
		}
		pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


initByValue::initByValue(intMatrix value, const uintList & loci, const uintList & ploidy,
	const floatList & proportions, const floatList & frequencies,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops,
	const stringList & infoFields)
	: baseOperator("", begin, end, step, at, reps, subPops, infoFields),
	m_value(value), m_proportion(proportions.elems()), m_frequencies(frequencies.elems()),
	m_loci(loci), m_ploidy(ploidy)
{
	DBG_FAILIF(m_value.empty(), ValueError,
		"Please specify an array of alleles in the order of chrom_1...chrom_n for all copies of chromosomes");

	DBG_FAILIF(!m_proportion.empty() && fcmp_ne(accumulate(m_proportion.begin(), m_proportion.end(), 0.0), 1),
		ValueError, "Proportion should add up to one.");

	DBG_FAILIF(!m_frequencies.empty() && fcmp_ne(accumulate(m_frequencies.begin(), m_frequencies.end(), 0.0), 1),
		ValueError, "Frequencies should add up to one.");

	DBG_FAILIF(!m_frequencies.empty() && !m_proportion.empty(), ValueError,
		"Please specify only one of parameters frequencies and proportions.");
}


bool initByValue::apply(population & pop)
{
#ifndef OPTIMIZED
	UINT gSz = m_value[0].size();
	for (size_t v = 1; v < m_value.size(); ++v)
		if (m_value[v].size() != gSz)
			throw ValueError("Given values should have the same length (either one copy of chromosomes or the whole genotype.");
#endif

	subPopList subPops = applicableSubPops();
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);

	vectoru loci = m_loci.elems();
	if (m_loci.allAvail())
		for (size_t i = 0 ; i < pop.totNumLoci(); ++i)
			loci.push_back(i);

	vectoru ploidy = m_ploidy.elems();
	if (m_ploidy.allAvail())
		for (size_t i = 0 ; i < pop.ploidy(); ++i)
			ploidy.push_back(i);

	pop.sortIndividuals();

	DBG_FAILIF(!m_proportion.empty() && m_proportion.size() != m_value.size(), ValueError,
		"If proportions are given, its length should match that of values.");

	DBG_FAILIF(!m_frequencies.empty() && m_frequencies.size() != m_value.size(), ValueError,
		"If frequenciess are given, its length should match that of values.");

	DBG_FAILIF(m_value.size() != 1 && m_value.size() != m_proportion.size()
		&& m_value.size() != m_frequencies.size() && m_value.size() != subPops.size(), ValueError,
		"If mutliple values are given, its length should match proportion or (virtual) subpopulations");

	DBG_FAILIF(m_value[0].size() != loci.size() && m_value[0].size() != loci.size() * ploidy.size(),
		ValueError, "Size of given value should be mutiples of number of loci.");

	subPopList::iterator sp = subPops.begin();
	subPopList::iterator sp_end = subPops.end();
	for (size_t idx = 0; sp != sp_end; ++sp, ++idx) {
		//
		weightedSampler ws(GetRNG());
		if (!m_proportion.empty())
			ws.set(m_proportion, pop.subPopSize(*sp) * ploidy.size());
		if (!m_frequencies.empty())
			ws.set(m_frequencies);

		pop.activateVirtualSubPop(*sp);
		// will go through virtual subpopulation if sp is virtual
		IndIterator it = pop.indIterator(sp->subPop());
		for (; it.valid(); ++it) {
			if (m_value[0].size() == loci.size()) { // for each ploidy
				for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p) {
					vectori & value = m_proportion.empty() && m_frequencies.empty() ?
					                  (m_value.size() == 1 ? m_value[0] : m_value[idx]) : m_value[ws.get()];
					for (size_t i = 0; i < value.size(); ++i)
						it->setAllele(ToAllele(value[i]), loci[i], *p);
				}
			} else {
				// (m_value[0].size() == loci.size() * ploidy.size())
				vectori & value = m_proportion.empty() && m_frequencies.empty() ?
				                  (m_value.size() == 1 ? m_value[0] : m_value[idx]) : m_value[ws.get()];
				size_t i = 0;
				for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
					for (vectoru::iterator loc = loci.begin(); loc != loci.end(); ++loc, ++i)
						it->setAllele(ToAllele(value[i]), *loc, *p);
			}
		}
		pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


}
