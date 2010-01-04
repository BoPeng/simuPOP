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

string InitSex::describe(bool format)
{
	string desc = "<simuPOP.InitSex> initialize sex ";

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


bool InitSex::apply(Population & pop)
{
	subPopList subPops = applicableSubPops();

	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);

	size_t idx = 0;
	subPopList::iterator sp = subPops.begin();
	subPopList::iterator sp_end = subPops.end();
	for (; sp != sp_end; ++sp) {
		Weightedsampler ws(getRNG());
		if (m_maleProp >= 0) {
			vectorf prop(2, m_maleProp);
			prop[1] = 1 - prop[0];
			ws.set(prop, pop.subPopSize(*sp));
		} else if (m_sex.empty()) {
			vectorf prop(2, m_maleFreq);
			prop[1] = 1 - prop[0];
			ws.set(prop);
		}
		pop.activateVirtualSubPop(*sp);
		IndIterator ind = pop.indIterator(sp->subPop());
		size_t sexSz = m_sex.size();
		if (!m_sex.empty())
			for (; ind.valid(); ++ind, ++idx)
				ind->setSex(m_sex[idx % sexSz] == 1 ? MALE : FEMALE);
		else
			for (; ind.valid(); ++ind)
				ind->setSex(ws.get() == 0 ? MALE : FEMALE);
		pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


string InitInfo::describe(bool format)
{
	string desc = "<simuPOP.InitInfo> initialize information field";

	if (infoSize() > 1)
		desc += "s";
	for (size_t i = 0; i < infoSize(); ++i) {
		desc += i == 0 ? " " : ", ";
		desc += infoField(i);
	}
	if (m_values.empty())
		desc += " using a Python function " + m_values.func().name();
	else
		desc += " using a list of values";
	return desc;
}


bool InitInfo::apply(Population & pop)
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


InitGenotype::InitGenotype(const vectorf & freq,
	const uintList & genotype, const vectorf & prop,
	const uintList & loci,
	const uintList & ploidy,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops,
	const stringList & infoFields)
	: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
	m_freq(freq), m_genotype(genotype.elems()), m_prop(prop),
	m_loci(loci), m_ploidy(ploidy)
{
	for (size_t i = 0; i < m_freq.size(); ++i)
		DBG_FAILIF(fcmp_lt(m_freq[i], 0) || fcmp_gt(m_freq[i], 1), ValueError,
			"Allele frequency should be between 0 and 1");
	for (size_t i = 0; i < m_prop.size(); ++i)
		DBG_FAILIF(fcmp_lt(m_prop[i], 0) || fcmp_gt(m_prop[i], 1), ValueError,
			"Allele proportion should be between 0 and 1");
	DBG_FAILIF(!m_freq.empty() && fcmp_ne(accumulate(m_freq.begin(), m_freq.end(), 0.), 1.0), ValueError,
		"Allele frequencies should add up to one.");
	DBG_FAILIF(!m_prop.empty() && fcmp_ne(accumulate(m_prop.begin(), m_prop.end(), 0.), 1.0), ValueError,
		"Allele proportions should add up to one.");
	DBG_ASSERT(!m_freq.empty() + !m_genotype.empty() + !m_prop.empty() == 1, ValueError,
		"Please specify one and only one of parameters freq, genotype and prop.");
}


string InitGenotype::describe(bool format)
{
	string desc = "<simuPOP.InitGenotype> initialize individual genotype ";

	if (!m_freq.empty())
		desc += "acccording to allele frequencies.";
	else if (m_prop.empty())
		desc += "according to proportion of alleles.";
	else
		desc += "using specified genotype.";
	return desc;
}


bool InitGenotype::apply(Population & pop)
{
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

	subPopList::iterator sp = subPops.begin();
	subPopList::iterator sp_end = subPops.end();
	size_t sz = m_genotype.size();
	for (size_t idx = 0; sp != sp_end; ++sp) {
		//
		Weightedsampler ws(getRNG());
		if (!m_prop.empty())
			ws.set(m_prop, pop.subPopSize(*sp));
		else if (!m_freq.empty())
			ws.set(m_freq);

		// will go through virtual subpopulation if sp is virtual
		pop.activateVirtualSubPop(*sp);
		IndIterator it = pop.indIterator(sp->subPop());
		if (!m_genotype.empty()) {
			for (; it.valid(); ++it)
				for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
					for (vectoru::iterator loc = loci.begin(); loc != loci.end(); ++loc, ++idx)
						it->setAllele(ToAllele(m_genotype[idx % sz]), *loc, *p);
		} else
			for (; it.valid(); ++it) {
				for (vectoru::iterator loc = loci.begin(); loc != loci.end(); ++loc)
					for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
						it->setAllele(ToAllele(ws.get()), *loc, *p);
			}
		pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


}
