/**
 *  $File: initializer.cpp $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
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

string InitSex::describe(bool /* format */) const
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


bool InitSex::apply(Population & pop) const
{
	const subPopList subPops = applicableSubPops(pop);

	size_t idx = 0;
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator sp_end = subPops.end();

	for (; sp != sp_end; ++sp) {
		WeightedSampler ws;
		if (m_maleProp >= 0) {
			vectorf prop(2, m_maleProp);
			prop[1] = 1 - prop[0];
			ws.set(prop.begin(), prop.end(), pop.subPopSize(*sp));
		} else if (m_sex.empty()) {
			vectorf prop(2, m_maleFreq);
			prop[1] = 1 - prop[0];
			ws.set(prop.begin(), prop.end());
		}
		pop.activateVirtualSubPop(*sp);
		IndIterator ind = pop.indIterator(sp->subPop());
		size_t sexSz = m_sex.size();
		if (!m_sex.empty())
			for (; ind.valid(); ++ind, ++idx)
				ind->setSex(m_sex[idx % sexSz] == 1 ? MALE : FEMALE);
		else
			for (; ind.valid(); ++ind)
				ind->setSex(ws.draw() == 0 ? MALE : FEMALE);
		pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


string InitInfo::describe(bool /* format */) const
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


bool InitInfo::apply(Population & pop) const
{
	vectoru infoIdx(infoSize());

	if (infoIdx.empty())
		return true;

	for (size_t i = 0; i < infoIdx.size(); ++i)
		infoIdx[i] = pop.infoIdx(infoField(i));

	const subPopList subPops = applicableSubPops(pop);
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator sp_end = subPops.end();
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
	const intMatrix & haplotypes,
	const lociList & loci,
	const uintList & ploidy,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops,
	const stringList & infoFields)
	: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
	m_freq(freq), m_genotype(genotype.elems()), m_prop(prop),
	m_haplotypes(haplotypes.elems()), m_loci(loci), m_ploidy(ploidy)
{
	for (size_t i = 0; i < m_freq.size(); ++i)
		PARAM_FAILIF(fcmp_lt(m_freq[i], 0) || fcmp_gt(m_freq[i], 1), ValueError,
			"Allele frequency should be between 0 and 1");
	for (size_t i = 0; i < m_prop.size(); ++i)
		PARAM_FAILIF(fcmp_lt(m_prop[i], 0) || fcmp_gt(m_prop[i], 1), ValueError,
			"Allele proportion should be between 0 and 1");
	PARAM_FAILIF(!m_haplotypes.empty() && !m_freq.empty() && m_haplotypes.size() != m_freq.size(),
		ValueError, "Haplotype frequency, if specified, should be specified for each haplotype.")
	PARAM_FAILIF(!m_haplotypes.empty() && !m_prop.empty() && m_haplotypes.size() != m_prop.size(),
		ValueError, "Haplotype proportion, if specified, should be specified for each haplotype.")
	PARAM_FAILIF(!m_freq.empty() && fcmp_ne(accumulate(m_freq.begin(), m_freq.end(), 0.), 1.0), ValueError,
		"Allele frequencies should add up to one.");
	PARAM_FAILIF(!m_prop.empty() && fcmp_ne(accumulate(m_prop.begin(), m_prop.end(), 0.), 1.0), ValueError,
		"Allele proportions should add up to one.");
	PARAM_ASSERT(!m_haplotypes.empty() || (!m_freq.empty() + !m_genotype.empty() + !m_prop.empty() == 1), ValueError,
		"Please specify one and only one of parameters freq, prop (or haplotypes) and genotype.");
}


string InitGenotype::describe(bool /* format */) const
{
	string desc = "<simuPOP.InitGenotype> initialize individual genotype ";

	if (!m_freq.empty())
		desc += string("acccording to ") + (m_haplotypes.empty() ? " allele " : " haplotype ") + "frequencies.";
	else if (m_prop.empty())
		desc += string("according to proportion of ") + (m_haplotypes.empty() ? " alleles." : " haplotypes.");
	else
		desc += "using specified genotype.";
	return desc;
}


bool InitGenotype::apply(Population & pop) const
{
	const subPopList subPops = applicableSubPops(pop);

	const vectoru & loci = m_loci.elems(&pop);

	for (size_t i = 0; i < m_haplotypes.size(); ++i) {
		DBG_WARNIF(m_haplotypes[i].size() != loci.size(),
			"Haplotype [" + shorten(toStr(m_haplotypes[i])) + "] specified in operator InitGenotype has "
			+ toStr(m_haplotypes[i].size()) + " alleles but " + toStr(loci.size())
			+ " alleles are needed. This haplotype will be truncated or repeated.");
	}
	vectoru ploidy = m_ploidy.elems();
	if (m_ploidy.allAvail())
		for (size_t i = 0 ; i < pop.ploidy(); ++i)
			ploidy.push_back(i);

	pop.syncIndPointers();

	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator sp_end = subPops.end();
	size_t sz = m_genotype.size();

	DBG_WARNIF(sz >= 1 && sz != loci.size() && sz != loci.size() * pop.ploidy(),
		"Genotype [" + shorten(toStr(m_genotype)) + "} specified in operator InitGenotype has "
		+ toStr(m_genotype.size()) + " alleles, which does not match number of loci of individuals. "
		                             "This sequence will be truncated or repeated and may lead to erroneous results.");
	for (size_t idx = 0; sp != sp_end; ++sp) {
		// will go through virtual subpopulation if sp is virtual
		pop.activateVirtualSubPop(*sp);
		if (!m_genotype.empty()) {
#pragma omp parallel firstprivate(idx) if(numThreads() > 1)
			{
#ifdef _OPENMP
				size_t id = omp_get_thread_num();
				IndIterator it = pop.indIterator(sp->subPop(), id);
				idx = idx + id * (pop.subPopSize(sp->subPop()) / numThreads()) * (ploidy.end() - ploidy.begin()) * (loci.end() - loci.begin());
#else
				IndIterator it = pop.indIterator(sp->subPop());
#endif
				for (; it.valid(); ++it)
					for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
						for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc, ++idx)
							it->setAllele(ToAllele(m_genotype[idx % sz]), *loc, static_cast<int>(*p));

			}
#ifdef _OPENMP
			idx = idx + pop.subPopSize(sp->subPop()) * (ploidy.end() - ploidy.begin()) * (loci.end() - loci.begin());
#endif
		} else if (!m_prop.empty()) {
			WeightedSampler ws;
			size_t sz = pop.subPopSize(*sp);
			if (m_haplotypes.empty()) {
				// initialize by allele. Gurantee proportion at each locus.
				for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc) {
					ws.set(m_prop.begin(), m_prop.end(), sz * ploidy.size());
					IndIterator it = pop.indIterator(sp->subPop());
					for (; it.valid(); ++it)
						for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
							it->setAllele(ToAllele(ws.draw()), *loc, static_cast<int>(*p));
				}
			} else {
				ws.set(m_prop.begin(), m_prop.end(), sz * ploidy.size());
#pragma omp parallel if(numThreads() > 1)
				{
#ifdef _OPENMP
					IndIterator it = pop.indIterator(sp->subPop(), omp_get_thread_num());
#else
					IndIterator it = pop.indIterator(sp->subPop());
#endif
					for (; it.valid(); ++it)
						for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p) {
							const vectori & haplotype = m_haplotypes[ws.draw()];
							size_t hapSz = haplotype.size();
							size_t j = 0;
							for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc, ++j)
								it->setAllele(ToAllele(haplotype[j % hapSz]), *loc, static_cast<int>(*p));
						}
				}
			}
		} else {
			// m_freq can be empty if ....
			WeightedSampler ws(m_freq.empty() ? vectorf(m_haplotypes.size(), 1. / m_haplotypes.size()) : m_freq);
#pragma omp parallel if(numThreads() > 1)
			{
#ifdef _OPENMP
				IndIterator it = pop.indIterator(sp->subPop(), omp_get_thread_num());
#else
				IndIterator it = pop.indIterator(sp->subPop());
#endif
				for (; it.valid(); ++it)
					for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p) {
						if (m_haplotypes.empty()) {
							for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc)
								it->setAllele(ToAllele(ws.draw()), *loc, static_cast<int>(*p));
						} else {
							const vectori & haplotype = m_haplotypes[ws.draw()];
							size_t hapSz = haplotype.size();
							size_t j = 0;
							for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc, ++j)
								it->setAllele(ToAllele(haplotype[j % hapSz]), *loc, static_cast<int>(*p));
						}
					}
			}
		}
		pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


}
