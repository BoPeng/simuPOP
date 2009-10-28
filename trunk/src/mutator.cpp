/**
 *  $File: baseMutator.cpp $
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

#include "mutator.h"

namespace simuPOP {

double baseMutator::mutRate(UINT loc)
{
	DBG_FAILIF(m_rates.empty(), ValueError, "Please specify mutation rate.");
	if (m_loci.allAvail()) {
		DBG_ASSERT(m_rates.size() == 1 || m_rates.size() > loc,
			IndexError, "Locus index out of range when retrieving mutation rate.");
		return m_rates.size() == 1 ? m_rates[0] : m_rates[loc];
	}
	const vectoru & loci = m_loci.elems();
	vectoru::const_iterator it = find(loci.begin(), loci.end(), loc);

	DBG_ASSERT(it != loci.end(), RuntimeError,
		"Failed to find mutation rate for locus " + toStr(loc));
	DBG_ASSERT(m_rates.size() == loci.size(), SystemError,
		"Incorrect rate and loci size");
	return m_rates[it - loci.begin()];
}


void baseMutator::fillContext(const population & pop, IndAlleleIterator ptr, UINT locus)
{
	// chromosome?
	UINT chrom = pop.chromLocusPair(locus).first;
	UINT beg = pop.chromBegin(chrom);
	UINT end = pop.chromEnd(chrom);
	UINT cnt = m_context.size() / 2;

	for (size_t i = 0; i < cnt; ++i) {
		if (locus >= beg + i)
			m_context[i] = *(ptr.ptr() - (cnt - i));
		else
			m_context[i] = -1;
	}
	for (size_t i = 0; i < cnt; ++i) {
		if (locus + i < end)
			m_context[cnt + i] = *(ptr.ptr() + i + 1);
		else
			m_context[cnt + i] = -1;
	}
	if (!m_mapIn.empty() || m_mapIn.func().isValid()) {
		for (size_t i = 0; i < m_context.size(); ++i) {
			if (m_context[i] == -1)
				continue;
			vectoru const & mapInList = m_mapIn.elems();
			if (mapInList.size() > 0) {
				if (static_cast<UINT>(m_context[i]) < mapInList.size())
					m_context[i] = mapInList[m_context[i]];
			} else {
				m_context[i] = m_mapIn.func() (PyObj_As_Int, "(i)",
					m_context[i]);
			}
		}
	}
}


bool baseMutator::apply(population & pop)
{
	DBG_DO(DBG_MUTATOR, cerr << "Mutate replicate " << pop.rep() << endl);

	// mapIn and mapOut
	bool mapIn = !m_mapIn.empty() || m_mapIn.func().isValid();
	vectoru const & mapInList = m_mapIn.elems();
	pyFunc mapInFunc = m_mapIn.func();
	UINT numMapInAllele = mapInList.size();
	bool mapOut = !m_mapOut.empty() || m_mapOut.func().isValid();
	vectoru const & mapOutList = m_mapOut.elems();
	UINT numMapOutAllele = mapOutList.size();
	pyFunc mapOutFunc = m_mapOut.func();
	// mutate each mutable locus

	subPopList subPops = applicableSubPops();
	// the usual whole population, easy case.
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);

	BernulliTrials bt(GetRNG());

	DBG_FAILIF(m_rates.empty(), ValueError, "Please specify mutation rate or rates.");
	// all use the same rate
	vectorf rates = m_rates;
	if (rates.size() == 1) {
		rates.resize(m_loci.allAvail() ? pop.totNumLoci() : m_loci.elems().size());
		fill(rates.begin() + 1, rates.end(), rates[0]);
	}
	// multiple (virtual) subpopulations
	for (UINT idx = 0; idx < subPops.size(); ++idx) {
		UINT sp = subPops[idx].subPop();

		// fromSubPops out of range....
		DBG_FAILIF(sp >= pop.numSubPop(), IndexError,
			"Subpopulation index " + toStr(sp) + " out of range");

		ULONG popSize = pop.subPopSize(subPops[idx]);
		DBG_DO(DBG_MUTATOR, cerr << "SP " << subPops[idx] << " size: " << popSize << endl);
		if (popSize == 0)
			continue;

		if (subPops[idx].isVirtual())
			pop.activateVirtualSubPop(subPops[idx]);

		bt.setParameter(rates, pop.ploidy() * popSize);

		bt.doTrial();
		const vectoru & loci = m_loci.elems();
		size_t iEnd = m_loci.allAvail() ? pop.totNumLoci() : loci.size();
		for (size_t i = 0; i < iEnd; ++i) {
			UINT locus = m_loci.allAvail() ? i : loci[i];
			DBG_DO(DBG_MUTATOR, cerr << "Mutate at locus " << locus << endl);
			size_t pos = bt.trialFirstSucc(i);
			size_t lastPos = 0;
			IndAlleleIterator ptr = pop.alleleIterator(locus, sp);
			if (pos != BernulliTrials::npos) {
				do {
					ptr += pos - lastPos;
					lastPos = pos;
					if (!ptr.valid())
						continue;
					DBG_DO(DBG_MUTATOR, cerr << "Allele " << int(*ptr) << " at locus " << locus);
					if (mapIn) {
						if (numMapInAllele > 0) {
							if (static_cast<size_t>(*ptr) < numMapInAllele)
								*ptr = ToAllele(mapInList[*ptr]);
						} else {
							*ptr = ToAllele(mapInFunc(PyObj_As_Int, "(i)",
									static_cast<int>(*ptr)));
						}
					}
					if (!m_context.empty())
						fillContext(pop, ptr, locus);
					// The virtual mutate functions in derived operators will be called.
					mutate(*ptr, locus);
					if (mapOut) {
						if (numMapOutAllele > 0) {
							if (static_cast<size_t>(*ptr) < numMapOutAllele)
								*ptr = ToAllele(mapOutList[*ptr]);
						} else {
							*ptr = ToAllele(mapOutFunc(PyObj_As_Int, "(i)",
									static_cast<int>(*ptr)));
						}
					}
					DBG_DO(DBG_MUTATOR, cerr << " is mutated to " << int(*ptr) << endl);
				} while ( (pos = bt.trialNextSucc(i, pos)) != BernulliTrials::npos);
			}                                                                                           // succ.any
		}

		if (subPops[idx].isVirtual())
			pop.deactivateVirtualSubPop(sp);
	}   // each subpopulation
	return true;
}


matrixMutator::matrixMutator(const matrix & rate,
	const uintList & loci, const uintListFunc & mapIn, const uintListFunc & mapOut,
	const stringFunc & output,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops,
	const stringList & infoFields)
	: baseMutator(vectorf(1, 0), loci, mapIn, mapOut, 0, output, begin, end, step,
	          at, reps, subPops, infoFields)
{
	matrix rateMatrix = rate;
	// step 0, determine mu
	double mu = 0;

	for (size_t i = 0; i < rateMatrix.size(); ++i) {
		DBG_ASSERT(rateMatrix[i].size() == rateMatrix.size(), ValueError,
			"A n by n matrix is required.");
		double sum = 0;
		for (size_t j = 0; j < rateMatrix[i].size(); ++j) {
			// ignore p_ii
			if (i == j)
				continue;
			DBG_FAILIF(rateMatrix[i][j] < 0 || rateMatrix[i][j] > 1, ValueError,
				"Elements in a mutation matrix must be between 0 and 1. " + toStr(rateMatrix[i][j]) + " observed.");
			sum += rateMatrix[i][j];
		}
		DBG_FAILIF(sum > 1, ValueError, "Sum of P_ij should not exceed 1");
		if (mu < sum)
			mu = sum;
	}
	DBG_DO(DBG_MUTATOR, cerr << "Mu " << mu << endl);
	setRate(vectorf(1, mu), loci);
	if (mu == 0.)
		return;
	// re-calculate probability
	m_sampler.clear();
	for (size_t i = 0; i < rateMatrix.size(); ++i) {
		double sum = 0;
		for (size_t j = 0; j < rateMatrix[i].size(); ++j) {
			if (i == j)
				continue;
			sum += rateMatrix[i][j];
			rateMatrix[i][j] /= mu;
		}
		rateMatrix[i][i] = 1 - sum / mu;
		DBG_DO(DBG_MUTATOR, cerr << "Setting weight for allele " << i << " to " << rateMatrix[i] << endl);
		m_sampler.push_back(weightedSampler(GetRNG(), rateMatrix[i]));
	}
}


void matrixMutator::mutate(AlleleRef allele, UINT)
{
	DBG_FAILIF(allele >= m_sampler.size(), IndexError,
		"Allele out of range of 1 ~ " + toStr(m_sampler.size() - 1)
		+ " (determined by the size of the mutation rate matrix).");
	allele = ToAllele(m_sampler[allele].get());
}


// mutate to a state other than current state with equal probability
void kamMutator::mutate(AlleleRef allele, UINT)
{
#ifdef BINARYALLELE
	allele = !allele;
#else
	Allele new_allele = static_cast<Allele>(GetRNG().randInt(m_k - 1));
	if (new_allele >= allele)
		allele = new_allele + 1;
	else
		allele = new_allele;
#endif
}


smmMutator::smmMutator(const floatList & rates, const uintList & loci,
	double incProb, UINT maxAllele, const floatListFunc & mutStep,
	const uintListFunc & mapIn, const uintListFunc & mapOut, const stringFunc & output,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops, const stringList & infoFields)
	: baseMutator(rates, loci, mapIn, mapOut, 0, output, begin, end, step, at, reps, subPops, infoFields),
	m_incProb(incProb), m_maxAllele(maxAllele), m_mutStep(mutStep)
{
#ifdef BINARYALLELE
	DBG_WARNING(true, "Symetric stepwise mutation does not work well on two state alleles.");
#endif
	DBG_ASSERT(fcmp_ge(m_incProb, 0.) && fcmp_le(m_incProb, 1.),
		ValueError, "Inc probability should be between [0,1], given " + toStr(m_incProb));

	if (m_maxAllele == 0)
		m_maxAllele = ModuleMaxAllele;
	if (static_cast<ULONG>(m_maxAllele) > ModuleMaxAllele)
		throw ValueError("maxAllele exceeds maximum allowed allele in this module.");

	DBG_FAILIF(!m_mutStep.func().isValid() && m_mutStep.empty(), ValueError,
		"Parameter mutStep must be a number, a list or a valid function.");

	DBG_FAILIF(m_mutStep.size() > 1 &&
		(fcmp_lt(m_mutStep[1], 0) || fcmp_gt(m_mutStep[1], 1)), ValueError,
		"Probability for the geometric distribution has to be between 0 and 1");

}


void smmMutator::mutate(AlleleRef allele, UINT)
{
	UINT step = 1;

	if (m_mutStep.size() == 1)
		step = static_cast<UINT>(m_mutStep[0]);
	else if (m_mutStep.size() == 2) {
		DBG_ASSERT(static_cast<int>(m_mutStep[0]) == GeometricDistribution, ValueError,
			"Incorrect mode for generating mutation step.");
		step = GetRNG().randGeometric(m_mutStep[1]);
	} else {
		DBG_ASSERT(m_mutStep.func().isValid(), ValueError,
			"Invalid Python function for smmMutator");
		step = m_mutStep.func() (PyObj_As_Int, "(i)", static_cast<int>(allele));
	}

	// increase
	if (GetRNG().randUniform() < m_incProb) {
#ifdef BINARYALLELE
		allele = 1;
#else
		if (static_cast<UINT>(allele + step) < m_maxAllele)
			AlleleAdd(allele, step);
		else
			allele = ToAllele(m_maxAllele);
#endif
	}
	// decrease
	else {
#ifdef BINARYALLELE
		allele = 0;
#else
		if (allele > step)
			AlleleMinus(allele, step);
		else
			allele = 0;
#endif
	}
}


void pyMutator::mutate(AlleleRef allele, UINT)
{
	int resInt = 0;

	if (m_contextObj == NULL && !context().empty())
		// this needs to be done only once
		m_contextObj = Int_Vec_As_NumArray(context().begin(), context().end());

	if (m_contextObj != NULL)
		resInt = m_func(PyObj_As_Int, "(iO)", static_cast<int>(allele), m_contextObj);
	else
		resInt = m_func(PyObj_As_Int, "(i)", static_cast<int>(allele));

#ifdef BINARYALLELE
	DBG_ASSERT(resInt == 0 || resInt == 1, ValueError,
		"Can only mutate to 0 or 1 in binary mode.");
	allele = resInt != 0;
#else
	DBG_ASSERT(static_cast<unsigned>(resInt) <= ModuleMaxAllele, ValueError,
		"Mutated to an allele greater than maximum allowed allele value");
	allele = static_cast<Allele>(resInt);
#endif
}


void mixedMutator::initialize(population & pop)
{
	baseMutator::initialize(pop);
	for (size_t i = 0; i < m_mutators.size(); ++i)
		reinterpret_cast<baseMutator *>(m_mutators[i])->initialize(pop);
}


void mixedMutator::mutate(AlleleRef allele, UINT locus)
{
	UINT idx = m_sampler.get();
	baseMutator * mut = reinterpret_cast<baseMutator *>(m_mutators[idx]);
	double mu = mut->mutRate(locus);

	if (mu == 1.0 || GetRNG().randUniform() < mu)
		mut->mutate(allele, locus);
}


void contextMutator::initialize(population & pop)
{
	baseMutator::initialize(pop);
	for (size_t i = 0; i < m_mutators.size(); ++i)
		reinterpret_cast<baseMutator *>(m_mutators[i])->initialize(pop);
}


void contextMutator::mutate(AlleleRef allele, UINT locus)
{
	const vectori & alleles = context();

	for (size_t i = 0; i < m_contexts.size(); ++i) {
		bool match = true;
		for (size_t j = 0; j < alleles.size(); ++j) {
			if (m_contexts[i][j] != alleles[j]) {
				match = false;
				break;
			}
		}
		if (match) {
			DBG_DO(DBG_MUTATOR, cerr << "Context " << alleles << " mutator " << i << endl);
			baseMutator * mut = reinterpret_cast<baseMutator *>(m_mutators[i]);
			if (GetRNG().randUniform() < mut->mutRate(locus))
				mut->mutate(allele, locus);
			return;
		}
	}
	if (m_contexts.size() + 1 == m_mutators.size()) {
		DBG_DO(DBG_MUTATOR, cerr << "No context found. Use last mutator." << endl);
		baseMutator * mut = reinterpret_cast<baseMutator *>(m_mutators[m_contexts.size()]);
		if (GetRNG().randUniform() < mut->mutRate(locus))
			mut->mutate(allele, locus);
	} else {
		cerr << "Failed to find context " << alleles << endl;
		throw RuntimeError("No match context is found and there is no default mutator");
	}
}


bool pointMutator::apply(population & pop)
{
	subPopList subPops = applicableSubPops();

	// the usual whole population, easy case.
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);

	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for (; sp != spEnd; ++sp) {
		if (sp->isVirtual())
			pop.activateVirtualSubPop(*sp);

		vectoru::iterator it = m_inds.begin();
		vectoru::iterator itEnd = m_inds.end();
		for (; it != itEnd; ++it) {
			IndIterator ind = pop.indIterator(sp->subPop(), VisibleInds) + *it;
			if (!ind.valid())
				continue;
			for (size_t p = 0; p < m_ploidy.size(); ++p) {
				if (m_loci.allAvail()) {
					for (size_t i = 0; i < pop.totNumLoci(); ++i) {
						ind->setAllele(m_allele, i, m_ploidy[p]);
						DBG_DO(DBG_MUTATOR,
							cerr << "Mutate locus " << i << " at ploidy " << m_ploidy[p]
							     << " to allele " << int(m_allele) << " at generation "
							     << pop.gen() << endl);
					}
				} else {
					const vectoru & loci = m_loci.elems();
					for (size_t i = 0; i < loci.size(); ++i) {
						ind->setAllele(m_allele, loci[i], m_ploidy[p]);
						DBG_DO(DBG_MUTATOR,
							cerr << "Mutate locus " << loci[i] << " at ploidy " << m_ploidy[p]
							     << " to allele " << int(m_allele) << " at generation "
							     << pop.gen() << endl);
					}
				}
			}   // ploidy
		}       // individual
		if (sp->isVirtual())
			pop.deactivateVirtualSubPop(sp->subPop());
	}       // subpopulation

	return true;
}


}
