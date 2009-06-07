/**
 *  $File: mutator.cpp $
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

void mutator::initialize(population & pop)
{

	DBG_DO(DBG_MUTATOR, cout << "initialize mutator" << endl);

	// deal with applicable loci
	if (m_loci.empty() ) {
		// all loci
		m_loci.resize(pop.totNumLoci() );
		for (UINT i = 0, iEnd = pop.totNumLoci(); i < iEnd;  ++i)
			m_loci[i] = i;
	}
#ifndef OPTIMIZED
	else {
		for (UINT i = 0; i < m_loci.size(); ++i) {
			if (m_loci[i] >= pop.totNumLoci())
				throw IndexError("Given loci is out of range");
		}
	}
#endif

	// all use the same rate
	if (m_rates.size() < m_loci.size() ) {
		m_rates.resize(m_loci.size());
		fill(m_rates.begin() + 1, m_rates.end(), m_rates[0]);
	}

	m_bt.setParameter(m_rates, pop.ploidy() * pop.popSize());

#ifndef OPTIMIZED
	for (size_t i = 0; i < m_rates.size(); ++i)
		if (fcmp_lt(m_rates[i], 0.) || fcmp_gt(m_rates[i], 1.) )
			throw ValueError("Migration rate should be between [0,1], given " + toStr(m_rates[i]));
#endif
	if (pop.totNumLoci() != m_mutCount.size())
		m_mutCount.resize(pop.totNumLoci(), 0);
	m_initialized = true;
}


bool mutator::apply(population & pop)
{
	if (!m_initialized || m_bt.trialSize() != pop.ploidy() * pop.popSize()) {
		initialize(pop);
		DBG_DO(DBG_MUTATOR, cout << "Reinitialize mutator at loci" << m_loci <<
			" at rate " << m_rates << endl);
	}

	DBG_DO(DBG_MUTATOR, cout << "Mutate replicate " << pop.rep() << endl);

	m_bt.doTrial();

	// mapIn and mapOut
	bool mapIn = !m_mapIn.empty() || m_mapIn.func().isValid();
	vectorlu const & mapInList = m_mapIn.elems();
	pyFunc mapInFunc = m_mapIn.func();
	UINT numMapInAllele = mapInList.size();
	bool mapOut = !m_mapOut.empty() || m_mapOut.func().isValid();
	vectorlu const & mapOutList = m_mapOut.elems();
	UINT numMapOutAllele = mapOutList.size();
	pyFunc mapOutFunc = m_mapOut.func();
	// mutate each mutable locus
	for (size_t i = 0, iEnd = m_loci.size(); i < iEnd; ++i) {
		int locus = m_loci[i];
		DBG_DO(DBG_MUTATOR, cout << "Mutate at locus " << locus << endl);
		size_t pos = m_bt.trialFirstSucc(i);
		if (pos != BernulliTrials::npos) {
			do {
				IndAlleleIterator ptr = pop.alleleIterator(locus) + pos;
				if (!ptr.valid())
					continue;
				DBG_DO(DBG_MUTATOR, cout << "Allele " << int(*ptr) << " at locus " << locus);
				if (mapIn) {
					if (numMapInAllele > 0) {
						if (static_cast<size_t>(*ptr) < numMapInAllele)
							*ptr = ToAllele(mapInList[*ptr]);
					} else {
						*ptr = ToAllele(mapInFunc(PyObj_As_Int, "(i)",
								static_cast<int>(*ptr)));
					}
				}
				// The virtual mutate functions in derived operators will be called.
				mutate(*ptr);
				if (mapOut) {
					if (numMapOutAllele > 0) {
						if (static_cast<size_t>(*ptr) < numMapOutAllele)
							*ptr = ToAllele(mapOutList[*ptr]);
					} else {
						*ptr = ToAllele(mapOutFunc(PyObj_As_Int, "(i)",
								static_cast<int>(*ptr)));
					}

				}
				DBG_DO(DBG_MUTATOR, cout << " is mutated to " << int(*ptr) << endl);
				m_mutCount[ locus ]++;
			} while ( (pos = m_bt.trialNextSucc(i, pos)) != BernulliTrials::npos);
		}                                                                                           // succ.any
	}                                                                                               // each applicable loci

	return true;
}


matrixMutator::matrixMutator(const matrix & rate,
	const uintList & loci, const uintListFunc & mapIn, const uintListFunc & mapOut,
	const stringFunc & output,
	int stage, int begin, int end, int step, const intList & at,
	const repList & rep, const subPopList & subPops,
	const stringList & infoFields)
	: mutator(vectorf(1, 0), loci, mapIn, mapOut, output, stage, begin, end, step,
	          at, rep, subPops, infoFields)
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
	DBG_DO(DBG_MUTATOR, cout << "Mu " << mu << endl);
	setRate(vectorf(1, mu), loci.elems());
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
		DBG_DO(DBG_MUTATOR, cout << "Setting weight for allele " << i << " to " << rateMatrix[i] << endl);
		m_sampler.push_back(weightedSampler(rng(), rateMatrix[i]));
	}
}


void matrixMutator::mutate(AlleleRef allele)
{
	DBG_FAILIF(allele >= m_sampler.size(), IndexError,
		"Allele out of range of 1 ~ " + toStr(m_sampler.size() - 1)
		+ " (determined by the size of the mutation rate matrix).");
	allele = ToAllele(m_sampler[allele].get());
}


// mutate to a state other than current state with equal probability
void kamMutator::mutate(AlleleRef allele)
{
#ifdef BINARYALLELE
	allele = !allele;
#else
	Allele new_allele = static_cast<Allele>(rng().randInt(m_k - 1));
	if (new_allele >= allele)
		allele = new_allele + 1;
	else
		allele = new_allele;
#endif
}


smmMutator::smmMutator(const floatList & rates, const uintList & loci,
	double incProb, UINT maxAllele, const floatListFunc & mutStep,
	const uintListFunc & mapIn, const uintListFunc & mapOut, const stringFunc & output,
	int stage, int begin, int end, int step, const intList & at,
	const repList & rep, const subPopList & subPops, const stringList & infoFields)
	: mutator(rates, loci, mapIn, mapOut, output, stage, begin, end, step, at, rep, subPops, infoFields),
	m_incProb(incProb), m_maxAllele(maxAllele), m_mutStep(mutStep)
{
#ifdef BINARYALLELE
	DBG_WARNING(true, "Symetric stepwise mutation does not work well on two state alleles.");
#endif
	DBG_ASSERT(fcmp_ge(m_incProb, 0.) && fcmp_le(m_incProb, 1.),
		ValueError, "Inc probability should be between [0,1], given " + toStr(m_incProb));
	
	if (m_maxAllele == 0)
		m_maxAllele = MaxAllele();
	if (m_maxAllele > MaxAllele())
		throw ValueError("maxAllele exceeds maximum allowed allele in this module.");
	
	DBG_FAILIF(! m_mutStep.func().isValid() && m_mutStep.empty(), ValueError,
		"Parameter mutStep must be a number, a list or a valid function.");
	
	DBG_FAILIF(m_mutStep.size() > 1 && 
		(fcmp_lt(m_mutStep[1], 0) || fcmp_gt(m_mutStep[1], 1)), ValueError,
		"Probability for the geometric distribution has to be between 0 and 1");
	
}


void smmMutator::mutate(AlleleRef allele)
{
	UINT step = 1;
	if (m_mutStep.size() == 1)
		step = static_cast<UINT>(m_mutStep[0]);
	else if (m_mutStep.size() == 2) {
		DBG_ASSERT(static_cast<int>(m_mutStep[0]) == GeometricDistribution, ValueError,
			"Incorrect mode for generating mutation step.");
		step = rng().randGeometric(m_mutStep[1]);
	} else {
		DBG_ASSERT(m_mutStep.func().isValid(), ValueError,
			"Invalid Python function for smmMutator");
		step = m_mutStep.func()(PyObj_As_Int, "(i)", static_cast<int>(allele));
	}

	// increase
	if (rng().randUniform01() < m_incProb) {
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


// mutate according to the mixed model
void pyMutator::mutate(AlleleRef allele)
{
	int resInt = m_func(PyObj_As_Int, "(i)", static_cast<int>(allele));

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


bool pointMutator::apply(population & pop)
{
	m_mutCount.resize(pop.totNumLoci(), 0);
	// mutate each mutable locus
	for (size_t i = 0, iEnd = m_loci.size(); i < iEnd; ++i) {
		for (vectorlu::iterator ind = m_inds.begin();
		     ind != m_inds.end(); ++ind) {
			for (size_t p = 0; p < m_ploidy.size(); ++p) {
				m_mutCount[m_loci[i]]++;
				*(pop.ind(*ind).genoBegin(m_ploidy[p]) + m_loci[i]) = m_allele;
				DBG_DO(DBG_MUTATOR, cout << "Mutate locus " << m_loci[i] <<
					" to allele " << toStr(m_allele) << " at generation " << pop.gen() << endl);
			}
		}
	}                                                                                 // each applicable loci

	return true;
}


}
