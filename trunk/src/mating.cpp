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

#include "mating.h"

namespace simuPOP {

offspringGenerator::offspringGenerator(const vectorop & ops,
	const floatListFunc & numOffspring, const floatList & sexMode) :
	m_numOffspring(numOffspring), m_sexMode(sexMode),
	m_transmitters(0), m_formOffGenotype(true), m_initialized(false)
{
	DBG_FAILIF(numOffspring.size() == 0 && !numOffspring.func().isValid(), ValueError,
		"Please specify either number of offspring or a function.");

	if (numOffspring.size() == 1) {
		DBG_FAILIF(static_cast<int>(numOffspring[0]) <= 0, ValueError,
			"Number of offspring has to be positive.");
	} else if (m_numOffspring.size() > 1) {
		int mode = static_cast<int>(m_numOffspring[0]);
		(void)mode;  // fix compiler warning.

		DBG_FAILIF(mode == BinomialDistribution
			&& (m_numOffspring.size() < 3 || static_cast<UINT>(m_numOffspring[2]) <= 1.),
			ValueError, "If mode is BinomialDistribution, the second parameter should be greater than 1");
		DBG_FAILIF(mode == UniformDistribution &&
			(m_numOffspring.size() < 3 || m_numOffspring[2] < static_cast<UINT>(m_numOffspring[1])),
			ValueError, "If mode is UniformDistribution, m_numOffspringParam should be greater than m_numOffspring");
		DBG_FAILIF(mode == GeometricDistribution &&
			(m_numOffspring.size() < 2 || fcmp_lt(m_numOffspring[1], 0) || fcmp_gt(m_numOffspring[1], 1.)),
			ValueError, "P for a geometric distribution should be within [0,1]");
		DBG_FAILIF(mode == BinomialDistribution &&
			(m_numOffspring.size() < 2 || fcmp_lt(m_numOffspring[1], 0) || fcmp_gt(m_numOffspring[1], 1.)),
			ValueError, "P for a Bionomial distribution should be within [0,1].");
		DBG_FAILIF(mode == BinomialDistribution &&
			(m_numOffspring.size() < 3 || m_numOffspring[2] < 1),
			ValueError, "Max number of offspring should be greater than 1.");
	}
	DBG_FAILIF(m_sexMode.empty(), ValueError, "Please specify one of the sex modes");
	DBG_FAILIF((static_cast<int>(m_sexMode[0]) == ProbOfMale ||
		        static_cast<int>(m_sexMode[0]) == NumOfMale ||
		        static_cast<int>(m_sexMode[0]) == NumOfFemale) && m_sexMode.size() < 2,
		ValueError, "A parameter is required for sex mode ProbOfMale, NumOfMale and NumOfFemale");

	DBG_FAILIF(static_cast<int>(m_sexMode[0]) == ProbOfMale &&
		(fcmp_lt(m_sexMode[1], 0) || fcmp_gt(m_sexMode[1], 1)),
		ValueError, "Probability of male has to be between 0 and 1");
	// the genotype transmitter that will be used when no during mating
	// operator is used to transmit genotype from parents to offspring.
	for (size_t i = 0; i < ops.size(); ++i)
		m_transmitters.push_back(ops[i]->clone());
}


offspringGenerator::offspringGenerator(const offspringGenerator & rhs)
	: m_numOffspring(rhs.m_numOffspring),
	m_sexMode(rhs.m_sexMode), m_transmitters(0),
	m_formOffGenotype(rhs.m_formOffGenotype),
	m_initialized(rhs.m_initialized)
{
	for (size_t i = 0; i < rhs.m_transmitters.size(); ++i)
		m_transmitters.push_back(rhs.m_transmitters[i]->clone());
}


ULONG offspringGenerator::numOffspring(int gen)
{
	if (m_numOffspring.size() == 1)
		return static_cast<UINT>(m_numOffspring[0]);

	if (m_numOffspring.func().isValid()) {
		int numOff = m_numOffspring.func() (PyObj_As_Int, "(i)", gen);
		DBG_FAILIF(numOff < 1, ValueError, "Need at least one offspring.");
		return numOff;
	}
	switch (static_cast<int>(m_numOffspring[0])) {
	case GeometricDistribution:
		return rng().randGeometric(m_numOffspring[1]);
	case PoissonDistribution:
		return rng().randPoisson(m_numOffspring[1]) + 1;
	case BinomialDistribution:
		return rng().randBinomial(static_cast<UINT>(m_numOffspring[2]) - 1, m_numOffspring[1]) + 1;
	case UniformDistribution:
		// max: 5
		// num: 2
		// randint(4)  ==> 0, 1, 2, 3
		// + 2 ==> 2, 3, 4, 5
		return rng().randInt(static_cast<UINT>(m_numOffspring[2]) - static_cast<UINT>(m_numOffspring[1]) + 1)
		       + static_cast<UINT>(m_numOffspring[1]);
	default:
		throw ValueError("Wrong mating numoffspring mode. Should be one of \n"
			             "NumOffspring, NumOffspringEachFamily and GEometricDistribution");
	}
	//
	DBG_ASSERT(false, SystemError, "This line should never be reached.");
	return 0;
}


Sex offspringGenerator::getSex(int count)
{
	int mode = static_cast<int>(m_sexMode[0]);

	if (mode == NoSex)
		return Male;
	else if (mode == RandomSex)
		return rng().randBit() ? Male : Female;
	else if (mode == ProbOfMale)
		return rng().randUniform01() < m_sexMode[1] ? Male : Female;
	else if (mode == NumOfMale)
		return count < static_cast<int>(m_sexMode[1]) ? Male : Female;
	else if (mode == NumOfFemale)
		return count < static_cast<int>(m_sexMode[1]) ? Female : Male;
	DBG_ASSERT(false, SystemError, "This line should not be reached.");
	return Male;
}


void offspringGenerator::initialize(const population & pop, SubPopID subPop, vector<baseOperator *> const & ops)
{
	m_formOffGenotype = true;
	vector<baseOperator *>::const_iterator iop = ops.begin();
	vector<baseOperator *>::const_iterator iop_end = ops.end();
	for (; iop != ops.end(); ++iop) {
		if ((*iop)->formOffGenotype())
			m_formOffGenotype = false;
	}
	m_initialized = true;
}


UINT offspringGenerator::generateOffspring(population & pop, individual * dad, individual * mom,
                                           RawIndIterator & it,
                                           RawIndIterator & itEnd,
                                           vector<baseOperator *> & ops)
{
	DBG_ASSERT(initialized(), ValueError,
		"Offspring generator is not initialized before used to generate offspring");

	// generate numOff offspring per mating, or until it  reaches offEnd
	UINT count = 0;
	bool accept = true;
	UINT numOff = numOffspring(pop.gen());
	UINT attempt = 0;
	while (count < numOff && attempt < numOff && it != itEnd) {
		++attempt;

		DBG_FAILIF(m_transmitters.empty() && ops.empty(), ValueError,
			"No valid genotype transmitter is defined.");

		// set sex, during mating operator will try to
		// follow the offspring sex (e.g. pass X or Y chromosome)
		it->setSex(getSex(count));
		//
		accept = true;
		vector<baseOperator *>::iterator iop = m_transmitters.begin();
		vector<baseOperator *>::iterator iopEnd = m_transmitters.end();
		for (; iop != iopEnd; ++iop) {
			try {
				// During mating operator might reject this offspring.
				if (!m_formOffGenotype && (*iop)->formOffGenotype())
					continue;

				if (!(*iop)->applyDuringMating(pop, it, dad, mom)) {
					accept = false;
					break;
				}
			} catch (...) {
				cout << "One of the transmitters " << (*iop)->__repr__()
				     << " throws an exception." << endl << endl;
				throw;
			}
		}

		if (!accept)
			continue;

		// apply during mating operators
		iop = ops.begin();
		iopEnd = ops.end();
		for (; iop != iopEnd; ++iop) {
			try {
				// During mating operator might reject this offspring.
				if (!(*iop)->applyDuringMating(pop, it, dad, mom)) {
					accept = false;
					break;
				}
			} catch (...) {
				cout << "DuringMating operator " << (*iop)->__repr__()
				     << " throws an exception." << endl << endl;
				throw;
			}
		}                                                                         // all during-mating operators
		if (accept) {
			it++;
			count++;
		}
	}
	return count;
}


controlledOffspringGenerator::controlledOffspringGenerator(
	const vectori & loci, const vectori & alleles, PyObject * freqFunc,
	const vectorop & ops, const floatListFunc & numOffspring,
	const floatList & sexMode)
	: offspringGenerator(ops, numOffspring, sexMode),
	m_loci(loci), m_alleles(alleles), m_freqFunc(freqFunc),
	m_expAlleles(), m_totAllele(), m_curAllele()
{
	if (!m_freqFunc.isValid())
		throw ValueError("Please specify a valid frequency function");
}


controlledOffspringGenerator::controlledOffspringGenerator(const controlledOffspringGenerator & rhs)
	: offspringGenerator(rhs),
	m_loci(rhs.m_loci),
	m_alleles(rhs.m_alleles),
	m_freqFunc(rhs.m_freqFunc)
{
}


void controlledOffspringGenerator::getExpectedAlleles(const population & pop,
                                                      vectorf & expFreq)
{
	// determine expected number of alleles of each allele
	// at each subpopulation.
	UINT nLoci = m_loci.size();
	UINT numSP = pop.numSubPop();

	DBG_ASSERT(expFreq.size() == nLoci || expFreq.size() == nLoci * numSP, SystemError,
		"Expect expected frequency for all loci (or all loci in all subpopulation)");

	for (size_t i = 0; i < expFreq.size(); ++i) {
		DBG_FAILIF(expFreq[i] < 0 || expFreq[i] > 1., ValueError,
			"Expected frequency out of range");
	}

	// in the order of Loc1: sp1, sp2, sp3, sloc2: p1, sp2, sp3
	m_expAlleles.resize(nLoci * numSP);
	if (numSP > 1 && expFreq.size() == nLoci) {
		// if exp frequencies in subpopulation is not specified.
		// I need to find the current allele frequencies
		// and use them as proportions for the next generation.
		for (size_t i = 0; i < nLoci; ++i) {
			int locus = m_loci[i];
			Allele allele = static_cast<Allele>(m_alleles[i]);

			// determine the number alleles at each subpopulation.
			vectorf curFreq(numSP, 0);
			bool hasAllele = false;
			for (size_t sp = 0, n = 0; sp < numSP; ++sp) {
				IndAlleleIterator a = const_cast<population &>(pop).alleleBegin(locus, sp);
				IndAlleleIterator aEnd = const_cast<population &>(pop).alleleEnd(locus, sp);
				for (; a != aEnd; ++a)
					if (AlleleUnsigned(*a) == allele)
						n++;
				hasAllele = hasAllele || n > 0;
				curFreq[sp] = double(n) / (pop.subPopSize(sp) * pop.ploidy());
			}

			DBG_DO(DBG_MATING, cout << "Current frequency at locus " << locus
				                    << " is " << curFreq << endl);

			// if there is no alleles
			if (!hasAllele && expFreq[i] > 0.)
				throw RuntimeError("No disease allele exists at generation "
					+ toStr(pop.gen()) + ", but exp allele frequency is greater than 0.");

			// calculate exp number of affected offspring in the next generation.
			//
			// step 1: totalsize*expFreq is the total number of disease alleles
			// step 2: assign these alleles to each subpopulation according to a multi-nomial
			// distribution with p_i beging allele frequency at each subpopulation.
			// assign these numbers to each subpopulation
			rng().randMultinomial(static_cast<unsigned int>(pop.popSize() * expFreq[i] * pop.ploidy()),
				curFreq, m_expAlleles.begin() + numSP * i);
		}
	} else {
		// simpler case, one subpopulation, or with gieven allele frequency
		for (size_t i = 0; i < nLoci; ++i) {
			for (size_t sp = 0; sp < numSP; ++sp) {
#ifndef OPTIMIZED
				int locus = m_loci[i];
				Allele allele = m_alleles[i];
				ULONG n = 0;
				// go through all alleles
				IndAlleleIterator a = const_cast<population &>(pop).alleleBegin(locus, sp);
				IndAlleleIterator aEnd = const_cast<population &>(pop).alleleEnd(locus, sp);
				for (; a != aEnd; ++a) {
					if (AlleleUnsigned(*a) == allele)
						n++;
				}

				// if there is no alleles
				if (n == 0 && expFreq[numSP * i + sp] > 0.)
					throw RuntimeError("No disease allele exists at generation " +
						toStr(pop.gen()) + " but exp allele frequency is greater than 0.");
#endif
				m_expAlleles[numSP * i + sp] = static_cast<UINT>(pop.subPopSize(sp) * pop.ploidy() * expFreq[numSP * i + sp]);
				if (expFreq[numSP * i + sp] > 0. && m_expAlleles[numSP * i + sp] == 0)
					m_expAlleles[numSP * i + sp] = 1;
			}
		}
	}
}


void controlledOffspringGenerator::initialize(const population & pop, SubPopID subPop, vector<baseOperator *> const & ops)
{
	offspringGenerator::initialize(pop, subPop, ops);

	// expected frequency at each locus
	if (subPop == 0) {
		vectorf expFreq = m_freqFunc(PyObj_As_Array, "(i)", pop.gen());
		DBG_DO(DBG_MATING, cout << "expected freq " << expFreq << endl);

		//
		// determine expected number of alleles of each allele
		// at each subpopulation.
		getExpectedAlleles(pop, expFreq);
		DBG_DO(DBG_MATING, cout << "expected alleles " << m_expAlleles << endl);
	}

	// now, for **this** subpopulation...

	// total allowed disease alleles.
	m_totAllele.resize(m_loci.size());
	fill(m_totAllele.begin(), m_totAllele.end(), 0);
	// currently available disease allele. (in the offspring generation)
	m_curAllele.resize(m_loci.size());
	fill(m_curAllele.begin(), m_curAllele.end(), 0);
	//
	m_flip.resize(m_loci.size());
	fill(m_flip.begin(), m_flip.end(), false);

	// We control allele 1 if expected allele frequency is less than 0.5
	for (size_t i = 0; i < m_loci.size(); ++i) {
		UINT maxCount = pop.subPopSize(subPop) * pop.ploidy();
		m_totAllele[i] = m_expAlleles[subPop + pop.numSubPop() * i];
		if (m_totAllele[i] > maxCount) {
			cout << "Warning: number of planned affected alleles exceed population size.";
			m_totAllele[i] = maxCount;
		}
		if (2 * m_totAllele[i] > maxCount) {
			m_flip[i] = true;
			m_totAllele[i] = maxCount - m_totAllele[i];
		}
	}
	//

	// it is possible that no disease allele is required
	// so everyone is accepted
	m_freqRequMet = true;
	//
	for (size_t i = 0; i < m_loci.size(); ++i) {
		if (m_totAllele[i] > 0) {
			m_freqRequMet = false;
			break;
		}
	}
	// if mor ethan noAAattempt times, no pure homo found,
	// accept non-homo cases.
	m_AAattempt = 200;
	m_aaAttempt = 200;
}


/// CPPONLY
UINT controlledOffspringGenerator::generateOffspring(population & pop, individual * dad, individual * mom,
                                                     RawIndIterator & offBegin,
                                                     RawIndIterator & offEnd,
                                                     vector<baseOperator *> & ops)
{
	UINT nLoci = m_loci.size();
	//
	// generate m_numOffspring offspring per mating
	// record family size (this may be wrong for the last family)
	//
	RawIndIterator itBegin = offBegin;
	UINT numOff = offspringGenerator::generateOffspring(pop,
		dad, mom, offBegin, offEnd, ops);

	//
	if (numOff == 0)
		return numOff;
	//
	// count number of alleles in the family.
	vectori na(nLoci, 0);
	bool hasAff = false;
	UINT totNumLoci = pop.totNumLoci();
	// we know that scratch population has ordered linear genotype
	for (size_t i = 0; i < nLoci; ++i) {
		GenoIterator ptr = itBegin->genoBegin() + m_loci[i];
		for (size_t j = 0; j < numOff * pop.ploidy(); ++j, ptr += totNumLoci) {
			if (m_flip[i] ? (*ptr != static_cast<Allele>(m_alleles[i]))
				: (*ptr == static_cast<Allele>(m_alleles[i]))) {
				na[i]++;
				hasAff = true;
			}
		}
	}
	// now check if this family is usable.
	bool accept = false;
	// all disease alleles have been satidfied.
	// only accept unaffected families
	// otherwise, accept any family that can help.
	if (m_freqRequMet) {
		if (!hasAff) {
			// has AA, so no need to compromise
			m_AAattempt = 10000;
			accept = true;
		}
		// tried 200 times, no AA is found.
		else if (m_AAattempt == 0) {
			m_AAattempt = 200;
			accept = true;
		}
		m_AAattempt--;
	} else {                                                          // do not use stack
		if (hasAff) {
			for (size_t i = 0; i < nLoci; ++i) {
				// accept the whole family, if we need this allele
				if (m_curAllele[i] < m_totAllele[i] && na[i] > 0) {
					accept = true;
					break;
				}
			}
			m_aaAttempt = 10000;
		} else if (m_aaAttempt == 0) {
			m_aaAttempt = 200;
			accept = true;
		}
		m_aaAttempt--;
	}
	//
	// reject this family
	if (!accept) {
		// it relocate to its begin point
		// DBG_DO(DBG_MATING, cout << "Reject " << na << endl);
		offBegin = itBegin;
		return 0;
	}
	// accpet this family, see if all done.
	for (size_t i = 0; i < nLoci; ++i)
		m_curAllele[i] += na[i];
	if (!m_freqRequMet) {
		m_freqRequMet = true;
		for (size_t i = 0; i < nLoci; ++i) {
			if (m_curAllele[i] < m_totAllele[i])
				m_freqRequMet = false;
		}
	}
	return numOff;
}


void sequentialParentChooser::initialize(population & pop, SubPopID sp)
{
	m_begin = pop.indBegin(sp);
	m_end = pop.indEnd(sp);
	m_ind = m_begin;
	m_initialized = true;
}


parentChooser::individualPair sequentialParentChooser::chooseParents(RawIndIterator)
{
	if (m_ind == m_end)
		m_ind = m_begin;
	return parentChooser::individualPair(& * m_ind++, NULL);
}


void sequentialParentsChooser::initialize(population & pop, SubPopID subPop)
{
	m_numMale = 0;
	m_numFemale = 0;
	m_curMale = 0;
	m_curFemale = 0;
	m_maleIndex.clear();
	m_femaleIndex.clear();

	IndIterator it = pop.indBegin(subPop);
	for (; it.valid(); ++it) {
		if (it->sex() == Male) {
			m_numMale++;
			m_maleIndex.push_back(it.rawIter());
		} else {
			m_numFemale++;
			m_femaleIndex.push_back(it.rawIter());
		}
	}
	m_initialized = true;
}


parentChooser::individualPair sequentialParentsChooser::chooseParents(RawIndIterator)
{
	DBG_ASSERT(initialized(), SystemError,
		"Please initialize this parent chooser before using it");

	individual * dad = NULL;
	individual * mom = NULL;

	if (m_curMale == m_numMale)
		m_curMale = 0;
	if (m_curFemale == m_numFemale)
		m_curFemale = 0;

	// using weighted sampler.
	if (m_numMale != 0)
		dad = & * (m_maleIndex[m_curMale++]);
	else
		dad = & * (m_femaleIndex[m_curFemale++]);

	if (m_numFemale != 0)
		mom = & * (m_femaleIndex[m_curFemale++]);
	else
		mom = & * (m_maleIndex[m_curMale++]);
	return std::make_pair(dad, mom);
}


void randomParentChooser::initialize(population & pop, SubPopID sp)
{
	m_selection = pop.selectionOn(sp);
	m_index.clear();

	// In a virtual subpopulation, because m_begin + ... is **really** slow
	// It is a good idea to cache IndIterators. This is however inefficient
	// for non-virtual populations
	if (pop.hasActivatedVirtualSubPop(sp) && !m_replacement) {
		IndIterator it = pop.indBegin(sp);
		for (; it.valid(); ++it)
			m_index.push_back(it.rawIter());
	}

	if (m_selection) {
		UINT fit_id = pop.infoIdx("fitness");
		// regardless of sex, get fitness for everyone.
		m_sampler.set(vectorf(pop.infoBegin(fit_id, sp),
				pop.infoEnd(fit_id, sp)));
	} else {
		m_size = m_index.size();
		if (m_size == 0)     // if m_index is not used (no VSP)
			m_size = pop.subPopSize(sp);
	}

	if (!m_replacement)
		std::random_shuffle(m_index.begin(), m_index.end());

	DBG_FAILIF(!m_replacement && m_selection, ValueError,
		"Selection is not allowed in random sample without replacement");

	m_initialized = true;
}


parentChooser::individualPair randomParentChooser::chooseParents(RawIndIterator basePtr)
{
	DBG_ASSERT(initialized(), SystemError,
		"Please initialize this parent chooser before using it");
	individual * ind = NULL;
	// choose a parent
	if (!m_replacement) {
		DBG_FAILIF(m_index.empty(), ValueError,
			"All parents have been chosen.");
		ind = & * m_index.back();
		m_index.pop_back();
	}
	if (m_index.empty()) {
		if (m_selection)
			ind = & * (basePtr + m_sampler.get());
		else
			ind = & * (basePtr + rng().randInt(m_size));
	} else {
		if (m_selection)
			ind = & * (m_index[m_sampler.get()]);
		else
			ind = & * (m_index[rng().randInt(m_size)]);
	}
	return individualPair(ind, NULL);
}


void randomParentsChooser::initialize(population & pop, SubPopID subPop)
{
	m_numMale = 0;
	m_numFemale = 0;

	IndIterator it = pop.indBegin(subPop);
	for (; it.valid(); ++it) {
		if (it->sex() == Male)
			m_numMale++;
		else
			m_numFemale++;
	}

	// allocate memory at first for performance reasons
	m_maleIndex.resize(m_numMale);
	m_femaleIndex.resize(m_numFemale);

	m_selection = pop.selectionOn(subPop);
	UINT fit_id = 0;
	if (m_selection) {
		fit_id = pop.infoIdx("fitness");
		m_maleFitness.resize(m_numMale);
		m_femaleFitness.resize(m_numFemale);
	}

	m_numMale = 0;
	m_numFemale = 0;

	it = pop.indBegin(subPop);
	for (; it.valid(); it++) {
		if (it->sex() == Male) {
			m_maleIndex[m_numMale] = it.rawIter();
			if (m_selection)
				m_maleFitness[m_numMale] = it->info(fit_id);
			m_numMale++;
		} else {
			m_femaleIndex[m_numFemale] = it.rawIter();
			if (m_selection)
				m_femaleFitness[m_numFemale] = it->info(fit_id);
			m_numFemale++;
		}
	}

	if (!m_replacement) {
		DBG_FAILIF(m_maleIndex.empty(), IndexError, "No male individual in this population");
		DBG_FAILIF(m_femaleIndex.empty(), IndexError, "No female individual in this population");
		std::random_shuffle(m_maleIndex.begin(), m_maleIndex.end());
		std::random_shuffle(m_femaleIndex.begin(), m_femaleIndex.end());
	}

	if (m_selection) {
		m_malesampler.set(m_maleFitness);
		m_femalesampler.set(m_femaleFitness);
		DBG_DO(DBG_DEVEL, cout << "Male fitness " << m_maleFitness << endl);
		DBG_DO(DBG_DEVEL, cout << "Female fitness " << m_femaleFitness << endl);
	}

	DBG_FAILIF(!m_replacement && m_selection, ValueError,
		"Selection is not allowed in random sample without replacement");

	m_initialized = true;
}


parentChooser::individualPair randomParentsChooser::chooseParents(RawIndIterator)
{
	DBG_ASSERT(initialized(), SystemError,
		"Please initialize this parent chooser before using it");

	individual * dad = NULL;
	individual * mom = NULL;

	if (!m_replacement) {
		DBG_FAILIF(m_femaleIndex.empty(), ValueError,
			"All females have been chosen.");
		mom = & * m_femaleIndex.back();
		m_femaleIndex.pop_back();

		DBG_FAILIF(m_maleIndex.empty(), ValueError,
			"All males have been chosen.");
		dad = & * m_maleIndex.back();
		m_maleIndex.pop_back();
		return std::make_pair(dad, mom);
	}

	// this exception should be raised also in optimized mode because the cause
	// can be random.
	if (m_numMale == 0)
		throw RuntimeError("randomParentsChooser fails because there is no male individual in a subpopulation.");
	if (m_numFemale == 0)
		throw RuntimeError("randomParentsChooser fails because there is no female individual in a subpopulation ");

	if (m_selection) {
		// using weighted sampler.
		dad = & * (m_maleIndex[m_malesampler.get()]);
		mom = & * (m_femaleIndex[m_femalesampler.get()]);
	} else {
		dad = & * (m_maleIndex[rng().randInt(m_numMale)]);
		mom = & * (m_femaleIndex[rng().randInt(m_numFemale)]);
	}
	return std::make_pair(dad, mom);
}


void polyParentsChooser::initialize(population & pop, SubPopID subPop)
{
	m_numMale = 0;
	m_numFemale = 0;
	m_polyCount = 0;

	IndIterator it = pop.indBegin(subPop);
	for (; it.valid(); ++it) {
		if (it->sex() == Male)
			m_numMale++;
		else
			m_numFemale++;
	}

	// allocate memory at first for performance reasons
	m_maleIndex.resize(m_numMale);
	m_femaleIndex.resize(m_numFemale);

	m_selection = pop.selectionOn(subPop);
	UINT fit_id = 0;
	if (m_selection) {
		fit_id = pop.infoIdx("fitness");
		m_maleFitness.resize(m_numMale);
		m_femaleFitness.resize(m_numFemale);
	}

	m_numMale = 0;
	m_numFemale = 0;

	it = pop.indBegin(subPop);
	for (; it.valid(); it++) {
		if (it->sex() == Male) {
			m_maleIndex[m_numMale] = it.rawIter();
			if (m_selection)
				m_maleFitness[m_numMale] = it->info(fit_id);
			m_numMale++;
		} else {
			m_femaleIndex[m_numFemale] = it.rawIter();
			if (m_selection)
				m_femaleFitness[m_numFemale] = it->info(fit_id);
			m_numFemale++;
		}
	}

	if (m_selection) {
		m_malesampler.set(m_maleFitness);
		m_femalesampler.set(m_femaleFitness);
		DBG_DO(DBG_DEVEL, cout << "Male fitness " << m_maleFitness << endl);
		DBG_DO(DBG_DEVEL, cout << "Female fitness " << m_femaleFitness << endl);
	}

	m_initialized = true;
}


parentChooser::individualPair polyParentsChooser::chooseParents(RawIndIterator)
{
	DBG_ASSERT(initialized(), SystemError,
		"Please initialize this parent chooser before using it");

	individual * dad = NULL;
	individual * mom = NULL;

	if (m_polyNum > 1 && m_polyCount > 0) {
		if (m_polySex == Male)
			dad = m_lastParent;
		else
			mom = m_lastParent;
		m_polyCount--;
	}

	// using weidhted sampler.
	if (dad == NULL) {
		if (m_numMale == 0)
			throw RuntimeError("polyParentsChooser fails because there is no male individual in a subpopulation.");

		if (m_selection)
			dad = & * (m_maleIndex[m_malesampler.get()]);
		else
			dad = & * (m_maleIndex[rng().randInt(m_numMale)]);

		if (m_polySex == Male && m_polyNum > 1) {
			m_polyCount = m_polyNum - 1;
			m_lastParent = dad;
		}
	}

	if (mom == NULL) {
		if (m_numFemale == 0)
			throw RuntimeError("polyParentsChooser fails because there is no female individual in a subpopulation.");

		if (m_selection)
			mom = & * (m_femaleIndex[m_femalesampler.get()]);
		else
			mom = & * (m_femaleIndex[rng().randInt(m_numFemale)]);

		if (m_polySex == Female && m_polyNum > 1) {
			m_polyCount = m_polyNum - 1;
			m_lastParent = mom;
		}
	}
	return std::make_pair(dad, mom);
}


void alphaParentsChooser::initialize(population & pop, SubPopID subPop)
{
	m_numMale = 0;
	m_numFemale = 0;

	UINT info_id = 0;
	bool hasAlphaMale = m_alphaNum != 0 || !m_alphaField.empty();
	bool useInfo = false;
	if (!m_alphaField.empty()) {
		info_id = pop.infoIdx(m_alphaField);
		useInfo = true;
	}

	IndIterator it = pop.indBegin(subPop);
	for (; it.valid(); ++it) {
		if (hasAlphaMale && useInfo && it->sex() == m_alphaSex
		    && it->info(info_id) == 0.)
			continue;
		if (it->sex() == Male)
			m_numMale++;
		else
			m_numFemale++;
	}

	DBG_FAILIF(hasAlphaMale && (m_numMale == 0 || m_numFemale == 0),
		ValueError, "No alpha individual or individual of opposite sex is found.");

	// allocate memory at first for performance reasons
	m_maleIndex.resize(m_numMale);
	m_femaleIndex.resize(m_numFemale);

	m_selection = pop.selectionOn(subPop);
	UINT fit_id = 0;
	if (m_selection) {
		fit_id = pop.infoIdx("fitness");
		m_maleFitness.resize(m_numMale);
		m_femaleFitness.resize(m_numFemale);
	}

	m_numMale = 0;
	m_numFemale = 0;

	it = pop.indBegin(subPop);
	for (; it.valid(); it++) {
		if (hasAlphaMale && useInfo && it->sex() == m_alphaSex
		    && it->info(info_id) == 0.)
			continue;
		if (it->sex() == Male) {
			m_maleIndex[m_numMale] = it.rawIter();
			if (m_selection)
				m_maleFitness[m_numMale] = it->info(fit_id);
			m_numMale++;
		} else {
			m_femaleIndex[m_numFemale] = it.rawIter();
			if (m_selection)
				m_femaleFitness[m_numFemale] = it->info(fit_id);
			m_numFemale++;
		}
	}

	if (m_selection) {
		m_malesampler.set(m_maleFitness);
		m_femalesampler.set(m_femaleFitness);
		DBG_DO(DBG_DEVEL, cout << "Male fitness " << m_maleFitness << endl);
		DBG_DO(DBG_DEVEL, cout << "Female fitness " << m_femaleFitness << endl);
	}

	if (!hasAlphaMale || useInfo || m_alphaNum >=
	    (m_alphaSex == Male ? m_numMale : m_numFemale)) {
		m_initialized = true;
		return;
	}

	// now, we need to choose a few alpha individuals
	vector<RawIndIterator> m_newAlphaIndex;
	vectorf m_newAlphaFitness;
	// select individuals
	for (size_t i = 0; i < m_alphaNum; ++i) {
		if (m_selection) {     // fix me, without replacement!
			// using weighted sampler.
			m_newAlphaIndex.push_back(
				m_alphaSex == Male ? m_maleIndex[m_malesampler.get()]
				: m_femaleIndex[m_femalesampler.get()]);
			m_newAlphaFitness.push_back(m_newAlphaIndex.back()->info(fit_id));
		} else     // fix me, without replacement!
			m_newAlphaIndex.push_back(
				m_alphaSex == Male ? m_maleIndex[rng().randInt(m_numMale)]
				: m_femaleIndex[rng().randInt(m_numFemale)]);
	}
	if (m_alphaSex == Male) {
		m_maleIndex.swap(m_newAlphaIndex);
		if (m_selection)
			m_malesampler.set(m_newAlphaFitness);
		m_numMale = m_maleIndex.size();
	} else {
		m_femaleIndex.swap(m_newAlphaIndex);
		if (m_selection)
			m_femalesampler.set(m_newAlphaFitness);
		m_numFemale = m_femaleIndex.size();
	}
	m_initialized = true;
}


parentChooser::individualPair alphaParentsChooser::chooseParents(RawIndIterator)
{
	DBG_ASSERT(initialized(), SystemError,
		"Please initialize this parent chooser before using it");

	individual * dad = NULL;
	individual * mom = NULL;

	// this exception should be raised also in optimized mode because the cause
	// can be random.
	if (m_numMale == 0)
		throw RuntimeError("alphaParentsChooser fails because there is no male individual in a subpopulation.");
	if (m_numFemale == 0)
		throw RuntimeError("alphaParentsChooser fails because there is no female individual in a subpopulation ");

	// using weidhted sampler.
	if (m_selection) {                                    // with selection
		dad = & * (m_maleIndex[m_malesampler.get()]);
		mom = & * (m_femaleIndex[m_femalesampler.get()]);
	} else {
		dad = & * (m_maleIndex[rng().randInt(m_numMale)]);
		mom = & * (m_femaleIndex[rng().randInt(m_numFemale)]);
	}

	return std::make_pair(dad, mom);
}


void infoParentsChooser::initialize(population & pop, SubPopID sp)
{
	if (m_func.isValid()) {
		PyObject * popObj = pyPopObj(static_cast<void *>(&pop));
		// if pop is valid?
		if (popObj == NULL)
			throw SystemError("Could not pass population to the provided function. \n"
				              "Compiled with the wrong version of SWIG?");

		// parammeter list, ref count increased
		bool resBool;
		if (m_param.isValid())
			resBool = m_func(PyObj_As_Bool, "(OO)", popObj, m_param.object());
		else
			resBool = m_func(PyObj_As_Bool, "(O)", popObj);

		Py_DECREF(popObj);
	}

	// indexes
	m_infoIdx.resize(m_infoFields.size());
	for (size_t i = 0; i < m_infoFields.size(); ++i)
		m_infoIdx[i] = pop.infoIdx(m_infoFields[i]);
	UINT infoSz = m_infoIdx.size();

	m_selection = pop.selectionOn(sp);
	m_index.clear();

	// In a virtual subpopulation, because m_begin + ... is **really** slow
	// It is a good idea to cache IndIterators. This is however inefficient
	// for non-virtual populations
	IndIterator it = pop.indBegin(sp);
	vectorf fitness;
	UINT fit_id = 0;
	if (m_selection)
		fit_id = pop.infoIdx("fitness");
	for (; it.valid(); ++it) {
		Sex mySex = it->sex();
		for (size_t i = 0; i < infoSz; ++i)
			// we only choose individual with an valid information field
			// and is of opposite sex
			if (it->info(m_infoIdx[i]) >= 0 && pop.ind(it->intInfo(m_infoIdx[i])).sex() != mySex) {
				m_index.push_back(it.rawIter());
				if (m_selection)
					fitness.push_back(it->info(fit_id));
				break;
			}
	}
	//
	m_degenerate = m_index.empty();
	DBG_WARNING(m_degenerate, "Parents are chosen randomly because there is no valid index.");
	if (m_degenerate) {
		for (it = pop.indBegin(sp); it.valid(); ++it) {
			m_index.push_back(it.rawIter());
			if (m_selection)
				fitness.push_back(it->info(fit_id));
		}
	}

	if (m_selection)
		m_sampler.set(fitness);
	else
		m_size = m_index.size();

	if (!m_replacement)
		std::random_shuffle(m_index.begin(), m_index.end());

	DBG_FAILIF(!m_replacement && m_selection, ValueError,
		"Selection is not allowed in random sample without replacement");

	m_initialized = true;
}


parentChooser::individualPair infoParentsChooser::chooseParents(RawIndIterator basePtr)
{
	DBG_ASSERT(initialized(), SystemError,
		"Please initialize this parent chooser before using it");
	individual * par1 = randomParentChooser::chooseParents(basePtr).first;
	Sex sex1 = par1->sex();
	// there is no valid information field value
	if (m_degenerate) {
		int attempt = 0;
		while (++attempt < 1000) {
			individual * par2 = randomParentChooser::chooseParents(basePtr).first;
			if (par2->sex() != sex1)
				return sex1 == Male ? std::make_pair(par1, par2) : std::make_pair(par2, par1);
		}
		throw RuntimeError("Can not locate any individual of opposite sex");
	}
	// the way this parent chooser is initialized guranttees that
	// theres is at lest one valid field.
	vector<individual *> validInds;
	for (size_t i = 0; i < m_infoIdx.size(); ++i) {
		int info = par1->intInfo(m_infoIdx[i]);
		if (info < 0)
			continue;
		RawIndIterator par2 = basePtr + info;
		if (par2->sex() != sex1)
			validInds.push_back(& * par2);
	}
	DBG_FAILIF(validInds.empty(), SystemError, "No valid relative is found");
	individual * par2 = validInds[rng().randInt(validInds.size())];
	DBG_DO(DBG_DEVEL, cout << "infoParentsChooser: par1: " << par1 - & * basePtr
		                   << " par2: " << par2 - & * basePtr << endl);
	return sex1 == Male ? std::make_pair(par1, par2) : std::make_pair(par2, par1);
}


pyParentsChooser::pyParentsChooser(PyObject * pc)
	: parentChooser(), m_func(pc), m_generator(NULL), m_parIterator(NULL)
{
}


void pyParentsChooser::initialize(population & pop, SubPopID sp)
{
#if PY_VERSION_HEX < 0x02040000
	throw SystemError("Your Python version does not have good support for generator"
		              " so this python parent chooser can not be used.");
#else
#  ifndef OPTIMIZED
	m_size = pop.subPopSize(sp);
#  endif
	m_begin = pop.indBegin(sp);

	PyObject * popObj = pyPopObj(static_cast<void *>(&pop));
	// if pop is valid?
	DBG_FAILIF(popObj == NULL, SystemError,
		"Could not pass population to the provided function. \n"
		"Compiled with the wrong version of SWIG?");

	m_generator = m_func("(Oi)", popObj, sp);

	// test if m_generator is a generator
	DBG_ASSERT(PyGen_Check(m_generator), ValueError,
		"Passed function is not a python generator");

	m_parIterator = PyObject_GetIter(m_generator);

	// test if m_parIterator is iteratable.
	DBG_FAILIF(m_parIterator == NULL, ValueError,
		"Can not iterate through parent generator");
#endif
	m_initialized = true;
}


parentChooser::individualPair pyParentsChooser::chooseParents(RawIndIterator)
{
	DBG_ASSERT(initialized(), SystemError,
		"Please initialize this parent chooser before using it");

	PyObject * item = PyIter_Next(m_parIterator);

#ifndef OPTIMIZED
	if (item == NULL && debug(DBG_GENERAL)) {
		PyErr_Print();
		PyErr_Clear();
	}
	DBG_FAILIF(item == NULL, ValueError,
		"User-defined function yield invalid value.");
#endif

	vectori parents;
	int parent;

	if (PySequence_Check(item)) {
		PyObj_As_IntArray(item, parents);
		DBG_ASSERT(parents.size() == 2, ValueError,
			"Returned parents indexes should have size 2");
#ifndef OPTIMIZED
		DBG_ASSERT(static_cast<unsigned>(parents[0]) < m_size
			&& static_cast<unsigned>(parents[1]) < m_size,
			ValueError, "Returned parent index (" + toStr(parents[0])
			+ ", and " + toStr(parents[1]) +
			") is greater than subpopulation size " + toStr(m_size));
#endif
		Py_DECREF(item);
		// FIXME: this can be really slow in a virtual population because
		// visibility between m_begin and m_begin+parents[x] need to
		// be checked. It should make sense to return individual directly.
		return std::make_pair(& * (m_begin + parents[0]),
			& * (m_begin + parents[1]));
	} else if (PyInt_Check(item) || PyLong_Check(item)) {
		PyObj_As_Int(item, parent);
#ifndef OPTIMIZED
		DBG_ASSERT(static_cast<unsigned>(parent) < m_size,
			ValueError, "Returned index (" + toStr(parent) +
			") is greater than subpopulation size " + toStr(m_size));
#endif
		Py_DECREF(item);
		// FIXME: this can be really slow in a virtual population because
		// visibility between m_begin and m_begin+parent need to
		// be checked. It should make sense to return individual directly.
		return parentChooser::individualPair(& * (m_begin + parent), NULL);
	} else
		DBG_ASSERT(false, ValueError,
			"Invalid type of returned parent index(es)");
	// this should not be reached
	return parentChooser::individualPair(NULL, NULL);
}


mating::mating(const uintListFunc & subPopSize)
	: m_subPopSize(subPopSize)
{
}


bool mating::prepareScratchPop(population & pop, population & scratch)
{
	// use population structure of pop
	if (m_subPopSize.empty() && !m_subPopSize.func().isValid())
		scratch.fitSubPopStru(pop.subPopSizes(), pop.subPopNames());
	else if (!m_subPopSize.empty())                                                     // set subPoplation size
		scratch.fitSubPopStru(m_subPopSize.elems(), pop.subPopNames());
	else {                                                                              // use m_subPopSizeFunc
		// get generation number
		int gen = pop.gen();
		// convert current pop size to a tuple
		PyObject * curSize = PyTuple_New(pop.numSubPop());

		DBG_ASSERT(curSize != NULL, SystemError, "Can not convert current pop size to a list");

		for (size_t i = 0; i < pop.numSubPop(); ++i)
			PyTuple_SetItem(curSize, i, PyInt_FromLong(pop.subPopSize(i)));

		vectori res = m_subPopSize.func() (PyObj_As_IntArray, "(iO)", gen, curSize);
		Py_XDECREF(curSize);

		vectorlu sz(res.size());

		for (size_t i = 0; i < res.size(); i++)
			sz[i] = static_cast<ULONG>(res[i]);

		// allow change of pop size of scratch
		scratch.fitSubPopStru(sz, pop.subPopNames());
	}
	// this is not absolutely necessary but will reduce confusions
	scratch.setVirtualSplitter(pop.virtualSplitter());

	DBG_DO(DBG_SIMULATOR, cout << "New subpop size " << scratch.subPopSizes() << endl);

	DBG_FAILIF(scratch.numSubPop() != pop.numSubPop(),
		ValueError, "number of subPopulaitons must agree.\n Pre: "
		+ toStr(pop.numSubPop()) + " now: " + toStr(scratch.numSubPop() ));
	return true;
}


bool mating::mate(population & pop, population & scratch,
                  vector<baseOperator * > & ops)
{
	// scrtach will have the right structure.
	if (!prepareScratchPop(pop, scratch))
		return false;

	for (SubPopID sp = 0; sp < static_cast<SubPopID>(pop.numSubPop()); ++sp)
		if (!mateSubPop(pop, sp, scratch.rawIndBegin(sp),
				scratch.rawIndEnd(sp), ops))
			return false;
	submitScratch(pop, scratch);
	return true;
}


void mating::submitScratch(population & pop, population & scratch)
{
	pop.turnOffSelection();
	// use scratch population,
	pop.push(scratch);
	scratch.validate("after push and discard");
}


pedigreeMating::pedigreeMating(const pedigree & ped,
	const offspringGenerator & generator, bool setSex, bool setAffection,
	const vectorstr & copyFields)
	: mating(uintListFunc()), m_ped(ped),
	m_setSex(setSex), m_setAffection(setAffection), m_copyFields(copyFields)
{
	m_generator = generator.clone();
}


pedigreeMating::pedigreeMating(const pedigreeMating & rhs)
	: mating(rhs), m_ped(rhs.m_ped), m_setSex(rhs.m_setSex),
	m_setAffection(rhs.m_setAffection), m_copyFields(rhs.m_copyFields)
{
	m_generator = rhs.m_generator->clone();
	DBG_FAILIF(m_ped.ancestralGens() == 0, ValueError,
		"Passed pedigree has no ancestral generation.");
	// scroll to the greatest generation, but this generation
	// should have no parental generation.
	m_ped.useAncestralGen(m_ped.ancestralGens());
}


pedigreeMating::~pedigreeMating()
{
	delete m_generator;
}


bool pedigreeMating::prepareScratchPop(population & pop, population & scratch)
{
	DBG_FAILIF(pop.numSubPop() != m_ped.numSubPop(), ValueError,
		"Evolving generation does not have the same number of subpopulation as the pedigree.");
	for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
		DBG_WARNING(pop.subPopSize(sp) > m_ped.subPopSize(sp),
			"Giving population has more individuals than the pedigree."
			"Some of the parents will be ignored");
		DBG_FAILIF(pop.subPopSize(sp) < m_ped.subPopSize(sp), ValueError,
			"Given population has less individuals in subpopulation " + toStr(sp)
			+ " than the pedigree. PedigreeMating cannot continue.");
	}
	if (m_ped.curAncestralGen() == 0)
		return false;

	// copy information to the greatest ancestral generation
	if (m_ped.curAncestralGen() == m_ped.ancestralGens() &&
	    (m_setSex || m_setAffection || !m_copyFields.empty())) {
		vectoru infoIdx;
		vectoru pedInfoIdx;
		for (size_t i = 0; i < m_copyFields.size(); ++i) {
			infoIdx.push_back(pop.infoIdx(m_copyFields[i]));
			pedInfoIdx.push_back(m_ped.infoIdx(m_copyFields[i]));
		}
		for (size_t it = 0; it < pop.popSize(); ++it) {
			individual & ind = pop.ind(it);
			individual & pedInd = m_ped.ind(it);
			if (m_setSex)
				ind.setSex(pedInd.sex());
			if (m_setAffection)
				ind.setAffected(pedInd.affected());
			for (size_t i = 0; i < m_copyFields.size(); ++i)
				ind.setInfo(pedInd.info(pedInfoIdx[i]), infoIdx[i]);
		}
	}
	m_parentalPopSize = m_ped.popSize();
	m_ped.useAncestralGen(m_ped.curAncestralGen() - 1);
	scratch.fitSubPopStru(m_ped.subPopSizes(), m_ped.subPopNames());
	return true;
}


bool pedigreeMating::mate(population & pop, population & scratch,
                          vector<baseOperator * > & ops)
{
	// scrtach will have the right structure.
	if (!prepareScratchPop(pop, scratch))
		return false;
	vectoru infoIdx;
	vectoru pedInfoIdx;
	for (size_t i = 0; i < m_copyFields.size(); ++i) {
		infoIdx.push_back(pop.infoIdx(m_copyFields[i]));
		pedInfoIdx.push_back(m_ped.infoIdx(m_copyFields[i]));
	}

	for (SubPopID sp = 0; sp < static_cast<SubPopID>(scratch.numSubPop()); ++sp) {
		if (!m_generator->initialized())
			m_generator->initialize(pop, sp, ops);

		RawIndIterator it = scratch.rawIndBegin(sp);
		RawIndIterator itEnd;
		for (size_t i = 0; i < scratch.subPopSize(sp); ++i) {
			int father_idx = m_ped.father(i, sp);
			DBG_FAILIF(father_idx > m_parentalPopSize, IndexError,
				"Parental index " + toStr(father_idx) + " out of range of 0 - "
				+ toStr(m_parentalPopSize - 1));
			individual * dad = father_idx >= 0 ? &pop.ind(father_idx) : NULL;

			int mother_idx = m_ped.mother(i, sp);
			DBG_FAILIF(mother_idx > m_parentalPopSize, IndexError,
				"Parental index " + toStr(mother_idx) + " out of range of 0 - "
				+ toStr(m_parentalPopSize - 1));
			individual * mom = mother_idx >= 0 ? &pop.ind(mother_idx) : NULL;

			if (m_setSex)
				it->setSex(m_ped.ind(i, sp).sex());
			if (m_setAffection)
				it->setAffected(m_ped.ind(i, sp).affected());
			for (size_t i = 0; i < m_copyFields.size(); ++i)
				it->setInfo(m_ped.ind(i, sp).info(pedInfoIdx[i]), infoIdx[i]);

			//
			itEnd = it + 1;
			// whatever the numOffspring function returns for this
			// offspring generator, only generate one offspring.
			UINT numOff = m_generator->generateOffspring(pop, dad, mom,
				it, itEnd, ops);
			(void)numOff;
			DBG_FAILIF(numOff != 1, RuntimeError,
				"Generation of offspring must succeed in pedigreeMating");
		}
		m_generator->finalize(pop);
	}

	submitScratch(pop, scratch);
	return true;
}


homoMating::homoMating(parentChooser & chooser,
	offspringGenerator & generator,
	const uintListFunc & subPopSize,
	vspID subPop, double weight)
	: mating(subPopSize), m_subPop(subPop), m_weight(weight)
{
	m_parentChooser = chooser.clone();
	m_offspringGenerator = generator.clone();
}


bool homoMating::mateSubPop(population & pop, SubPopID subPop,
                            RawIndIterator offBegin, RawIndIterator offEnd,
                            vector<baseOperator * > & ops)
{
	// nothing to do.
	if (offBegin == offEnd)
		return true;

	if (!m_parentChooser->initialized())
		m_parentChooser->initialize(pop, subPop);

	if (!m_offspringGenerator->initialized())
		m_offspringGenerator->initialize(pop, subPop, ops);

	// generate scratch.subPopSize(sp) individuals.
	RawIndIterator it = offBegin;
	while (it != offEnd) {
		individual * dad = NULL;
		individual * mom = NULL;
		parentChooser::individualPair const parents = m_parentChooser->chooseParents(pop.rawIndBegin());
		dad = parents.first;
		mom = parents.second;

		//
		UINT numOff = m_offspringGenerator->generateOffspring(pop, dad, mom, it, offEnd, ops);
		(void)numOff;             // silent warning about unused variable.
	}
	m_parentChooser->finalize(pop, subPop);
	m_offspringGenerator->finalize(pop);
	return true;
}


heteroMating::heteroMating(const vectormating & matingSchemes,
	const uintListFunc & subPopSize,
	bool shuffleOffspring)
	: mating(subPopSize),
	m_shuffleOffspring(shuffleOffspring)
{
	vectormating::const_iterator it = matingSchemes.begin();
	vectormating::const_iterator it_end = matingSchemes.end();

	for (; it != it_end; ++it)
		m_matingSchemes.push_back(dynamic_cast<homoMating *>((*it)->clone()));
}


heteroMating::~heteroMating()
{
	vectormating::iterator it = m_matingSchemes.begin();
	vectormating::iterator it_end = m_matingSchemes.end();

	for (; it != it_end; ++it)
		delete *it;
}


heteroMating::heteroMating(const heteroMating & rhs) :
	mating(rhs), m_shuffleOffspring(rhs.m_shuffleOffspring)
{
	vectormating::const_iterator it = rhs.m_matingSchemes.begin();
	vectormating::const_iterator it_end = rhs.m_matingSchemes.end();

	for (; it != it_end; ++it)
		m_matingSchemes.push_back(dynamic_cast<homoMating *>((*it)->clone()));
}


bool heteroMating::mate(population & pop, population & scratch,
                        vector<baseOperator * > & ops)
{
	// scrtach will have the right structure.
	if (!prepareScratchPop(pop, scratch))
		return false;
	vectormating::const_iterator it = m_matingSchemes.begin();
	vectormating::const_iterator it_end = m_matingSchemes.end();

	for (SubPopID sp = 0; sp < static_cast<SubPopID>(pop.numSubPop()); ++sp) {
		vectormating m;
		vectorf w_pos;                          // positive weights
		vectorf w_neg;                          // negative weights
		vectormating::iterator it = m_matingSchemes.begin();
		vectormating::iterator it_end = m_matingSchemes.end();
		for (; it != it_end; ++it) {
			DBG_FAILIF((*it)->subPop() == InvalidSubPopID, ValueError,
				"Please specify which subpopulation each homogeneous mating scheme is applied to");
			// if it is used for this subpop,
			// for use for all subPops ...
			if ((*it)->subPop() == sp) {
				m.push_back(*it);
				double w = (*it)->weight();
				// less than zero...
				if (fcmp_lt(w, 0.)) {
					w_pos.push_back(0);
					w_neg.push_back(-w);
				} else {
					w_pos.push_back(w);
					w_neg.push_back(0);
				}
			}
		}
		DBG_FAILIF(m.empty(), ValueError,
			"No mating scheme is available for subpopulation " + toStr(sp));
		// determine the weight
		if (m.size() == 1) {
			w_pos[0] = 1.;
			w_neg[0] = 0.;
		}

		// the default case (all zero)
		if (fcmp_eq(std::accumulate(w_pos.begin(), w_pos.end(), 0.), 0.)) {
			// weight is subpopulation size
			for (size_t i = 0; i < m.size(); ++i)
				// if there is no negative weight...
				if (w_neg[i] == 0)
					w_pos[i] = pop.subPopSize(vspID(sp, m[i]->virtualSubPop()));
		}
		DBG_DO(DBG_DEVEL, cout << "Positive mating scheme weights: " << w_pos << '\n'
			                   << "Negative mating scheme weights: " << w_neg << endl);

		// weight.
		double overall_pos = std::accumulate(w_pos.begin(), w_pos.end(), 0.);
		double overall_neg = std::accumulate(w_neg.begin(), w_neg.end(), 0.);
		(void)overall_neg;                         // silent warning about unused variable.
		DBG_FAILIF(fcmp_eq(overall_pos, 0.) && fcmp_eq(overall_neg, 0.), ValueError,
			"Overall weight is zero");
		//
		vectorlu vspSize(m.size());
		ULONG all = scratch.subPopSize(sp);
		// first count negative ones
		for (size_t i = 0; i < m.size(); ++i) {
			if (fcmp_gt(w_neg[i], 0.)) {
				vspSize[i] = static_cast<ULONG>(pop.subPopSize(vspID(sp, m[i]->virtualSubPop())) * w_neg[i]);
				DBG_ASSERT(all >= vspSize[i], ValueError,
					"Not enough offspring to accommodate specified weight scheme. "
					"Current item is " + toStr(i) + " requires " + toStr(vspSize[i])
					+ " available " + toStr(all) + ".");
				all -= vspSize[i];
			}
		}
		// then count positive ones
		ULONG all_pos = all;
		for (size_t i = 0; i < m.size(); ++i) {
			if (all > 0 && fcmp_gt(w_pos[i], 0.)) {
				vspSize[i] = static_cast<ULONG>(all_pos * w_pos[i] / overall_pos);
				DBG_ASSERT(all >= vspSize[i], ValueError,
					"Not enough offspring to accommodate specified weight scheme. "
					"Current item is " + toStr(i) + " requires " + toStr(vspSize[i])
					+ " available " + toStr(all_pos) + ".");
				all -= vspSize[i];
			}
		}
		DBG_FAILIF(fcmp_eq(overall_pos, 0) && all > 0, ValueError,
			"An exact (all negative) weight system is used, but does not fill offspring subpopulation.");

		// individuals left by floating point calculation is added to
		// the last non-zero, positive weight virtual subpopulation.
		if (all > 0) {
			for (size_t i = m.size() - 1; i >= 0; --i)
				if (vspSize[i] != 0 && w_pos[i] > 0) {
					vspSize[i] += all;
					break;
				}
		}
		DBG_DO(DBG_DEVEL, cout << "VSP sizes in subpop " << sp << " is "
			                   << vspSize << endl);

		// it points to the first mating scheme.
		it = m.begin();
		it_end = m.end();
		vectorlu::iterator itSize = vspSize.begin();
		RawIndIterator ind = scratch.rawIndBegin(sp);
		DBG_FAILIF(pop.hasActivatedVirtualSubPop(sp), ValueError,
			"Subpopulation " + toStr(sp) + " has activated virtual subpopulation.");
		for (; it != it_end; ++it, ++itSize) {
			if ((*it)->virtualSubPop() != InvalidSubPopID)
				pop.activateVirtualSubPop(vspID(sp, (*it)->virtualSubPop()));
			// if previous mating scheme works on a virtual subpop,
			// and the current one is not. deactivate it.
			else if (pop.hasActivatedVirtualSubPop(sp))
				pop.deactivateVirtualSubPop(sp);

			// real mating
			try {
				if (!(*it)->mateSubPop(pop, sp, ind, ind + *itSize, ops))
					return false;
			} catch (...) {
				cout << "Mating in subpopulation " + toStr(sp) + " failed" << endl;
				throw;
			}
			ind += *itSize;
		}
		DBG_ASSERT(ind == scratch.rawIndEnd(sp), SystemError,
			"Mating scheme somehow does not fill the whole offspring population.");
		pop.deactivateVirtualSubPop(sp);
		// if more than two mating schemes working on the same subpopulation,
		// it is better to shuffle offspring afterwards,
		if (m.size() > 1 && m_shuffleOffspring) {
			DBG_DO(DBG_MATING, cout << "Random shuffle individuals in the offspring generation." << endl);
			std::random_shuffle(scratch.rawIndBegin(sp), scratch.rawIndEnd(sp));
			scratch.setIndOrdered(false);
		}
	}                         // each subpopulation.
	submitScratch(pop, scratch);
	return true;
}


}
