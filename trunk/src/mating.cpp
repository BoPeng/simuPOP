/**
 *  $File: mating.cpp $
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

#include "mating.h"

#if TR1_SUPPORT == 0
#  include <map>
typedef std::map<ULONG, simuPOP::Individual *> IdMap;
#elif TR1_SUPPORT == 1
#  include <unordered_map>
typedef std::tr1::unordered_map<ULONG, simuPOP::Individual *> IdMap;
#else
#  include <tr1/unordered_map>
typedef std::tr1::unordered_map<ULONG, simuPOP::Individual *> IdMap;
#endif

namespace simuPOP {


OffspringGenerator::OffspringGenerator(const opList & ops,
	const floatListFunc & numOffspring, const floatListFunc & sexMode) :
	m_numOffModel(NULL), m_sexModel(NULL), m_transmitters(ops), m_initialized(false)
{

	if (numOffspring.size() == 0) {
		DBG_FAILIF(!numOffspring.func().isValid(), ValueError,
			"Please specify either number of offspring or a function.");
		m_numOffModel = new FuncNumOffModel(numOffspring.func());
	} else if (numOffspring.size() == 1) {
		DBG_FAILIF(static_cast<int>(numOffspring[0]) <= 0, ValueError,
			"Number of offspring has to be positive.");
		m_numOffModel = new ConstNumOffModel(static_cast<UINT>(numOffspring[0]));
	} else if (numOffspring.size() > 1) {
		int mode = static_cast<int>(numOffspring[0]);

		if (mode == GEOMETRIC_DISTRIBUTION) {
			DBG_FAILIF(numOffspring.size() < 2 || fcmp_lt(numOffspring[1], 0) || fcmp_gt(numOffspring[1], 1.),
				ValueError, "P for a geometric distribution should be within [0,1]");
			m_numOffModel = new GeometricNumOffModel(numOffspring[1]);
		} else if (mode == POISSON_DISTRIBUTION) {
			m_numOffModel = new PoissonNumOffModel(numOffspring[1]);
		} else if (mode == BINOMIAL_DISTRIBUTION) {
			DBG_FAILIF(numOffspring.size() < 3 || static_cast<UINT>(numOffspring[2]) <= 1.,
				ValueError, "If mode is BINOMIAL_DISTRIBUTION, the second parameter should be greater than 1");
			DBG_FAILIF(numOffspring.size() < 2 || fcmp_le(numOffspring[1], 0) || fcmp_gt(numOffspring[1], 1.),
				ValueError, "P for a Bionomial distribution should be within (0,1].");
			DBG_FAILIF(numOffspring.size() < 3 || numOffspring[2] < 1,
				ValueError, "Max number of offspring in a binomial distribution should be greater than or equal to 1.");
			m_numOffModel = new BinomialNumOffModel(numOffspring[2], numOffspring[1]);
		} else if (mode == UNIFORM_DISTRIBUTION) {
			DBG_FAILIF(numOffspring.size() < 3 || numOffspring[2] < static_cast<UINT>(numOffspring[1]),
				ValueError, "If mode is UNIFORM_DISTRIBUTION, numOffspringParam should be greater than numOffspring");
			m_numOffModel = new UniformNumOffModel(numOffspring[1], numOffspring[2]);
		} else {
			throw ValueError("Wrong mating numoffspring mode. Should be one of \n"
				             "NumOffspring, NumOffspringEachFamily and GEometricDistribution");
		}
	}
	if (sexMode.size() == 0) {
		DBG_FAILIF(!sexMode.func().isValid(), ValueError,
			"Please specify one of the sex modes or a Python function.");
		m_sexModel = new FuncSexModel(sexMode.func());
	} else {
		int mode = static_cast<int>(sexMode[0]);

		if (mode == NO_SEX) {
			DBG_FAILIF(sexMode.size() != 1, ValueError, "No parameter is allowed for NO_SEX mode");
			m_sexModel = new NoSexModel();
		} else if (mode == RANDOM_SEX) {
			DBG_FAILIF(sexMode.size() != 1, ValueError, "No parameter is allowed for RANDOM_SEX mode");
			m_sexModel = new RandomSexModel();
		} else if (mode == PROB_OF_MALES) {
			DBG_FAILIF(sexMode.size() != 2, ValueError, "A parameter is required for mode PROB_OF_MALES.");
			DBG_FAILIF(fcmp_lt(sexMode[1], 0) || fcmp_gt(sexMode[1], 1),
				ValueError, "Probability of male has to be between 0 and 1");
			m_sexModel = new ProbOfMalesSexModel(sexMode[1]);
		} else if (mode == NUM_OF_MALES) {
			DBG_FAILIF(sexMode.size() != 2, ValueError, "A parameter is required for mode NUM_OF_MALES.");
			m_sexModel = new NumOfMalesSexModel(static_cast<UINT>(sexMode[1]));
		} else if (mode == NUM_OF_FEMALES) {
			DBG_FAILIF(sexMode.size() != 2, ValueError, "A parameter is required for mode NUM_OF_FEMALES.");
			m_sexModel = new NumOfFemalesSexModel(static_cast<UINT>(sexMode[1]));
		} else if (mode == SEQUENCE_OF_SEX) {
			DBG_FAILIF(sexMode.size() <= 2, ValueError, "A sequence of sex is required for mode SEQUENCE_OF_SEX.");
			m_sexModel = new SeqSexModel(sexMode.elems());
		} else if (mode == GLOBAL_SEQUENCE_OF_SEX) {
			DBG_FAILIF(sexMode.size() <= 2, ValueError, "A sequence of sex is required for mode GLOBAL_SEQUENCE_OF_SEX.");
			m_sexModel = new GlobalSeqSexModel(sexMode.elems());
		} else {
			DBG_FAILIF(true, ValueError, "Unrecognized sexMode.");
		}
	}
}


ULONG OffspringGenerator::numOffspring(int gen)
{
	return m_numOffModel->getNumOff(gen);

}


Sex OffspringGenerator::getSex(UINT count)
{
	return m_sexModel->getSex(count);
}


void OffspringGenerator::initialize(const Population & pop, SubPopID subPop)
{
	m_initialized = true;
}


string OffspringGenerator::describe(bool format) const
{
	string desc = "<simuPOP.OffspringGenerator> produces offspring using operators\n<ul>\n";
	opList::const_iterator iop = m_transmitters.begin();
	opList::const_iterator iopEnd = m_transmitters.end();

	for (; iop != iopEnd; ++iop)
		desc += "<li>" + (*iop)->describe(false) + " " + (*iop)->applicability() + "\n";
	desc += "</ul>\n";
	return format ? formatDescription(desc) : desc;
}


UINT OffspringGenerator::generateOffspring(Population & pop, Individual * dad, Individual * mom,
                                           RawIndIterator & it,
                                           RawIndIterator & itEnd)
{
	DBG_ASSERT(initialized(), ValueError,
		"Offspring generator is not initialized before used to generate offspring");

	// generate numOff offspring per mating, or until it  reaches offEnd
	UINT count = 0;
	bool accept = true;
	UINT numOff = numOffspring(pop.gen());
	UINT attempt = 0;
	while (attempt < numOff && it != itEnd) {
		// not all families have the same size because some offspring
		// may be discarded (count).
		++attempt;

		// set sex, during mating operator will try to
		// follow the offspring sex (e.g. pass X or Y chromosome)
		it->setSex(getSex(count));
		//
		accept = true;
		opList::const_iterator iop = m_transmitters.begin();
		opList::const_iterator iopEnd = m_transmitters.end();
		for (; iop != iopEnd; ++iop) {
			try {
				if (!(*iop)->isActive(pop.rep(), pop.gen()))
					continue;
				if (!(*iop)->applyDuringMating(pop, it, dad, mom)) {
					accept = false;
					break;
				}
			} catch (Exception e) {
				cerr << "One of the transmitters " << (*iop)->describe()
				     << " throws an exception.\n" << e.message() << "\n" << endl;
				throw e;
			}
		}

		if (accept) {
			++it;
			++count;
		}
	}
	return count;
}


ControlledOffspringGenerator::ControlledOffspringGenerator(
	const uintList & loci, const uintList & alleles, PyObject * freqFunc,
	const opList & ops, const floatListFunc & numOffspring,
	const floatListFunc & sexMode)
	: OffspringGenerator(ops, numOffspring, sexMode),
	m_loci(loci.elems()), m_alleles(alleles.elems()), m_freqFunc(freqFunc),
	m_expAlleles(), m_totAllele(), m_curAllele()

{
	if (!m_freqFunc.isValid())
		throw ValueError("Please specify a valid frequency function");
}


ControlledOffspringGenerator::ControlledOffspringGenerator(const ControlledOffspringGenerator & rhs)
	: OffspringGenerator(rhs),
	m_loci(rhs.m_loci),
	m_alleles(rhs.m_alleles),
	m_freqFunc(rhs.m_freqFunc)
{
}


string ControlledOffspringGenerator::describe(bool format) const
{
	string desc = "<simuPOP.ControlledOffspringGenerator> produces offspring using operators\n<ul>\n";
	opList::const_iterator iop = m_transmitters.begin();
	opList::const_iterator iopEnd = m_transmitters.end();

	for (; iop != iopEnd; ++iop)
		desc += "<li>" + (*iop)->describe(false) + " " + (*iop)->applicability() + "\n";
	desc += "</ul>\nwhile controlling allele frequency";
	return format ? formatDescription(desc) : desc;
}


void ControlledOffspringGenerator::getExpectedAlleles(const Population & pop,
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
			UINT locus = m_loci[i];
			Allele allele = ToAllele(m_alleles[i]);

			// determine the number alleles at each subpopulation.
			vectorf curFreq(numSP, 0);
			bool hasAllele = false;
			for (size_t sp = 0, n = 0; sp < numSP; ++sp) {
				IndAlleleIterator a = const_cast<Population &>(pop).alleleIterator(locus, sp);
				for (; a.valid(); ++a)
					if (AlleleUnsigned(*a) == allele)
						n++;
				hasAllele = hasAllele || n > 0;
				curFreq[sp] = double(n) / (pop.subPopSize(sp) * pop.ploidy());
			}

			DBG_DO(DBG_MATING, cerr << "Current frequency at locus " << locus
				                    << " is " << curFreq << endl);

			// if there is no alleles
			if (!hasAllele && expFreq[i] > 0.)
				throw RuntimeError("No disease allele exists at generation "
					+ toStr(pop.gen()) + ", but expected allele frequency at locus " + toStr(locus) + " is greater than 0.");

			DBG_WARNIF(hasAllele && fcmp_eq(expFreq[i], 0.), "Disease allele exists at generation "
				+ toStr(pop.gen()) + ", but expected allele frequency is zero.");

			// calculate exp number of affected offspring in the next generation.
			//
			// step 1: totalsize*expFreq is the total number of disease alleles
			// step 2: assign these alleles to each subpopulation according to a multi-nomial
			// distribution with p_i beging allele frequency at each subpopulation.
			// assign these numbers to each subpopulation
			vectoru counts = getRNG().randMultinomial(static_cast<unsigned int>(pop.popSize() * expFreq[i] * pop.ploidy()),
				curFreq);
			for (size_t sp = 0; sp < counts.size(); ++sp)
				m_expAlleles[numSP * i + sp] = counts[sp];
		}
	} else {
		// simpler case, one subpopulation, or with gieven allele frequency
		for (size_t i = 0; i < nLoci; ++i) {
			for (size_t sp = 0; sp < numSP; ++sp) {
#ifndef OPTIMIZED
				UINT locus = m_loci[i];
				Allele allele = ToAllele(m_alleles[i]);
				ULONG n = 0;
				// go through all alleles
				IndAlleleIterator a = const_cast<Population &>(pop).alleleIterator(locus, sp);
				for (; a.valid(); ++a) {
					if (AlleleUnsigned(*a) == allele)
						n++;
				}

				// if there is no alleles
				if (n == 0 && expFreq[sp * nLoci + i] > 0.)
					throw RuntimeError("No disease allele exists at generation " +
						toStr(pop.gen()) + " but exp allele frequency at locus " + toStr(locus) +
						" in subpopulation " + toStr(sp) + " is greater than 0.");
#endif
				m_expAlleles[numSP * i + sp] = static_cast<UINT>(pop.subPopSize(sp) * pop.ploidy() * expFreq[sp * nLoci + i]);
				if (expFreq[sp * nLoci + i] > 0. && m_expAlleles[numSP * i + sp] == 0)
					m_expAlleles[numSP * i + sp] = 1;
			}
		}
	}
}


void ControlledOffspringGenerator::initialize(const Population & pop, SubPopID subPop)
{
	OffspringGenerator::initialize(pop, subPop);

	// expected frequency at each locus
	if (subPop == 0) {
		vectorf expFreq = m_freqFunc(PyObj_As_Array, "(i)", pop.gen());
		DBG_DO(DBG_MATING, cerr << "expected freq " << expFreq << endl);

		//
		// determine expected number of alleles of each allele
		// at each subpopulation.
		getExpectedAlleles(pop, expFreq);
		DBG_DO(DBG_MATING, cerr << "expected alleles " << m_expAlleles << endl);
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
			cerr << "Warning: number of planned affected alleles exceed population size.";
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
UINT ControlledOffspringGenerator::generateOffspring(Population & pop, Individual * dad, Individual * mom,
                                                     RawIndIterator & offBegin,
                                                     RawIndIterator & offEnd)
{
	UINT nLoci = m_loci.size();
	//
	// generate m_numOffspring offspring per mating
	// record family size (this may be wrong for the last family)
	//
	RawIndIterator itBegin = offBegin;
	UINT numOff = OffspringGenerator::generateOffspring(pop,
		dad, mom, offBegin, offEnd);

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
			if (m_flip[i] ? (*ptr != ToAllele(m_alleles[i]))
				: (*ptr == ToAllele(m_alleles[i]))) {
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
		// DBG_DO(DBG_MATING, cerr << "Reject " << na << endl);
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


void SequentialParentChooser::initialize(Population & pop, SubPopID sp)
{
	m_begin = pop.indIterator(sp);
	m_ind = m_begin;
	m_initialized = true;
}


ParentChooser::IndividualPair SequentialParentChooser::chooseParents(RawIndIterator)
{
	if (!m_ind.valid())
		m_ind = m_begin;
	return ParentChooser::IndividualPair(&*m_ind++, NULL);
}


void SequentialParentsChooser::initialize(Population & pop, SubPopID subPop)
{
	m_numMale = 0;
	m_numFemale = 0;
	m_curMale = 0;
	m_curFemale = 0;
	m_maleIndex.clear();
	m_femaleIndex.clear();

	IndIterator it = pop.indIterator(subPop);
	for (; it.valid(); ++it) {
		if (it->sex() == MALE) {
			m_numMale++;
			m_maleIndex.push_back(it.rawIter());
		} else {
			m_numFemale++;
			m_femaleIndex.push_back(it.rawIter());
		}
	}
	m_initialized = true;
}


ParentChooser::IndividualPair SequentialParentsChooser::chooseParents(RawIndIterator)
{
	DBG_ASSERT(initialized(), SystemError,
		"Please initialize this parent chooser before using it");

	Individual * dad = NULL;
	Individual * mom = NULL;

	if (m_curMale == m_numMale)
		m_curMale = 0;
	if (m_curFemale == m_numFemale)
		m_curFemale = 0;

	// using weighted sampler.
	if (m_numMale != 0)
		dad = &*(m_maleIndex[m_curMale++]);
	else
		dad = &*(m_femaleIndex[m_curFemale++]);

	if (m_numFemale != 0)
		mom = &*(m_femaleIndex[m_curFemale++]);
	else
		mom = &*(m_maleIndex[m_curMale++]);
	return std::make_pair(dad, mom);
}


void RandomParentChooser::initialize(Population & pop, SubPopID sp)
{
	m_index.clear();

	// In a virtual subpopulation, because m_begin + ... is **really** slow
	// It is a good idea to cache IndIterators. This is however inefficient
	// for non-virtual populations
	if (pop.hasActivatedVirtualSubPop(sp) && !m_replacement) {
		IndIterator it = pop.indIterator(sp);
		for (; it.valid(); ++it)
			m_index.push_back(it.rawIter());
	}

	m_selection = m_replacement && pop.hasInfoField(m_selectionField);
	if (m_selection) {
		UINT fit_id = pop.infoIdx(m_selectionField);
		// regardless of sex, get fitness for everyone.
		m_sampler.set(vectorf(pop.infoBegin(fit_id, sp),
				pop.infoEnd(fit_id, sp)));
	} else {
		m_size = m_index.size();
		if (m_size == 0)         // if m_index is not used (no VSP)
			m_size = pop.subPopSize(sp);
	}

	if (!m_replacement)
		getRNG().randomShuffle(m_index.begin(), m_index.end());

	m_shift = pop.subPopBegin(sp);
	m_initialized = true;
}


ParentChooser::IndividualPair RandomParentChooser::chooseParents(RawIndIterator basePtr)
{
	DBG_ASSERT(initialized(), SystemError,
		"Please initialize this parent chooser before using it");
	// choose a parent
	if (!m_replacement) {
		if (m_index.empty())
			throw RuntimeError("All parents have been chosen.");
		Individual * ind = &*m_index.back();
		m_index.pop_back();
		return IndividualPair(ind, NULL);
	}
	Individual * ind = NULL;
	if (m_index.empty()) {
		if (m_selection)
			// basePtr points to the beginning of the population, not subpopulation
			ind = &*(basePtr + m_shift + m_sampler.draw());
		else
			ind = &*(basePtr + m_shift + getRNG().randInt(m_size));
	} else {
		if (m_selection)
			ind = &*(m_index[m_sampler.draw()]);
		else
			ind = &*(m_index[getRNG().randInt(m_size)]);
	}
	return IndividualPair(ind, NULL);
}


void RandomParentsChooser::initialize(Population & pop, SubPopID subPop)
{
	m_numMale = 0;
	m_numFemale = 0;

	IndIterator it = pop.indIterator(subPop);
	for (; it.valid(); ++it) {
		if (it->sex() == MALE)
			m_numMale++;
		else
			m_numFemale++;
	}

	// allocate memory at first for performance reasons
	m_maleIndex.resize(m_numMale);
	m_femaleIndex.resize(m_numFemale);

	m_selection = m_replacement && pop.hasInfoField(m_selectionField);
	UINT fit_id = 0;
	if (m_selection) {
		fit_id = pop.infoIdx(m_selectionField);
		m_maleFitness.resize(m_numMale);
		m_femaleFitness.resize(m_numFemale);
	}

	m_numMale = 0;
	m_numFemale = 0;

	it = pop.indIterator(subPop);
	for (; it.valid(); it++) {
		if (it->sex() == MALE) {
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
		getRNG().randomShuffle(m_maleIndex.begin(), m_maleIndex.end());
		getRNG().randomShuffle(m_femaleIndex.begin(), m_femaleIndex.end());
	}

	if (m_selection) {
		m_malesampler.set(m_maleFitness);
		m_femalesampler.set(m_femaleFitness);
		DBG_DO(DBG_DEVEL, cerr << "Male fitness " << m_maleFitness << endl);
		DBG_DO(DBG_DEVEL, cerr << "Female fitness " << m_femaleFitness << endl);
	}

	DBG_FAILIF(!m_replacement && m_selection, ValueError,
		"Selection is not allowed in random sample without replacement");

	m_initialized = true;
}


ParentChooser::IndividualPair RandomParentsChooser::chooseParents(RawIndIterator)
{
	DBG_ASSERT(initialized(), SystemError,
		"Please initialize this parent chooser before using it");

	Individual * dad = NULL;
	Individual * mom = NULL;

	if (!m_replacement) {
		if (m_femaleIndex.empty())
			throw ValueError("All females have been chosen.");
		mom = &*m_femaleIndex.back();
		m_femaleIndex.pop_back();

		if (m_maleIndex.empty())
			throw ValueError("All males have been chosen.");
		dad = &*m_maleIndex.back();
		m_maleIndex.pop_back();
		return std::make_pair(dad, mom);
	}

	// this exception should be raised also in optimized mode because the cause
	// can be random.
	if (m_numMale == 0)
		throw RuntimeError("RandomParentsChooser fails because there is no male individual in a subpopulation.");
	if (m_numFemale == 0)
		throw RuntimeError("RandomParentsChooser fails because there is no female individual in a subpopulation ");

	if (m_selection) {
		// using weighted sampler.
		dad = &*(m_maleIndex[m_malesampler.draw()]);
		mom = &*(m_femaleIndex[m_femalesampler.draw()]);
	} else {
		dad = &*(m_maleIndex[getRNG().randInt(m_numMale)]);
		mom = &*(m_femaleIndex[getRNG().randInt(m_numFemale)]);
	}
	return std::make_pair(dad, mom);
}


void PolyParentsChooser::initialize(Population & pop, SubPopID subPop)
{
	m_numMale = 0;
	m_numFemale = 0;
	m_polyCount = 0;

	IndIterator it = pop.indIterator(subPop);
	for (; it.valid(); ++it) {
		if (it->sex() == MALE)
			m_numMale++;
		else
			m_numFemale++;
	}

	// allocate memory at first for performance reasons
	m_maleIndex.resize(m_numMale);
	m_femaleIndex.resize(m_numFemale);

	m_selection = pop.hasInfoField(m_selectionField);

	UINT fit_id = 0;
	if (m_selection) {
		fit_id = pop.infoIdx(m_selectionField);
		m_maleFitness.resize(m_numMale);
		m_femaleFitness.resize(m_numFemale);
	}

	m_numMale = 0;
	m_numFemale = 0;

	it = pop.indIterator(subPop);
	for (; it.valid(); it++) {
		if (it->sex() == MALE) {
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
		DBG_DO(DBG_DEVEL, cerr << "Male fitness " << m_maleFitness << endl);
		DBG_DO(DBG_DEVEL, cerr << "Female fitness " << m_femaleFitness << endl);
	}

	m_initialized = true;
}


ParentChooser::IndividualPair PolyParentsChooser::chooseParents(RawIndIterator)
{
	DBG_ASSERT(initialized(), SystemError,
		"Please initialize this parent chooser before using it");

	Individual * dad = NULL;
	Individual * mom = NULL;

	if (m_polyNum > 1 && m_polyCount > 0) {
		if (m_polySex == MALE)
			dad = m_lastParent;
		else
			mom = m_lastParent;
		m_polyCount--;
	}

	// using weidhted sampler.
	if (dad == NULL) {
		if (m_numMale == 0)
			throw RuntimeError("PolyParentsChooser fails because there is no male individual in a subpopulation.");

		if (m_selection)
			dad = &*(m_maleIndex[m_malesampler.draw()]);
		else
			dad = &*(m_maleIndex[getRNG().randInt(m_numMale)]);

		if (m_polySex == MALE && m_polyNum > 1) {
			m_polyCount = m_polyNum - 1;
			m_lastParent = dad;
		}
	}

	if (mom == NULL) {
		if (m_numFemale == 0)
			throw RuntimeError("PolyParentsChooser fails because there is no female individual in a subpopulation.");

		if (m_selection)
			mom = &*(m_femaleIndex[m_femalesampler.draw()]);
		else
			mom = &*(m_femaleIndex[getRNG().randInt(m_numFemale)]);

		if (m_polySex == FEMALE && m_polyNum > 1) {
			m_polyCount = m_polyNum - 1;
			m_lastParent = mom;
		}
	}
	return std::make_pair(dad, mom);
}


/*
   void infoParentsChooser::initialize(Population & pop, SubPopID sp)
   {
    if (m_func.isValid()) {
        PyObject * popObj = pyPopObj(static_cast<void *>(&pop));
        // if pop is valid?
        if (popObj == NULL)
            throw SystemError("Could not pass Population to the provided function. \n"
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

    m_selection = pop.hasInfoField(m_selectionField);

    m_index.clear();

    // In a virtual subpopulation, because m_begin + ... is **really** slow
    // It is a good idea to cache IndIterators. This is however inefficient
    // for non-virtual populations
    IndIterator it = pop.indIterator(sp);
    vectorf fitness;
    UINT fit_id = 0;
    if (m_selection)
        fit_id = pop.infoIdx(m_selectionField);
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
    DBG_WARNIF(m_degenerate, "Parents are chosen randomly because there is no valid index.");
    if (m_degenerate) {
        for (it = pop.indIterator(sp); it.valid(); ++it) {
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
        getRNG().randomShuffle(m_index.begin(), m_index.end());

    DBG_FAILIF(!m_replacement && m_selection, ValueError,
        "Selection is not allowed in random sample without replacement");

    m_shift = pop.subPopBegin(sp);
    m_initialized = true;
   }
 */


/*
   ParentChooser::IndividualPair infoParentsChooser::chooseParents(RawIndIterator basePtr)
   {
    DBG_ASSERT(initialized(), SystemError,
        "Please initialize this parent chooser before using it");
    Individual * par1 = RandomParentChooser::chooseParents(basePtr).first;
    Sex sex1 = par1->sex();
    // there is no valid information field value
    if (m_degenerate) {
        int attempt = 0;
        while (++attempt < 1000) {
            Individual * par2 = RandomParentChooser::chooseParents(basePtr).first;
            if (par2->sex() != sex1)
                return sex1 == MALE ? std::make_pair(par1, par2) : std::make_pair(par2, par1);
        }
        throw RuntimeError("Can not locate any individual of opposite sex");
    }
    // the way this parent chooser is initialized guranttees that
    // theres is at lest one valid field.
    vector<Individual *> validInds;
    for (size_t i = 0; i < m_infoIdx.size(); ++i) {
        int info = par1->intInfo(m_infoIdx[i]);
        if (info < 0)
            continue;
        if (m_idIdx < 0) {
            RawIndIterator par2 = basePtr + info;
            if (par2->sex() != sex1)
                validInds.push_back(&*par2);
        }
    }
    DBG_FAILIF(validInds.empty(), SystemError, "No valid relative is found");
    Individual * par2 = validInds[getRNG().randInt(validInds.size())];
    DBG_DO(DBG_DEVEL, cerr	<< "infoParentsChooser: par1: " << par1 - &*basePtr
                            << " par2: " << par2 - &*basePtr << endl);
    return sex1 == MALE ? std::make_pair(par1, par2) : std::make_pair(par2, par1);
   }
 */


PyParentsChooser::PyParentsChooser(PyObject * pc)
	: ParentChooser(), m_func(pc), m_popObj(NULL),
	m_generator(NULL)
{
}


void PyParentsChooser::initialize(Population & pop, SubPopID sp)
{
#if PY_VERSION_HEX < 0x02040000
	throw SystemError("Your Python version does not have good support for generator"
		              " so this python parent chooser can not be used.");
#endif

	DBG_FAILIF(pop.hasActivatedVirtualSubPop(sp), ValueError,
		"Python parent chooser can not be used in a virtual subpopulation.");

#ifndef OPTIMIZED
	m_size = pop.subPopSize(sp);
#endif
	m_begin = pop.indIterator(sp);

	m_popObj = pyPopObj(static_cast<void *>(&pop));

	PyObject * args = PyTuple_New(m_func.numArgs());
	DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");
	for (int i = 0; i < m_func.numArgs(); ++i) {
		const string & arg = m_func.arg(i);
		if (arg == "pop")
			PyTuple_SET_ITEM(args, i, m_popObj);
		else if (arg == "subPop")
			PyTuple_SET_ITEM(args, i, PyInt_FromLong(sp));
		else {
			DBG_FAILIF(true, ValueError,
				"Only parameters 'pop' and 'subPop' are acceptable in a generator function.");
		}
	}
	m_generator.set(m_func(args));
	Py_DECREF(args);
	m_initialized = true;
}


ParentChooser::IndividualPair PyParentsChooser::chooseParents(RawIndIterator)
{
	DBG_ASSERT(initialized(), SystemError,
		"Please initialize this parent chooser before using it");

	PyObject * item = m_generator.next();

#ifndef OPTIMIZED
	if (item == NULL && debug(DBG_GENERAL)) {
		PyErr_Print();
		PyErr_Clear();
	}
	DBG_FAILIF(item == NULL, ValueError,
		"User-defined function yield invalid value.");
#endif


	if (PyInt_Check(item) || PyLong_Check(item)) {
		long int parent;
		PyObj_As_Int(item, parent);
#ifndef OPTIMIZED
		DBG_ASSERT(static_cast<unsigned>(parent) < m_size,
			ValueError, "Returned index (" + toStr(parent) +
			") is greater than subpopulation size " + toStr(m_size));
#endif
		Py_DECREF(item);
		return ParentChooser::IndividualPair(&*(m_begin + parent), NULL);
	} else if (PySequence_Check(item)) {
		DBG_ASSERT(PySequence_Size(item) == 2, RuntimeError,
			"Parents should be returned in the form of a sequence of two elements");

		Individual * parents[2];
		for (size_t i = 0; i < 2; ++i) {
			PyObject * v = PySequence_GetItem(item, i);
			if (PyInt_Check(v) || PyLong_Check(v)) {
				ULONG idx = PyInt_AsLong(v);
				DBG_ASSERT(idx < m_size, ValueError, "Returned parent index (" + toStr(idx) +
					") is greater than subpopulation size " + toStr(m_size));
				parents[i] = &*(m_begin + idx);
			} else {
				void * ind = pyIndPointer(v);
				if (ind)
					parents[i] = reinterpret_cast<Individual *>(ind);
				else
					DBG_ASSERT(false, ValueError, "Invalid type of returned parent.");
			}
			Py_DECREF(v);
		}
		Py_DECREF(item);
		return ParentChooser::IndividualPair(parents[0], parents[1]);
	} else {
		// is an individual object is returned?
		void * ind = pyIndPointer(item);
		if (ind)
			return ParentChooser::IndividualPair(reinterpret_cast<Individual *>(pyIndPointer(item)), NULL);
		else
			DBG_ASSERT(false, ValueError, "Invalid type of returned parent.");
	}
	// this should not be reached
	return ParentChooser::IndividualPair(NULL, NULL);
}


void PyParentsChooser::finalize(Population & pop, SubPopID sp)
{
	DBG_FAILIF(m_popObj == NULL, SystemError, "Python generator is not properly initialized.");
	Py_DECREF(m_popObj);
	m_generator.set(NULL);
	m_popObj = NULL;
	m_initialized = false;
}


MatingScheme::MatingScheme(const uintListFunc & subPopSize)
	: m_subPopSize(subPopSize)
{
}


bool MatingScheme::prepareScratchPop(Population & pop, Population & scratch)
{
	if (scratch.genoStruIdx() != pop.genoStruIdx())
		scratch.fitGenoStru(pop.genoStruIdx());

	// use population structure of pop
	if (m_subPopSize.empty() && !m_subPopSize.func().isValid())
		scratch.fitSubPopStru(pop.subPopSizes(), pop.subPopNames());
	else if (!m_subPopSize.empty())                                                     // set subPoplation size
		scratch.fitSubPopStru(m_subPopSize.elems(), pop.subPopNames());
	else {                                                                              // use m_subPopSizeFunc
		const pyFunc & func = m_subPopSize.func();
		PyObject * args = PyTuple_New(func.numArgs());
		DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

		for (int i = 0; i < func.numArgs(); ++i) {
			const string & arg = func.arg(i);
			if (arg == "gen")
				PyTuple_SET_ITEM(args, i, PyInt_FromLong(pop.gen()));
			else if (arg == "pop")
				PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(&pop)));
			else {
				DBG_FAILIF(true, ValueError,
					"Only parameters 'gen' and 'pop' are acceptable in a demographic function.");
			}
		}
		vectori res = func(PyObj_As_IntArray, args);
		Py_XDECREF(args);

		vectoru sz(res.size());
		for (size_t i = 0; i < res.size(); i++)
			sz[i] = static_cast<ULONG>(res[i]);

		// allow change of pop size of scratch
		scratch.fitSubPopStru(sz, pop.subPopNames());
	}
	// this is not absolutely necessary but will reduce confusions
	scratch.setVirtualSplitter(pop.virtualSplitter());
	scratch.clearInfo();
	DBG_DO(DBG_SIMULATOR, cerr << "New subpop size " << scratch.subPopSizes() << endl);

	DBG_FAILIF(scratch.numSubPop() != pop.numSubPop(),
		ValueError, "number of subPopulaitons must agree.\n Pre: "
		+ toStr(pop.numSubPop()) + " now: " + toStr(scratch.numSubPop()));
	return true;
}


bool MatingScheme::mate(Population & pop, Population & scratch)
{
	// scrtach will have the right structure.
	if (!prepareScratchPop(pop, scratch))
		return false;

	for (SubPopID sp = 0; sp < static_cast<SubPopID>(pop.numSubPop()); ++sp)
		if (!mateSubPop(pop, sp, scratch.rawIndBegin(sp), scratch.rawIndEnd(sp)))
			return false;
	submitScratch(pop, scratch);
	return true;
}


void MatingScheme::submitScratch(Population & pop, Population & scratch)
{
	// use scratch population,
	pop.push(scratch);
	scratch.validate("after push and discard");
}


HomoMating::HomoMating(ParentChooser & chooser,
	OffspringGenerator & generator,
	const uintListFunc & subPopSize,
	subPopList subPops, double weight)
	: MatingScheme(subPopSize), m_subPops(subPops), m_weight(weight)
{
	m_ParentChooser = chooser.clone();
	m_OffspringGenerator = generator.clone();
}


string HomoMating::describe(bool format) const
{
	string desc = "<simuPOP.HomoMating> a homogeneous mating scheme that uses\n<ul>\n<li>"
	              + m_ParentChooser->describe(false) + "\n<li>"
	              + m_OffspringGenerator->describe(false) + "</ul>\n";

	return format ? formatDescription(desc) : desc;
}


bool HomoMating::mateSubPop(Population & pop, SubPopID subPop,
                            RawIndIterator offBegin, RawIndIterator offEnd)
{
	// nothing to do.
	if (offBegin == offEnd)
		return true;

	if (!m_ParentChooser->initialized())
		m_ParentChooser->initialize(pop, subPop);

	if (!m_OffspringGenerator->initialized())
		m_OffspringGenerator->initialize(pop, subPop);

	// generate scratch.subPopSize(sp) individuals.
	RawIndIterator it = offBegin;
	while (it != offEnd) {
		Individual * dad = NULL;
		Individual * mom = NULL;
		ParentChooser::IndividualPair const parents = m_ParentChooser->chooseParents(pop.rawIndBegin());
		dad = parents.first;
		mom = parents.second;

		//
		UINT numOff = m_OffspringGenerator->generateOffspring(pop, dad, mom, it, offEnd);
		(void)numOff;             // silent warning about unused variable.
	}
	m_ParentChooser->finalize(pop, subPop);
	m_OffspringGenerator->finalize(pop);
	return true;
}


bool PedigreeMating::mate(Population & pop, Population & scratch)
{
	if (m_gen == -1)
		return false;

	// scrtach will have the right structure.
	if (scratch.genoStruIdx() != pop.genoStruIdx())
		scratch.fitGenoStru(pop.genoStruIdx());
	//
	UINT oldGen = m_ped.curAncestralGen();
	const_cast<Pedigree &>(m_ped).useAncestralGen(m_gen);
	DBG_DO(DBG_MATING, cerr << "Producing offspring generation of size " << m_ped.subPopSizes() <<
		" using generation " << m_gen << " of the pedigree." << endl);
	scratch.fitSubPopStru(m_ped.subPopSizes(), m_ped.subPopNames());
	scratch.setVirtualSplitter(pop.virtualSplitter());
	scratch.clearInfo();

	// build an index for parents
	IdMap idMap;
	UINT idIdx = pop.infoIdx(m_idField);
	RawIndIterator it = pop.rawIndBegin();
	RawIndIterator it_end = pop.rawIndEnd();
	for (; it != it_end; ++it)
		idMap[toID(it->info(idIdx))] = &*it;

	it = scratch.rawIndBegin();
	it_end = scratch.rawIndEnd();
	for (size_t i = 0; it != it_end; ++it, ++i) {
		const Individual & pedInd = m_ped.individual(i);

		ULONG my_id = toID(pedInd.info(m_ped.idIdx()));
		ULONG father_id = m_ped.fatherOf(my_id);
		ULONG mother_id = m_ped.motherOf(my_id);
		Individual * dad = NULL;
		Individual * mom = NULL;

		if (father_id) {
			IdMap::iterator dad_it = idMap.find(father_id);
			DBG_FAILIF(dad_it == idMap.end(), RuntimeError,
				"Could not locate individual with ID " + toStr(father_id));
			dad = &*(dad_it->second);
		}
		if (mother_id) {
			IdMap::iterator mom_it = idMap.find(mother_id);
			DBG_FAILIF(mom_it == idMap.end(), RuntimeError,
				"Could not locate individual with ID " + toStr(mother_id));
			mom = &*(mom_it->second);
		}
		DBG_DO(DBG_MATING, cerr << "Choosing parents " << father_id << " and "
			                    << mother_id << " for offspring " << my_id << endl);

		// copy sex
		it->setSex(pedInd.sex());
		// copy id
		it->setInfo(my_id, m_idField);
		//
		opList::const_iterator iop = m_transmitters.begin();
		opList::const_iterator iopEnd = m_transmitters.end();
		for (; iop != iopEnd; ++iop) {
			try {
				if ((*iop)->isActive(pop.rep(), pop.gen()))
					(*iop)->applyDuringMating(pop, it, dad, mom);
			} catch (Exception e) {
				cerr << "One of the transmitters " << (*iop)->describe()
				     << " throws an exception.\n" << e.message() << "\n" << endl;
				throw e;
			}
		}
		// copy individual ID again, just to make sure that even if during mating operators
		// changes ID, pedigree mating could proceed normally.
		it->setInfo(my_id, m_idField);
	}
	const_cast<Pedigree &>(m_ped).useAncestralGen(oldGen);
	submitScratch(pop, scratch);
	--m_gen;
	return true;
}


string PedigreeMating::describe(bool format) const
{
	string desc = "<simuPOP.PedigreeMating> evolves a population following a pedigree, using operators\n<ul>\n";
	opList::const_iterator iop = m_transmitters.begin();
	opList::const_iterator iopEnd = m_transmitters.end();

	for (; iop != iopEnd; ++iop)
		desc += "<li>" + (*iop)->describe(false) + " " + (*iop)->applicability() + "\n";
	desc += "</ul>\n";
	return format ? formatDescription(desc) : desc;
}


HeteroMating::HeteroMating(const vectormating & matingSchemes,
	const uintListFunc & subPopSize,
	bool shuffleOffspring)
	: MatingScheme(subPopSize),
	m_shuffleOffspring(shuffleOffspring)
{
	vectormating::const_iterator it = matingSchemes.begin();
	vectormating::const_iterator it_end = matingSchemes.end();

	for (; it != it_end; ++it)
		m_matingSchemes.push_back(dynamic_cast<HomoMating *>((*it)->clone()));
}


string HeteroMating::describe(bool format) const
{
	string desc = "<simuPOP.HeteroMating> a heterogeneous mating scheme with " +
	              toStr(m_matingSchemes.size()) + " homogeneous mating schemes:\n<ul>\n";
	vectormating::const_iterator it = m_matingSchemes.begin();
	vectormating::const_iterator it_end = m_matingSchemes.end();

	for (; it != it_end; ++it) {
		desc += "<li>" + dynamic_cast<HomoMating *>(*it)->describe(false) +
		        "<indent>in ";
		subPopList subPops = (*it)->subPops();
		if (subPops.allAvail())
			desc += "all subpopulations.\n";
		else {
			desc += "subpopulations ";
			for (size_t i = 0; i < subPops.size(); ++i) {
				vspID sp = subPops[i];
				if (i != 0)
					desc += ", ";
				if (sp.isVirtual()) {
					desc += "(" + (sp.allAvailSP() ? "ALL_AVAIL" : toStr(sp.subPop())) + ", " +
					        (sp.allAvailVSP() ? "ALL_AVAIL" : toStr(sp.virtualSubPop())) + ")";
				} else
					desc += toStr(sp.subPop());
			}
			desc += ".\n";
		}
	}
	desc += "</ul>\n";
	return format ? formatDescription(desc) : desc;
}


HeteroMating::~HeteroMating()
{
	vectormating::iterator it = m_matingSchemes.begin();
	vectormating::iterator it_end = m_matingSchemes.end();

	for (; it != it_end; ++it)
		delete *it;
}


HeteroMating::HeteroMating(const HeteroMating & rhs) :
	MatingScheme(rhs), m_shuffleOffspring(rhs.m_shuffleOffspring)
{
	vectormating::const_iterator it = rhs.m_matingSchemes.begin();
	vectormating::const_iterator it_end = rhs.m_matingSchemes.end();

	for (; it != it_end; ++it) {
		m_matingSchemes.push_back(dynamic_cast<HomoMating *>((*it)->clone()));
		DBG_WARNIF(dynamic_cast<HomoMating *>(*it)->subPopSizeSpecified(),
			"Parameter subPopSize of a HomoMating is ignored when this mating"
			" scheme is used in a heterogeneous mating scheme.");
	}
}


bool HeteroMating::mate(Population & pop, Population & scratch)
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
		subPopList sps;                         // each subpopulations
		//
		vectormating::iterator it = m_matingSchemes.begin();
		vectormating::iterator it_end = m_matingSchemes.end();
		for (; it != it_end; ++it) {
			subPopList subPops = (*it)->subPops().expandFrom(pop);
			subPopList::const_iterator vsp = subPops.begin();
			subPopList::const_iterator vspEnd = subPops.end();
			for (; vsp != vspEnd; ++vsp) {
				if (vsp->subPop() != sp)
					continue;
				// if it is used for this subpop, or all subpopulations
				m.push_back(*it);
				sps.push_back(*vsp);
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
					w_pos[i] = pop.subPopSize(sps[i]);
		}
		DBG_DO(DBG_DEVEL, cerr << "Positive mating scheme weights: " << w_pos << '\n'
			                   << "Negative mating scheme weights: " << w_neg << endl);

		// weight.
		double overall_pos = std::accumulate(w_pos.begin(), w_pos.end(), 0.);
		double overall_neg = std::accumulate(w_neg.begin(), w_neg.end(), 0.);
		(void)overall_neg;                         // silent warning about unused variable.
		DBG_FAILIF(fcmp_eq(overall_pos, 0.) && fcmp_eq(overall_neg, 0.), ValueError,
			"Overall weight is zero");
		//
		vectoru vspSize(m.size());
		ULONG all = scratch.subPopSize(sp);
		// first count negative ones
		for (size_t i = 0; i < m.size(); ++i) {
			if (fcmp_gt(w_neg[i], 0.)) {
				vspSize[i] = static_cast<ULONG>(pop.subPopSize(sps[i]) * w_neg[i]);
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
		DBG_DO(DBG_DEVEL, cerr << "VSP sizes in subpop " << sp << " is "
			                   << vspSize << endl);

		DBG_ASSERT(vspSize.size() == m.size() && m.size() == sps.size(),
			SystemError, "Failed to determine subpopulation size");
		// it points to the first mating scheme.
		vectoru::iterator itSize = vspSize.begin();
		RawIndIterator ind = scratch.rawIndBegin(sp);
		DBG_FAILIF(pop.hasActivatedVirtualSubPop(sp), ValueError,
			"SubPopulation " + toStr(sp) + " has activated virtual subpopulation.");
		for (UINT idx = 0; idx < m.size(); ++idx, ++itSize) {
			if (sps[idx].isVirtual())
				pop.activateVirtualSubPop(sps[idx]);
			// if previous mating scheme works on a virtual subpop,
			// and the current one is not. deactivate it.
			else if (pop.hasActivatedVirtualSubPop(sp))
				pop.deactivateVirtualSubPop(sp);

			// real mating
			try {
				if (!m[idx]->mateSubPop(pop, sp, ind, ind + *itSize))
					return false;
			} catch (...) {
				cerr << "Mating in subpopulation " + toStr(sp) + " failed" << endl;
				throw;
			}
			ind += *itSize;
		}
		DBG_ASSERT(ind == scratch.rawIndEnd(sp), SystemError,
			"Mating scheme somehow does not fill the whole offspring population.");
		// we do not deactivate each time to save some time
		if (pop.hasActivatedVirtualSubPop(sp))
			pop.deactivateVirtualSubPop(sp);
		// if more than two mating schemes working on the same subpopulation,
		// it is better to shuffle offspring afterwards,
		if (m.size() > 1 && m_shuffleOffspring) {
			DBG_DO(DBG_MATING, cerr << "Random shuffle individuals in the offspring generation." << endl);
			getRNG().randomShuffle(scratch.rawIndBegin(sp), scratch.rawIndEnd(sp));
			scratch.setIndOrdered(false);
		}
	}                         // each subpopulation.
	submitScratch(pop, scratch);
	return true;
}


}


