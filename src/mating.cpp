/**
 *  $File: mating.cpp $
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

#include "mating.h"

#if PY_VERSION_HEX >= 0x03000000
#define PyInt_Check(x) PyLong_Check(x)
#define PyInt_AsLong(x) PyLong_AsLong(x)
#define PyInt_FromLong(x) PyLong_FromLong(x)
#endif

#if TR1_SUPPORT == 0
#include <map>
typedef std::map<ULONG, simuPOP::Individual *> IdMap;
#elif TR1_SUPPORT == 1
#include <unordered_map>
typedef std::unordered_map<size_t, simuPOP::Individual *> IdMap;
#else
#include <tr1/unordered_map>
typedef std::tr1::unordered_map<size_t, simuPOP::Individual *> IdMap;
#endif

namespace simuPOP
{

SeqSexModel::SeqSexModel(const vectorf &sex) : m_sex()
{
	DBG_FAILIF(sex.empty(), ValueError, "A sequence of sex is needed.");
	vectorf::const_iterator it = sex.begin() + 1;
	vectorf::const_iterator it_end = sex.end();
	for (; it != it_end; ++it)
	{
		int s = static_cast<int>(*it);
		m_sex.push_back(static_cast<Sex>(s));
	}
}

GlobalSeqSexModel::GlobalSeqSexModel(const vectorf &sex) : m_sex(), m_index(0)
{
	DBG_FAILIF(sex.empty(), ValueError, "A sequence of sex is needed.");
	vectorf::const_iterator it = sex.begin() + 1;
	vectorf::const_iterator it_end = sex.end();
	for (; it != it_end; ++it)
	{
		int s = static_cast<int>(*it);
		m_sex.push_back(static_cast<Sex>(s));
	}
}

Sex FuncSexModel::getSex(UINT count)
{
	if (m_generator.isValid())
	{
		long val;
		PyObject *obj = m_generator.next();
		PyObj_As_Int(obj, val);
		Py_DECREF(obj);
		return static_cast<Sex>(val);
	}
	else
	{
		PyObject *obj = NULL;
		try
		{
			obj = m_func("()");
			long val;
			PyObj_As_Int(obj, val);
			Py_DECREF(obj);
			return static_cast<Sex>(val);
		}
		catch (ValueError &)
		{
			if (PyGen_Check(obj))
			{
				PyErr_Clear();
				m_generator.set(obj);
				return getSex(count);
			}
			else
				throw ValueError("Passed function should be a function that return MALE/FEMALE or a generator.");
		}
		// this should not be reached.
		return MALE;
	}
	return MALE;
}

UINT FuncNumOffModel::getNumOff(ssize_t gen)
{
	if (m_generator.isValid())
	{
		int attempts = 0;
		long numOff = 0;
		while (++attempts < 50)
		{
			PyObject *obj = m_generator.next();
			PyObj_As_Int(obj, numOff);
			Py_DECREF(obj);
			DBG_DO(DBG_DEVEL, cerr << "Number of offspring produced from a generator: " << numOff << endl);
			if (numOff > 0)
				return numOff;
		}
		DBG_WARNIF(true, "One offspring is returned because user provided function returns 0 (#offspring) for more than 50 times.");
		return 1;
	}
	else
	{
		DBG_FAILIF(m_func.numArgs() > 1 || (m_func.numArgs() == 1 && m_func.arg(0) != "gen"),
				   ValueError, "Function passed to parameter numOffspring should have no parameter or a parameter named gen");
		int attempts = 0;
		long numOff = 0;
		while (++attempts < 50)
		{
			PyObject *obj = NULL;
			try
			{
				if (m_func.numArgs() == 0)
					obj = m_func("()");
				else
					obj = m_func("(i)", gen);
				PyObj_As_Int(obj, numOff);
				Py_DECREF(obj);
				if (numOff > 0)
					return numOff;
			}
			catch (ValueError &)
			{
				if (PyGen_Check(obj))
				{
					PyErr_Clear();
					m_generator.set(obj);
					return getNumOff(gen);
				}
				else
				{
					throw ValueError("Function should return a number or a generator.");
				}
			}
		}
		DBG_WARNIF(true, "One offspring is returned because user provided function returns 0 (#offspring) for more than 50 times.");
		return 1;
	}
	return 1;
}

OffspringGenerator::OffspringGenerator(const opList &ops,
									   const floatListFunc &numOffspring, const floatListFunc &sexMode) : m_numOffModel(NULL), m_sexModel(NULL), m_transmitters(ops), m_initialized(false)
{

	if (numOffspring.size() == 0)
	{
		DBG_FAILIF(!numOffspring.func().isValid(), ValueError,
				   "Please specify either number of offspring or a function.");
		m_numOffModel = new FuncNumOffModel(numOffspring.func());
	}
	else if (numOffspring.size() == 1)
	{
		DBG_FAILIF(static_cast<int>(numOffspring[0]) <= 0, ValueError,
				   "Number of offspring has to be positive.");
		m_numOffModel = new ConstNumOffModel(static_cast<UINT>(numOffspring[0]));
	}
	else if (numOffspring.size() > 1)
	{
		int mode = static_cast<int>(numOffspring[0]);

		if (mode == GEOMETRIC_DISTRIBUTION)
		{
			DBG_FAILIF(numOffspring.size() < 2 || fcmp_lt(numOffspring[1], 0) || fcmp_gt(numOffspring[1], 1.),
					   ValueError, "P for a geometric distribution should be within [0,1]");
			m_numOffModel = new GeometricNumOffModel(numOffspring[1]);
		}
		else if (mode == POISSON_DISTRIBUTION)
		{
			m_numOffModel = new PoissonNumOffModel(numOffspring[1]);
		}
		else if (mode == BINOMIAL_DISTRIBUTION)
		{
			DBG_FAILIF(numOffspring.size() < 3 || static_cast<UINT>(numOffspring[2]) <= 1.,
					   ValueError, "If mode is BINOMIAL_DISTRIBUTION, the second parameter should be greater than 1");
			DBG_FAILIF(numOffspring.size() < 2 || fcmp_le(numOffspring[1], 0) || fcmp_gt(numOffspring[1], 1.),
					   ValueError, "P for a Bionomial distribution should be within (0,1].");
			DBG_FAILIF(numOffspring.size() < 3 || numOffspring[2] < 1,
					   ValueError, "Max number of offspring in a binomial distribution should be greater than or equal to 1.");
			m_numOffModel = new BinomialNumOffModel(static_cast<UINT>(numOffspring[2]), numOffspring[1]);
		}
		else if (mode == UNIFORM_DISTRIBUTION)
		{
			DBG_FAILIF(numOffspring.size() < 3 || numOffspring[2] < static_cast<UINT>(numOffspring[1]),
					   ValueError, "If mode is UNIFORM_DISTRIBUTION, numOffspringParam should be greater than numOffspring");
			m_numOffModel = new UniformNumOffModel(static_cast<UINT>(numOffspring[1]),
												   static_cast<UINT>(numOffspring[2]));
		}
		else
		{
			throw ValueError("Wrong mating numoffspring mode. Should be one of \n"
							 "NumOffspring, NumOffspringEachFamily and GEometricDistribution");
		}
	}
	if (sexMode.size() == 0)
	{
		DBG_FAILIF(!sexMode.func().isValid(), ValueError,
				   "Please specify one of the sex modes or a Python function.");
		m_sexModel = new FuncSexModel(sexMode.func());
	}
	else
	{
		int mode = static_cast<int>(sexMode[0]);

		if (mode == NO_SEX)
		{
			DBG_FAILIF(sexMode.size() != 1, ValueError, "No parameter is allowed for NO_SEX mode");
			m_sexModel = new NoSexModel();
		}
		else if (mode == RANDOM_SEX)
		{
			DBG_FAILIF(sexMode.size() != 1, ValueError, "No parameter is allowed for RANDOM_SEX mode");
			m_sexModel = new RandomSexModel();
		}
		else if (mode == PROB_OF_MALES)
		{
			DBG_FAILIF(sexMode.size() != 2, ValueError, "A parameter is required for mode PROB_OF_MALES.");
			DBG_FAILIF(fcmp_lt(sexMode[1], 0) || fcmp_gt(sexMode[1], 1),
					   ValueError, "Probability of male has to be between 0 and 1");
			m_sexModel = new ProbOfMalesSexModel(sexMode[1]);
		}
		else if (mode == NUM_OF_MALES)
		{
			DBG_FAILIF(sexMode.size() != 2, ValueError, "A parameter is required for mode NUM_OF_MALES.");
			m_sexModel = new NumOfMalesSexModel(static_cast<UINT>(sexMode[1]));
		}
		else if (mode == NUM_OF_FEMALES)
		{
			DBG_FAILIF(sexMode.size() != 2, ValueError, "A parameter is required for mode NUM_OF_FEMALES.");
			m_sexModel = new NumOfFemalesSexModel(static_cast<UINT>(sexMode[1]));
		}
		else if (mode == SEQUENCE_OF_SEX)
		{
			DBG_FAILIF(sexMode.size() <= 2, ValueError, "A sequence of sex is required for mode SEQUENCE_OF_SEX.");
			m_sexModel = new SeqSexModel(sexMode.elems());
		}
		else if (mode == GLOBAL_SEQUENCE_OF_SEX)
		{
			DBG_FAILIF(sexMode.size() <= 2, ValueError, "A sequence of sex is required for mode GLOBAL_SEQUENCE_OF_SEX.");
			m_sexModel = new GlobalSeqSexModel(sexMode.elems());
		}
		else
		{
			DBG_FAILIF(true, ValueError, "Unrecognized sexMode.");
		}
	}
}

UINT OffspringGenerator::numOffspring(ssize_t gen)
{
	return m_numOffModel->getNumOff(gen);
}

bool OffspringGenerator::parallelizable() const
{
#ifdef MUTANTALLELE
	// mutant allele model cannot generate offspring in parallele
	// due to the use of std::map
	return false;
#else
	if (!m_sexModel->parallelizable())
		return false;
	if (!m_numOffModel->parallelizable())
		return false;
	opList::const_iterator iop = m_transmitters.begin();
	opList::const_iterator iopEnd = m_transmitters.end();
	for (; iop != iopEnd; ++iop)
	{
		if (!(*iop)->parallelizable())
			return false;
	}
	return true;
#endif
}

Sex OffspringGenerator::getSex(UINT count)
{
	return m_sexModel->getSex(count);
}

void OffspringGenerator::initialize(const Population &pop, size_t /* subPop */)
{
	opList::const_iterator iop = m_transmitters.begin();
	opList::const_iterator iopEnd = m_transmitters.end();

	for (; iop != iopEnd; ++iop)
		(*iop)->initializeIfNeeded(*pop.rawIndBegin());

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

UINT OffspringGenerator::generateOffspring(Population &pop, Population &offPop, Individual *dad, Individual *mom,
										   RawIndIterator &it,
										   RawIndIterator &itEnd)
{
	DBG_ASSERT(initialized(), ValueError,
			   "Offspring generator is not initialized before used to generate offspring");

	// generate numOff offspring per mating, or until it  reaches offEnd
	UINT count = 0;
	bool accept = true;
	UINT numOff = numOffspring(pop.gen());
	UINT attempt = 0;
	while (attempt < numOff && it != itEnd)
	{
		// not all families have the same size because some offspring
		// may be discarded (count).
		++attempt;

		// set sex, during mating operator will try to
		// follow the offspring sex (e.g. pass X or Y chromosome)
		it->setSex(getSex(count));
		// set first offspring
		it->setFirstOffspring(count == 0);
		//
		accept = true;
		opList::const_iterator iop = m_transmitters.begin();
		opList::const_iterator iopEnd = m_transmitters.end();
		for (; iop != iopEnd; ++iop)
		{
			if (!(*iop)->isActive(pop.rep(), pop.gen()))
				continue;
			if (!(*iop)->applyDuringMating(pop, offPop, it, dad, mom))
			{
				accept = false;
				break;
			}
		}

		if (accept)
		{
			++it;
			++count;
		}
	}
	return count;
}

ControlledOffspringGenerator::ControlledOffspringGenerator(
	const lociList &loci, const uintList &alleles, PyObject *freqFunc,
	const opList &ops, const floatListFunc &numOffspring,
	const floatListFunc &sexMode)
	: OffspringGenerator(ops, numOffspring, sexMode),
	  m_loci(loci), m_alleles(alleles.elems()), m_freqFunc(freqFunc),
	  m_expAlleles(), m_totAllele(), m_curAllele()

{
	if (!m_freqFunc.isValid())
		throw ValueError("Please specify a valid frequency function");
}

ControlledOffspringGenerator::ControlledOffspringGenerator(const ControlledOffspringGenerator &rhs)
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

void ControlledOffspringGenerator::getExpectedAlleles(const Population &pop,
													  vectorf &expFreq)
{
	// determine expected number of alleles of each allele
	// at each subpopulation.
	const vectoru &loci = m_loci.elems(&pop);
	size_t nLoci = loci.size();
	size_t numSP = pop.numSubPop();

	DBG_ASSERT(expFreq.size() == nLoci || expFreq.size() == nLoci * numSP, SystemError,
			   "Expect expected frequency for all loci (or all loci in all subpopulation)");

	for (size_t i = 0; i < expFreq.size(); ++i)
	{
		DBG_FAILIF(expFreq[i] < 0 || expFreq[i] > 1., ValueError,
				   "Expected frequency out of range");
	}

	// in the order of Loc1: sp1, sp2, sp3, sloc2: p1, sp2, sp3
	m_expAlleles.resize(nLoci * numSP);
	if (numSP > 1 && expFreq.size() == nLoci)
	{
		// if exp frequencies in subpopulation is not specified.
		// I need to find the current allele frequencies
		// and use them as proportions for the next generation.
		for (size_t i = 0; i < nLoci; ++i)
		{
			size_t locus = loci[i];
			Allele allele = TO_ALLELE(m_alleles[i]);

			// determine the number alleles at each subpopulation.
			vectorf curFreq(numSP, 0);
			bool hasAllele = false;
			for (size_t sp = 0, n = 0; sp < numSP; ++sp)
			{
				IndAlleleIterator a = const_cast<Population &>(pop).alleleIterator(locus, sp);
				for (; a.valid(); ++a)
					if (ALLELE_AS_UNSINGED(DEREF_ALLELE(a)) == allele)
						n++;
				hasAllele = hasAllele || n > 0;
				curFreq[sp] = double(n) / (pop.subPopSize(sp) * pop.ploidy());
			}

			DBG_DO(DBG_MATING, cerr << "Current frequency at locus " << locus
									<< " is " << curFreq << endl);

			// if there is no alleles
			if (!hasAllele && expFreq[i] > 0.)
				throw RuntimeError((boost::format("No disease allele exists at generation %1%"
												  ", but expected allele frequency at locus %2% is greater than 0.") %
									pop.gen() % locus)
									   .str());

			DBG_WARNIF(hasAllele && fcmp_eq(expFreq[i], 0.), (boost::format("Disease allele exists at generation %1%"
																			", but expected allele frequency is zero.") %
															  pop.gen())
																 .str());

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
	}
	else
	{
		// simpler case, one subpopulation, or with gieven allele frequency
		for (size_t i = 0; i < nLoci; ++i)
		{
			for (size_t sp = 0; sp < numSP; ++sp)
			{
#ifndef OPTIMIZED
				size_t locus = loci[i];
				Allele allele = TO_ALLELE(m_alleles[i]);
				ULONG n = 0;
				// go through all alleles
				IndAlleleIterator a = const_cast<Population &>(pop).alleleIterator(locus, sp);
				for (; a.valid(); ++a)
				{
					if (ALLELE_AS_UNSINGED(DEREF_ALLELE(a)) == allele)
						n++;
				}

				// if there is no alleles
				if (n == 0 && expFreq[sp * nLoci + i] > 0.)
					throw RuntimeError((boost::format("No disease allele exists at generation %1% but exp allele frequency at locus %2%"
													  " in subpopulation %3% is greater than 0.") %
										pop.gen() % locus % sp)
										   .str());
#endif
				m_expAlleles[numSP * i + sp] = static_cast<UINT>(pop.subPopSize(sp) * pop.ploidy() * expFreq[sp * nLoci + i]);
				if (expFreq[sp * nLoci + i] > 0. && m_expAlleles[numSP * i + sp] == 0)
					m_expAlleles[numSP * i + sp] = 1;
			}
		}
	}
}

void ControlledOffspringGenerator::initialize(const Population &pop, size_t subPop)
{
	OffspringGenerator::initialize(pop, subPop);

	// expected frequency at each locus
	if (subPop == 0)
	{
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
	for (size_t i = 0; i < m_loci.size(); ++i)
	{
		size_t maxCount = pop.subPopSize(subPop) * pop.ploidy();
		m_totAllele[i] = m_expAlleles[subPop + pop.numSubPop() * i];
		if (m_totAllele[i] > maxCount)
		{
			cerr << "Warning: number of planned affected alleles exceed population size.";
			m_totAllele[i] = maxCount;
		}
		if (2 * m_totAllele[i] > maxCount)
		{
			m_flip[i] = true;
			m_totAllele[i] = maxCount - m_totAllele[i];
		}
	}
	//

	// it is possible that no disease allele is required
	// so everyone is accepted
	m_freqRequMet = true;
	//
	for (size_t i = 0; i < m_loci.size(); ++i)
	{
		if (m_totAllele[i] > 0)
		{
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
UINT ControlledOffspringGenerator::generateOffspring(Population &pop, Population &offPop, Individual *dad, Individual *mom,
													 RawIndIterator &offBegin,
													 RawIndIterator &offEnd)
{
	const vectoru &loci = m_loci.elems(&pop);
	size_t nLoci = loci.size();
	//
	// generate m_numOffspring offspring per mating
	// record family size (this may be wrong for the last family)
	//
	RawIndIterator itBegin = offBegin;
	UINT numOff = OffspringGenerator::generateOffspring(pop, offPop,
														dad, mom, offBegin, offEnd);

	//
	if (numOff == 0)
		return numOff;
	//
	// count number of alleles in the family.
	vectori na(nLoci, 0);
	bool hasAff = false;
	size_t totNumLoci = pop.totNumLoci();
	// we know that scratch population has ordered linear genotype
	for (size_t i = 0; i < nLoci; ++i)
	{
		GenoIterator ptr = itBegin->genoBegin();
		for (size_t j = 0; j < numOff * pop.ploidy(); ++j, ptr += totNumLoci)
		{
			if (m_flip[i] ? (DEREF_ALLELE(ptr + loci[i]) != TO_ALLELE(m_alleles[i]))
						  : (DEREF_ALLELE(ptr + loci[i]) == TO_ALLELE(m_alleles[i])))
			{
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
	if (m_freqRequMet)
	{
		if (!hasAff)
		{
			// has AA, so no need to compromise
			m_AAattempt = 10000;
			accept = true;
		}
		// tried 200 times, no AA is found.
		else if (m_AAattempt == 0)
		{
			m_AAattempt = 200;
			accept = true;
		}
		m_AAattempt--;
	}
	else
	{ // do not use stack
		if (hasAff)
		{
			for (size_t i = 0; i < nLoci; ++i)
			{
				// accept the whole family, if we need this allele
				if (m_curAllele[i] < m_totAllele[i] && na[i] > 0)
				{
					accept = true;
					break;
				}
			}
			m_aaAttempt = 10000;
		}
		else if (m_aaAttempt == 0)
		{
			m_aaAttempt = 200;
			accept = true;
		}
		m_aaAttempt--;
	}
	//
	// reject this family
	if (!accept)
	{
		// it relocate to its begin point
		// DBG_DO(DBG_MATING, cerr << "Reject " << na << endl);
		offBegin = itBegin;
		return 0;
	}
	// accpet this family, see if all done.
	for (size_t i = 0; i < nLoci; ++i)
		m_curAllele[i] += na[i];
	if (!m_freqRequMet)
	{
		m_freqRequMet = true;
		for (size_t i = 0; i < nLoci; ++i)
		{
			if (m_curAllele[i] < m_totAllele[i])
				m_freqRequMet = false;
		}
	}
	return numOff;
}

void SequentialParentChooser::initialize(Population &pop, size_t sp)
{
	//
	m_curInd = 0;
	m_index.clear();
	if (m_choice == ANY_SEX)
	{
		m_begin = pop.indIterator(sp);
		m_ind = m_begin;
	}
	else
	{
		IndIterator it = pop.indIterator(sp);
		Sex s = m_choice == MALE_ONLY ? MALE : FEMALE;
		for (; it.valid(); ++it)
		{
			if (it->sex() == s)
				m_index.push_back(it.rawIter());
		}
		DBG_FAILIF(m_index.empty(), RuntimeError,
				   string("No ") + (s == MALE ? "male" : "female") + " individual exists in a population.");
	}
	m_initialized = true;
}

ParentChooser::IndividualPair SequentialParentChooser::chooseParents()
{
	DBG_ASSERT(initialized(), SystemError,
			   "Please initialize this parent chooser before using it");
	if (m_choice == ANY_SEX)
	{
		if (!m_ind.valid())
		{
			m_ind = m_begin;
			DBG_ASSERT(m_ind.valid(), RuntimeError, "No valid individual if found.")
		}
		return ParentChooser::IndividualPair(&*m_ind++, (Individual *)(0));
	}
	else
	{
		if (m_curInd == m_index.size())
			m_curInd = 0;

		return std::make_pair(&*m_index[m_curInd++], static_cast<Individual *>(0));
	}
}

void RandomParentChooser::initialize(Population &pop, size_t sp)
{
	m_basePtr = pop.rawIndBegin();

	m_index.clear();

	m_selection = m_replacement && pop.hasInfoField(m_selectionField);
	size_t fit_id = m_selection ? pop.infoIdx(m_selectionField) : 0;
	// In a virtual subpopulation, because m_begin + ... is **really** slow
	// It is a good idea to cache IndIterators. This is however inefficient
	// for non-virtual populations
	vectorf fitness;
	if (pop.hasActivatedVirtualSubPop(sp) || m_choice != ANY_SEX)
	{
		IndIterator it = pop.indIterator(sp);
		if (m_choice == ANY_SEX)
		{
			for (; it.valid(); ++it)
			{
				m_index.push_back(it.rawIter());
				if (m_selection)
					fitness.push_back(it->info(fit_id));
			}
			DBG_FAILIF(m_index.empty(), RuntimeError, "Can not select parent from an empty subpopulation.");
		}
		else
		{
			Sex s = m_choice == MALE_ONLY ? MALE : FEMALE;
			for (; it.valid(); ++it)
			{
				if (it->sex() == s)
				{
					m_index.push_back(it.rawIter());
					if (m_selection)
						fitness.push_back(it->info(fit_id));
				}
			}
			DBG_FAILIF(m_index.empty(), RuntimeError,
					   string("No ") + (s == MALE ? "male" : "female") + " individual exists in a population.");
		}
	}
	else
	{
		if (!m_replacement)
			for (IndIterator it = pop.indIterator(sp); it.valid(); ++it)
				m_index.push_back(it.rawIter());
		if (m_selection) {
			IndInfoIterator it = pop.infoBegin(fit_id, sp);
			IndInfoIterator it_end = pop.infoEnd(fit_id, sp);
			for (; it != it_end; ++it) {
				fitness.push_back(*it);
			}
		}
	}

	if (m_selection)
		m_sampler.set(fitness.begin(), fitness.end());
	else
	{
		m_size = m_index.size();
		if (m_size == 0) // if m_index is not used (no VSP)
			m_size = pop.subPopSize(sp);
	}

	if (!m_replacement)
		getRNG().randomShuffle(m_index.begin(), m_index.end());

	m_shift = pop.subPopBegin(sp);
	m_initialized = true;
}

ParentChooser::IndividualPair RandomParentChooser::chooseParents()
{
	DBG_ASSERT(initialized(), SystemError,
			   "Please initialize this parent chooser before using it");
	// choose a parent
	if (!m_replacement)
	{
		if (m_index.empty())
			throw RuntimeError("All parents have been chosen.");
		Individual *ind = &*m_index.back();
		m_index.pop_back();
		return IndividualPair(ind, (Individual *)(0));
	}
	Individual *ind = NULL;
	if (m_index.empty())
	{
		if (m_selection)
			// basePtr points to the beginning of the population, not subpopulation
			ind = &*(m_basePtr + m_shift + m_sampler.draw());
		else
			ind = &*(m_basePtr + m_shift + getRNG().randInt(static_cast<ULONG>(m_size)));
	}
	else
	{
		if (m_selection)
			ind = &*(m_index[m_sampler.draw()]);
		else
			ind = &*(m_index[getRNG().randInt(static_cast<ULONG>(m_size))]);
	}
	return IndividualPair(ind, (Individual *)(0));
}

void RandomParentsChooser::initialize(Population &pop, size_t subPop)
{
	m_numMale = 0;
	m_numFemale = 0;

	m_maleIndex = 0;
	m_femaleIndex = 0;

	// index might not all be used because of virtual subpopulation
	m_index.resize(pop.subPopSize(subPop));

	m_selection = m_replacement && pop.hasInfoField(m_selectionField);
	size_t fit_id = 0;
	if (m_selection)
	{
		fit_id = pop.infoIdx(m_selectionField);
		m_fitness.resize(pop.subPopSize(subPop));
	}

	vector<RawIndIterator>::iterator maleIndex = m_index.begin();
	vector<RawIndIterator>::reverse_iterator femaleIndex = m_index.rbegin();
	vectorf::iterator maleFitness;
	vectorf::reverse_iterator femaleFitness;
	if (m_selection)
	{
		maleFitness = m_fitness.begin();
		femaleFitness = m_fitness.rbegin();
	}
	IndIterator it = pop.indIterator(subPop);
	for (; it.valid(); it++)
	{
		if (it->sex() == MALE)
		{
			*maleIndex++ = it.rawIter();
			if (m_selection)
				*maleFitness++ = it->info(fit_id);
		}
		else
		{
			*femaleIndex++ = it.rawIter();
			if (m_selection)
				*femaleFitness++ = it->info(fit_id);
		}
	}
	// m_numMale + m_numFemale might not be pop.subPopSize because of virtual subpopulation
	m_numMale = maleIndex - m_index.begin();
	m_numFemale = femaleIndex - m_index.rbegin();

	if (!m_replacement)
	{
		DBG_FAILIF(m_numMale == 0, IndexError, "No male individual in this population");
		DBG_FAILIF(m_numFemale == 0, IndexError, "No female individual in this population");
		getRNG().randomShuffle(m_index.begin(), m_index.begin() + m_numMale);
		getRNG().randomShuffle(m_index.rbegin(), m_index.rbegin() + m_numFemale);
	}

	if (m_selection)
	{
		m_malesampler.set(m_fitness.begin(), m_fitness.begin() + m_numMale);
		m_femalesampler.set(m_fitness.rbegin(), m_fitness.rbegin() + m_numFemale);
		DBG_DO(DBG_DEVEL, cerr << "Male and Female fitness " << m_fitness << endl);
	}

	DBG_FAILIF(!m_replacement && m_selection, ValueError,
			   "Selection is not allowed in random sample without replacement");

	m_initialized = true;
}

ParentChooser::IndividualPair RandomParentsChooser::chooseParents()
{
	DBG_ASSERT(initialized(), SystemError,
			   "Please initialize this parent chooser before using it");

	Individual *dad = NULL;
	Individual *mom = NULL;

	if (!m_replacement)
	{
		if (m_femaleIndex >= m_numFemale)
			throw RuntimeError("All females have been chosen.");
		mom = &**(m_index.rbegin() + m_femaleIndex++);

		if (m_maleIndex >= m_numMale)
			throw RuntimeError("All males have been chosen.");
		dad = &**(m_index.begin() + m_maleIndex++);
		return std::make_pair(dad, mom);
	}

	// this exceptio should be raised also in optimized mode because the cause
	// can be random.
	if (m_numMale == 0)
		throw RuntimeError("RandomParentsChooser fails because there is no male individual in a subpopulation.");
	if (m_numFemale == 0)
		throw RuntimeError("RandomParentsChooser fails because there is no female individual in a subpopulation ");

	if (m_selection)
	{
		// using weighted sampler.
		dad = &**(m_index.begin() + m_malesampler.draw());
		mom = &**(m_index.rbegin() + m_femalesampler.draw());
	}
	else
	{
		dad = &**(m_index.begin() + getRNG().randInt(static_cast<ULONG>(m_numMale)));
		mom = &**(m_index.rbegin() + getRNG().randInt(static_cast<ULONG>(m_numFemale)));
	}
	return std::make_pair(dad, mom);
}

void PolyParentsChooser::initialize(Population &pop, size_t subPop)
{
	m_numMale = 0;
	m_numFemale = 0;
#ifdef _OPENMP
	for (size_t i = 0; i < m_polyCount.size(); i++)
		m_polyCount[i] = 0;
#else
	m_polyCount = 0;
#endif

	IndIterator it = pop.indIterator(subPop);
	for (; it.valid(); ++it)
	{
		if (it->sex() == MALE)
			m_numMale++;
		else
			m_numFemale++;
	}

	// allocate memory at first for performance reasons
	m_maleIndex.resize(m_numMale);
	m_femaleIndex.resize(m_numFemale);

	m_selection = pop.hasInfoField(m_selectionField);

	size_t fit_id = 0;
	if (m_selection)
	{
		fit_id = pop.infoIdx(m_selectionField);
		m_maleFitness.resize(m_numMale);
		m_femaleFitness.resize(m_numFemale);
	}

	m_numMale = 0;
	m_numFemale = 0;

	it = pop.indIterator(subPop);
	for (; it.valid(); it++)
	{
		if (it->sex() == MALE)
		{
			m_maleIndex[m_numMale] = it.rawIter();
			if (m_selection)
				m_maleFitness[m_numMale] = it->info(fit_id);
			m_numMale++;
		}
		else
		{
			m_femaleIndex[m_numFemale] = it.rawIter();
			if (m_selection)
				m_femaleFitness[m_numFemale] = it->info(fit_id);
			m_numFemale++;
		}
	}

	if (m_selection)
	{
		m_malesampler.set(m_maleFitness.begin(), m_maleFitness.end());
		m_femalesampler.set(m_femaleFitness.begin(), m_femaleFitness.end());
		DBG_DO(DBG_DEVEL, cerr << "Male fitness " << m_maleFitness << endl);
		DBG_DO(DBG_DEVEL, cerr << "Female fitness " << m_femaleFitness << endl);
	}

	m_initialized = true;
}

ParentChooser::IndividualPair PolyParentsChooser::chooseParents()
{
	DBG_ASSERT(initialized(), SystemError,
			   "Please initialize this parent chooser before using it");

	Individual *dad = NULL;
	Individual *mom = NULL;
#ifdef _OPENMP
	int threadID = omp_get_thread_num();
#endif
#ifdef _OPENMP
	if (m_polyNum > 1 && m_polyCount[threadID] > 0)
	{
		if (m_polySex == MALE)
			dad = m_lastParent[threadID];
		else
			mom = m_lastParent[threadID];
		m_polyCount[threadID]--;
	}
#else
	if (m_polyNum > 1 && m_polyCount > 0)
	{
		if (m_polySex == MALE)
			dad = m_lastParent;
		else
			mom = m_lastParent;
		m_polyCount--;
	}
#endif

	// using weidhted sampler.
	if (dad == NULL)
	{
		if (m_numMale == 0)
			throw RuntimeError("PolyParentsChooser fails because there is no male individual in a subpopulation.");

		if (m_selection)
			dad = &*(m_maleIndex[m_malesampler.draw()]);
		else
			dad = &*(m_maleIndex[getRNG().randInt(static_cast<ULONG>(m_numMale))]);

		if (m_polySex == MALE && m_polyNum > 1)
		{
#ifdef _OPENMP
			m_polyCount[threadID] = m_polyNum - 1;
			m_lastParent[threadID] = dad;
#else
			m_polyCount = m_polyNum - 1;
			m_lastParent = dad;
#endif
		}
	}

	if (mom == NULL)
	{
		if (m_numFemale == 0)
			throw RuntimeError("PolyParentsChooser fails because there is no female individual in a subpopulation.");

		if (m_selection)
			mom = &*(m_femaleIndex[m_femalesampler.draw()]);
		else
			mom = &*(m_femaleIndex[getRNG().randInt(static_cast<ULONG>(m_numFemale))]);

		if (m_polySex == FEMALE && m_polyNum > 1)
		{
#ifdef _OPENMP
			m_polyCount[threadID] = m_polyNum - 1;
			m_lastParent[threadID] = mom;
#else
			m_polyCount = m_polyNum - 1;
			m_lastParent = mom;
#endif
		}
	}
	return std::make_pair(dad, mom);
}

/*
   void infoParentsChooser::initialize(Population & pop, size_t sp)
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
   ParentChooser::IndividualPair infoParentsChooser::chooseParents()
   {
    DBG_ASSERT(initialized(), SystemError,
        "Please initialize this parent chooser before using it");
    Individual * par1 = RandomParentChooser::chooseParents().first;
    Sex sex1 = par1->sex();
    // there is no valid information field value
    if (m_degenerate) {
        int attempt = 0;
        while (++attempt < 1000) {
            Individual * par2 = RandomParentChooser::chooseParents().first;
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

CombinedParentsChooser::CombinedParentsChooser(const ParentChooser &fatherChooser,
											   const ParentChooser &motherChooser, bool allowSelfing)
	: m_fatherChooser(fatherChooser.clone()),
	  m_motherChooser(motherChooser.clone()),
	  m_allowSelfing(allowSelfing)
{
}

void CombinedParentsChooser::initialize(Population &pop, size_t sp)
{
	m_fatherChooser->initialize(pop, sp);
	m_motherChooser->initialize(pop, sp);
}

ParentChooser::IndividualPair CombinedParentsChooser::chooseParents()
{
	size_t attempts = 0;
	while (attempts < 100)
	{
		ParentChooser::IndividualPair p1 = m_fatherChooser->chooseParents();
		ParentChooser::IndividualPair p2 = m_motherChooser->chooseParents();
		Individual *dad = p1.first != NULL ? p1.first : p1.second;
		Individual *mom = p2.second != NULL ? p2.second : p2.first;

		if (dad != mom || m_allowSelfing)
			return ParentChooser::IndividualPair(dad, mom);
		if (attempts++ == 100)
			throw RuntimeError("Failed to select distinct parents using CombinedParentsChooser.");
	}
	// just make the compiler happy
	return ParentChooser::IndividualPair((Individual *)(0), (Individual *)(0));
}

void CombinedParentsChooser::finalize()
{
	m_fatherChooser->finalize();
	m_motherChooser->finalize();
}

PyParentsChooser::PyParentsChooser(PyObject *pc)
	: ParentChooser(), m_func(pc), m_popObj(NULL),
	  m_generator(NULL)
{
}

void PyParentsChooser::initialize(Population &pop, size_t sp)
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

	PyObject *args = PyTuple_New(m_func.numArgs());
	DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");
	for (size_t i = 0; i < m_func.numArgs(); ++i)
	{
		const string &arg = m_func.arg(i);
		if (arg == "pop")
			PyTuple_SET_ITEM(args, i, m_popObj);
		else if (arg == "subPop")
			PyTuple_SET_ITEM(args, i, PyInt_FromLong(static_cast<long>(sp)));
		else
		{
			DBG_FAILIF(true, ValueError,
					   "Only parameters 'pop' and 'subPop' are acceptable in a generator function.");
		}
	}
	m_generator.set(m_func(args));
	Py_DECREF(args);
	m_initialized = true;
}

ParentChooser::IndividualPair PyParentsChooser::chooseParents()
{
	DBG_ASSERT(initialized(), SystemError,
			   "Please initialize this parent chooser before using it");

	PyObject *item = m_generator.next();

#ifndef OPTIMIZED
	if (item == NULL && debug(DBG_GENERAL))
	{
		PyErr_Print();
		PyErr_Clear();
	}
	DBG_FAILIF(item == NULL, ValueError,
			   "User-defined function yield invalid value.");
#endif

	if (PyInt_Check(item) || PyLong_Check(item))
	{
		long parent;
		PyObj_As_Int(item, parent);
#ifndef OPTIMIZED
		DBG_ASSERT(static_cast<unsigned>(parent) < m_size,
				   ValueError, (boost::format("Returned index (%1%) is greater than subpopulation size %2%") % parent % m_size).str());
#endif
		Py_DECREF(item);
		return ParentChooser::IndividualPair(&*(m_begin + parent), (Individual *)(0));
	}
	else if (PySequence_Check(item))
	{
		DBG_ASSERT(PySequence_Size(item) == 2, RuntimeError,
				   "Parents should be returned in the form of a sequence of two elements");

		Individual *parents[2];
		for (size_t i = 0; i < 2; ++i)
		{
			PyObject *v = PySequence_GetItem(item, i);
			if (PyInt_Check(v) || PyLong_Check(v))
			{
				ULONG idx = PyInt_AsLong(v);
				DBG_ASSERT(idx < m_size, ValueError, (boost::format("Returned parent index (%1%) is greater than subpopulation size %2%") % idx % m_size).str());
				parents[i] = &*(m_begin + idx);
			}
			else
			{
				void *ind = pyIndPointer(v);
				DBG_ASSERT(ind, ValueError, "Invalid type of returned parent.");
				parents[i] = reinterpret_cast<Individual *>(ind);
			}
			Py_DECREF(v);
		}
		// PySequence_GetItem gets a new reference so we should DECREF
		Py_DECREF(item);
		return ParentChooser::IndividualPair(parents[0], parents[1]);
	}
	else
	{
		// is an individual object is returned?
		void *ind = pyIndPointer(item);
		DBG_ASSERT(ind, ValueError, "Invalid type of returned parent.");
		return ParentChooser::IndividualPair(reinterpret_cast<Individual *>(ind), (Individual *)(0));
	}
	// this should not be reached
	return ParentChooser::IndividualPair((Individual *)(0), (Individual *)(0));
}

void PyParentsChooser::finalize()
{
	DBG_FAILIF(m_popObj == NULL, SystemError, "Python generator is not properly initialized.");
	Py_DECREF(m_popObj);
	m_generator.set(NULL);
	m_popObj = NULL;
	m_initialized = false;
}

MatingScheme::MatingScheme(const uintListFunc &subPopSize)
	: m_subPopSize(subPopSize)
{
}

bool MatingScheme::prepareScratchPop(Population &pop, Population &scratch)
{
	if (scratch.genoStruIdx() != pop.genoStruIdx())
		scratch.fitGenoStru(pop.genoStruIdx());

	// use population structure of pop
	if (m_subPopSize.empty() && !m_subPopSize.func().isValid())
		scratch.fitSubPopStru(pop.subPopSizes(), pop.subPopNames());
	else if (!m_subPopSize.empty()) // set subPoplation size
		scratch.fitSubPopStru(m_subPopSize.elems(), pop.subPopNames());
	else
	{ // use m_subPopSizeFunc
		const pyFunc &func = m_subPopSize.func();
		PyObject *args = PyTuple_New(func.numArgs());
		DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

		for (size_t i = 0; i < func.numArgs(); ++i)
		{
			const string &arg = func.arg(i);
			if (arg == "gen")
				PyTuple_SET_ITEM(args, i, PyInt_FromLong(static_cast<long>(pop.gen())));
			else if (arg == "pop")
				PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(&pop)));
			else
			{
				DBG_FAILIF(true, ValueError,
						   "Only parameters 'gen' and 'pop' are acceptable in a demographic function.");
			}
		}
		vectori res = func(PyObj_As_IntArray, args);
		Py_XDECREF(args);

		if (res.empty())
		{
			DBG_DO(DBG_SIMULATOR, cerr << "Stop iteration due to empty offspring population size." << endl);
			return false;
		}
		vectoru sz(res.size());
		for (size_t i = 0; i < res.size(); i++)
		{
			if (res[i] < 0)
				throw ValueError((boost::format("Negative population size %1% returned for subpopulation %2%") % res[i] % i).str());
			sz[i] = static_cast<ULONG>(res[i]);
		}

		// allow change of pop size of scratch
		scratch.fitSubPopStru(sz, pop.subPopNames());
	}
	// this is not absolutely necessary but will reduce confusions
	scratch.setVirtualSplitter(pop.virtualSplitter());
	// the scratch population has the same generation and rep number as the parent population.
	// The numbers might be used by during mating operator
	scratch.setGen(pop.gen());
	scratch.setRep(pop.rep());
	scratch.clearInfo();
#ifdef MUTANTALLELE
	// for mutant allele, clearing all existing genotype will make subsequent
	// copyChromosomes much faster ...
	scratch.setGenotype(vectoru(1, 0));
#endif
	DBG_DO(DBG_SIMULATOR, cerr << "New subpop size " << scratch.subPopSizes() << endl);

	DBG_FAILIF(scratch.numSubPop() != pop.numSubPop(),
			   ValueError, (boost::format("number of subPopulaitons must agree.\n Pre: %1% now: %2%") % pop.numSubPop() % scratch.numSubPop()).str());
	return true;
}

bool MatingScheme::mate(Population &pop, Population &scratch)
{
	// scrtach will have the right structure.
	if (!prepareScratchPop(pop, scratch))
		return false;
	for (size_t sp = 0; sp < static_cast<size_t>(pop.numSubPop()); ++sp)
		if (!mateSubPop(pop, scratch, sp, scratch.rawIndBegin(sp), scratch.rawIndEnd(sp)))
			return false;
	submitScratch(pop, scratch);
	return true;
}

void MatingScheme::submitScratch(Population &pop, Population &scratch)
{
	// use scratch population,
	pop.push(scratch);
	scratch.validate("after push and discard");
}

HomoMating::HomoMating(ParentChooser &chooser,
					   OffspringGenerator &generator,
					   const uintListFunc &subPopSize,
					   subPopList subPops, double weight)
	: MatingScheme(subPopSize), m_subPops(subPops), m_weight(weight)
{
	m_ParentChooser = chooser.clone();
	m_OffspringGenerator = generator.clone();
}

string HomoMating::describe(bool format) const
{
	string desc = "<simuPOP.HomoMating> a homogeneous mating scheme that uses\n<ul>\n<li>" + m_ParentChooser->describe(false) + "\n<li>" + m_OffspringGenerator->describe(false) + "</ul>\n";

	return format ? formatDescription(desc) : desc;
}

bool HomoMating::mateSubPop(Population &pop, Population &offPop, size_t subPop,
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
	// If the parent chooser is not parallelizable, or if openMP is not supported
	// or if number of thread is set to 1, use the sequential method.
	if (!m_ParentChooser->parallelizable() || numThreads() == 1 || !m_OffspringGenerator->parallelizable())
	{
		DBG_DO(DBG_MATING, cerr << "Mating is done in single-thread mode" << endl);
		while (it != offEnd)
		{
			Individual *dad = NULL;
			Individual *mom = NULL;
			ParentChooser::IndividualPair const parents = m_ParentChooser->chooseParents();
			dad = parents.first;
			mom = parents.second;

			m_OffspringGenerator->generateOffspring(pop, offPop, dad, mom, it, offEnd);
		}
	}
	else
	{
		DBG_DO(DBG_MATING, cerr << "Mating is done in " << numThreads() << " threads" << endl);
		// in this case, openMP must have been supported with numThreads() > 1
#ifdef _OPENMP
		size_t offPopSize = offEnd - offBegin;
		ssize_t nBlocks = numThreads() * 2;
		ssize_t numOffspring = m_OffspringGenerator->numOffspring(pop.gen());
		int except = 0;
		string msg;
#pragma omp parallel for
		for (int i = 0; i < nBlocks; i++)
		{
			try
			{
				RawIndIterator local_it = offBegin + i * (offPopSize / nBlocks / numOffspring) * numOffspring;
				RawIndIterator local_offEnd = i == nBlocks - 1 ? offEnd : local_it + (offPopSize / nBlocks / numOffspring) * numOffspring;

				while (local_it != local_offEnd)
				{
					if (except)
						break;
					Individual *dad = NULL;
					Individual *mom = NULL;
					ParentChooser::IndividualPair const parents = m_ParentChooser->chooseParents();
					dad = parents.first;
					mom = parents.second;
					m_OffspringGenerator->generateOffspring(pop, offPop, dad, mom, local_it, local_offEnd);
				}
			}
			catch (StopEvolution e)
			{
				if (!except)
				{
					except = 1;
					msg = e.message();
				}
			}
			catch (ValueError e)
			{
				if (!except)
				{
					except = 2;
					msg = e.message();
				}
			}
			catch (RuntimeError e)
			{
				if (!except)
				{
					except = 3;
					msg = e.message();
				}
			}
			catch (Exception e)
			{
				if (!except)
				{
					except = 4;
					msg = e.message();
				}
			}
			catch (...)
			{
				if (!except)
					except = -1;
			}
		}

		if (except == 1)
			throw StopEvolution(msg);
		else if (except == 2)
			throw ValueError(msg);
		else if (except == 3)
			throw RuntimeError(msg);
		else if (except == 4)
			throw Exception(msg);
		else if (except == -1)
			throw Exception("Unexpected error from openMP parallel region");
#endif
	}
	m_ParentChooser->finalize();
	m_OffspringGenerator->finalize(pop);
	return true;
}

bool PedigreeMating::mate(Population &pop, Population &scratch)
{
	if (m_gen == -1)
		return false;

	// scrtach will have the right structure.
	if (scratch.genoStruIdx() != pop.genoStruIdx())
		scratch.fitGenoStru(pop.genoStruIdx());
	//
	size_t oldGen = m_ped.curAncestralGen();
	const_cast<Pedigree &>(m_ped).useAncestralGen(m_gen);
	DBG_DO(DBG_MATING, cerr << "Producing offspring generation of size " << m_ped.subPopSizes() << " using generation " << m_gen << " of the pedigree." << endl);
	scratch.fitSubPopStru(m_ped.subPopSizes(), m_ped.subPopNames());
	scratch.setVirtualSplitter(pop.virtualSplitter());
	scratch.clearInfo();

	// build an index for parents
	IdMap idMap;
	size_t idIdx = pop.infoIdx(m_idField);
	RawIndIterator it = pop.rawIndBegin();
	RawIndIterator it_end = pop.rawIndEnd();
	for (; it != it_end; ++it)
		idMap[toID(it->info(idIdx))] = &*it;

	// initialize operator before entering parallel region in order to avoid race condition
	opList::const_iterator iop = m_transmitters.begin();
	opList::const_iterator iopEnd = m_transmitters.end();
	for (; iop != iopEnd; ++iop)
		(*iop)->initializeIfNeeded(*pop.rawIndBegin());

#pragma omp parallel private(it, it_end) if (numThreads() > 1 && parallelizable())
	{
#ifdef _OPENMP
		size_t id = omp_get_thread_num();
		size_t offPopSize = scratch.rawIndEnd() - scratch.rawIndBegin();
		it = scratch.rawIndBegin() + id * (offPopSize / numThreads());
		it_end = id == numThreads() - 1 ? scratch.rawIndEnd() : it + (offPopSize / numThreads());
		size_t i = id * (offPopSize / numThreads());
#else
		it = scratch.rawIndBegin();
		it_end = scratch.rawIndEnd();
		size_t i = 0;
#endif
		for (; it != it_end; ++it, ++i)
		{
			const Individual &pedInd = m_ped.individual(static_cast<double>(i));

			size_t my_id = toID(pedInd.info(m_ped.idIdx()));
			size_t father_id = m_ped.fatherOf(my_id);
			size_t mother_id = m_ped.motherOf(my_id);
			Individual *dad = NULL;
			Individual *mom = NULL;

			if (father_id)
			{
				IdMap::iterator dad_it = idMap.find(father_id);
				DBG_FAILIF(dad_it == idMap.end(), RuntimeError,
						   (boost::format("Could not locate individual with ID %1%") % father_id).str());
				dad = &*(dad_it->second);
			}
			if (mother_id)
			{
				IdMap::iterator mom_it = idMap.find(mother_id);
				DBG_FAILIF(mom_it == idMap.end(), RuntimeError,
						   (boost::format("Could not locate individual with ID %1%") % mother_id).str());
				mom = &*(mom_it->second);
			}
			DBG_DO(DBG_MATING, cerr << "Choosing parents " << father_id << " and "
									<< mother_id << " for offspring " << my_id << endl);

			// copy sex
			it->setSex(pedInd.sex());
			// copy id
			it->setInfo(static_cast<double>(my_id), m_idField);
			//
			opList::const_iterator iop = m_transmitters.begin();
			opList::const_iterator iopEnd = m_transmitters.end();
			for (; iop != iopEnd; ++iop)
			{
				if ((*iop)->isActive(pop.rep(), pop.gen()))
					(*iop)->applyDuringMating(pop, scratch, it, dad, mom);
			}
			// copy individual ID again, just to make sure that even if during mating operators
			// changes ID, pedigree mating could proceed normally.
			it->setInfo(static_cast<double>(my_id), m_idField);
		}
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

bool PedigreeMating::parallelizable() const
{
	opList::const_iterator iop = m_transmitters.begin();
	opList::const_iterator iopEnd = m_transmitters.end();

	for (; iop != iopEnd; ++iop)
	{
		if (!(*iop)->parallelizable())
			return false;
	}
	return true;
}

HeteroMating::HeteroMating(const vectormating &matingSchemes,
						   const uintListFunc &subPopSize,
						   bool shuffleOffspring, SexChoice weightBy)
	: MatingScheme(subPopSize),
	  m_shuffleOffspring(shuffleOffspring),
	  m_weightBy(weightBy)
{
	vectormating::const_iterator it = matingSchemes.begin();
	vectormating::const_iterator it_end = matingSchemes.end();

	for (; it != it_end; ++it)
		m_matingSchemes.push_back(dynamic_cast<HomoMating *>((*it)->clone()));
}

string HeteroMating::describe(bool format) const
{
	string desc = (boost::format("<simuPOP.HeteroMating> a heterogeneous mating scheme with %1% homogeneous mating schemes:\n<ul>\n") % m_matingSchemes.size()).str();
	vectormating::const_iterator it = m_matingSchemes.begin();
	vectormating::const_iterator it_end = m_matingSchemes.end();

	for (; it != it_end; ++it)
	{
		desc += "<li>" + dynamic_cast<HomoMating *>(*it)->describe(false) +
				"<indent>in ";
		subPopList subPops = (*it)->subPops();
		if (subPops.allAvail())
			desc += "all subpopulations.\n";
		else
		{
			desc += "subpopulations ";
			for (size_t i = 0; i < subPops.size(); ++i)
			{
				vspID sp = subPops[i];
				if (i != 0)
					desc += ", ";
				if (sp.isVirtual())
				{
					desc += (boost::format("(%1%, %2%)") % (sp.allAvailSP() ? "ALL_AVAIL" : (boost::format("%1%") % sp.subPop()).str()) %
							 (sp.allAvailVSP() ? "ALL_AVAIL" : (boost::format("%1%") % sp.virtualSubPop()).str()))
								.str();
				}
				else
					desc += (boost::format("%1%") % sp.subPop()).str();
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

HeteroMating::HeteroMating(const HeteroMating &rhs) : MatingScheme(rhs), m_shuffleOffspring(rhs.m_shuffleOffspring), m_weightBy(rhs.m_weightBy)
{
	vectormating::const_iterator it = rhs.m_matingSchemes.begin();
	vectormating::const_iterator it_end = rhs.m_matingSchemes.end();

	for (; it != it_end; ++it)
	{
		m_matingSchemes.push_back(dynamic_cast<HomoMating *>((*it)->clone()));
		DBG_WARNIF(dynamic_cast<HomoMating *>(*it)->subPopSizeSpecified(),
				   "Parameter subPopSize of a HomoMating is ignored when this mating"
				   " scheme is used in a heterogeneous mating scheme.");
	}
}

bool HeteroMating::mate(Population &pop, Population &scratch)
{
	// scrtach will have the right structure.
	if (!prepareScratchPop(pop, scratch))
		return false;

	for (size_t sp = 0; sp < static_cast<size_t>(pop.numSubPop()); ++sp)
	{
		vectormating m;
		vectorf w_pos;  // positive weights
		vectorf w_neg;  // negative weights
		subPopList sps; // each subpopulations
		//
		vectormating::iterator it = m_matingSchemes.begin();
		vectormating::iterator it_end = m_matingSchemes.end();
		for (; it != it_end; ++it)
		{
			subPopList subPops = (*it)->subPops().expandFrom(pop);
			subPopList::const_iterator vsp = subPops.begin();
			subPopList::const_iterator vspEnd = subPops.end();
			for (; vsp != vspEnd; ++vsp)
			{
				if (vsp->subPop() != sp)
					continue;
				// if it is used for this subpop, or all subpopulations
				m.push_back(*it);
				sps.push_back(*vsp);
				double w = (*it)->weight();
				// less than zero...
				if (fcmp_lt(w, 0.))
				{
					w_pos.push_back(0);
					w_neg.push_back(-w);
				}
				else
				{
					w_pos.push_back(w);
					w_neg.push_back(0);
				}
			}
		}
		DBG_FAILIF(m.empty(), ValueError,
				   (boost::format("No mating scheme is available for subpopulation %1%") % sp).str());
		// determine the weight
		if (m.size() == 1)
		{
			w_pos[0] = 1.;
			w_neg[0] = 0.;
		}

		vectoru vspSize(m.size());
		vectoru parentSize(m.size());
		// the default case (all zero)
		bool all_zero = fcmp_eq(std::accumulate(w_pos.begin(), w_pos.end(), 0.), 0.);
		for (size_t i = 0; i < m.size(); ++i)
		{
			parentSize[i] = pop.subPopSize(sps[i], -1, m_weightBy);
			if (all_zero)
			{
				// if there is no negative weight, use population size as weight
				if (w_neg[i] == 0)
					w_pos[i] = static_cast<double>(parentSize[i]);
			}
			else
			{
				if (parentSize[i] == 0)
				{
					DBG_WARNIF(parentSize[i] == 0, "WARNING: One of the parental (virtual) subpopulation is empty and will not "
												   "produce any offspring.");
					w_pos[i] = 0;
					w_neg[i] = 0;
				}
			}
		}
		DBG_DO(DBG_DEVEL, cerr << "Parental Population Size: " << parentSize << "\n"
							   << "Positive mating scheme weights: " << w_pos << '\n'
							   << "Negative mating scheme weights: " << w_neg << endl);

		// weight.
		double overall_pos = std::accumulate(w_pos.begin(), w_pos.end(), 0.);

		// if it happens to be the case that positive weight happens at empty population
		// we need to distribute the weight to other populations.
		if (!all_zero && overall_pos == 0) {
			for (size_t i = 0; i < m.size(); ++i)
			{
				// if there is no negative weight, use population size as weight
				if (w_neg[i] == 0)
					w_pos[i] = static_cast<double>(parentSize[i]);
			}
		}
		// re-calculate overall pos
		overall_pos = std::accumulate(w_pos.begin(), w_pos.end(), 0.);
		double overall_neg = std::accumulate(w_neg.begin(), w_neg.end(), 0.);
		(void)overall_neg; // silent warning about unused variable.
		//
		size_t all = scratch.subPopSize(sp);

		DBG_FAILIF(all > 0 && fcmp_eq(overall_pos, 0.) && fcmp_eq(overall_neg, 0.), ValueError,
				   "No valid parents to produce an non-empty offspring population.");
		// first count negative ones
		for (size_t i = 0; i < m.size(); ++i)
		{
			if (fcmp_gt(w_neg[i], 0.))
			{
				vspSize[i] = static_cast<ULONG>(parentSize[i] * w_neg[i]);
				DBG_ASSERT(all >= vspSize[i], ValueError,
						   (boost::format("Mating scheme with a negative weight of %1% would like to produce %2%"
										  " offspring, but there are only %3% unclaimed offspring left.") %
							w_neg[i] % vspSize[i] % all)
							   .str());
				all -= vspSize[i];
			}
		}
		// then count positive ones
		size_t all_pos = all;
		for (size_t i = 0; i < m.size(); ++i)
		{
			if (all > 0 && fcmp_gt(w_pos[i], 0.))
			{
				vspSize[i] = static_cast<ULONG>(all_pos * w_pos[i] / overall_pos);
				DBG_ASSERT(all >= vspSize[i], ValueError,
						   (boost::format("Mating scheme with a positive weight of %1% would like to produce %2%"
										  " offspring, but there are only %1% unclaimed offspring left.") %
							w_pos[i] % vspSize[i] % all)
							   .str());
				all -= vspSize[i];
			}
		}
		DBG_FAILIF(fcmp_eq(overall_pos, 0) && all > 0, ValueError,
				   "An exact (all negative) weight system is used, but does not fill offspring subpopulation.");

		// individuals left by floating point calculation is added to
		// the last non-zero, positive weight virtual subpopulation.
		if (all > 0)
		{
			for (ssize_t i = m.size() - 1; i >= 0; --i)
				if (vspSize[i] != 0 && w_pos[i] > 0)
				{
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
				   (boost::format("SubPopulation %1% has activated virtual subpopulation.") % sp).str());
		for (UINT idx = 0; idx < m.size(); ++idx, ++itSize)
		{
			DBG_WARNIF(*itSize == 0, "WARNING: One of the mating schemes has zero weight and produces no offspring. "
									 "Because the default weight of a mating scheme is 0, which is handled differently "
									 "when all weights are zero (proportion to sizes of parental subpopulations or virtual "
									 "subpopulations) and when there is a positive weight (weight zero, no offspring). You "
									 "might have forgotten to assign a weight to a mating scheme when you change the weight "
									 "of another mating scheme.");
			if (*itSize == 0)
				continue;
			if (sps[idx].isVirtual())
				pop.activateVirtualSubPop(sps[idx]);
			// if previous mating scheme works on a virtual subpop,
			// and the current one is not. deactivate it.
			else if (pop.hasActivatedVirtualSubPop(sp))
				pop.deactivateVirtualSubPop(sp);

			// real mating
			try
			{
				if (!m[idx]->mateSubPop(pop, scratch, sp, ind, ind + *itSize))
					return false;
			}
			catch (Exception &)
			{
				cerr << "Mating scheme " << idx << " in subpopulation " << sp << " failed to produce " << (*itSize) << " offspring." << endl;
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
		if (m.size() > 1 && m_shuffleOffspring)
		{
			DBG_DO(DBG_MATING, cerr << "Random shuffle individuals in the offspring generation." << endl);
			getRNG().randomShuffle(scratch.rawIndBegin(sp), scratch.rawIndEnd(sp));
			scratch.setIndOrdered(false);
		}
	} // each subpopulation.
	submitScratch(pop, scratch);
	return true;
}

ConditionalMating::ConditionalMating(PyObject *cond, const MatingScheme &ifMatingScheme,
									 const MatingScheme &elseMatingScheme)
#if PY_VERSION_HEX >= 0x03000000
	: m_cond(PyUnicode_Check(cond) ? PyObj_AsString(cond) : string()),
#else
	: m_cond((PyString_Check(cond) || PyUnicode_Check(cond)) ? PyObj_AsString(cond) : string()),
#endif
	  m_func(PyCallable_Check(cond) ? cond : NULL),
	  m_fixedCond(-1), m_ifMS(NULL), m_elseMS(NULL)
{
#if PY_VERSION_HEX >= 0x03000000
	if (!PyUnicode_Check(cond) && !PyCallable_Check(cond))
	{
#else
	if (!PyString_Check(cond) && !PyUnicode_Check(cond) && !PyCallable_Check(cond))
	{
#endif
		bool c;
		PyObj_As_Bool(cond, c);
		m_fixedCond = c ? 1 : 0;
	}

	m_ifMS = ifMatingScheme.clone();
	m_elseMS = elseMatingScheme.clone();
} // namespace simuPOP

string ConditionalMating::describe(bool format) const
{
	string ifDesc = m_ifMS->describe(format);
	string elseDesc = m_elseMS->describe(format);
	string desc = "<simuPOP.ConditionalMating> a conditional mating scheme that ";

	if (m_fixedCond != -1)
		desc += "always applies mating scheme \n" + (m_fixedCond == 1 ? ifDesc : elseDesc);
	else if (m_func.isValid())
		desc += "applies mating scheme \n" + ifDesc + "\n<indent>if function " + m_func.name() + " returns True, and otherwise apply \n" + elseDesc + "\n";
	else
	{
		desc += "applies mating scheme \n" + ifDesc + "\n<indent>if " + m_cond.expr() + " returns True, and otherwise apply\n" + elseDesc + "\n";
	}

	return format ? formatDescription(desc) : desc;
}

ConditionalMating::~ConditionalMating()
{
	delete m_ifMS;
	delete m_elseMS;
}

ConditionalMating::ConditionalMating(const ConditionalMating &rhs) : m_cond(rhs.m_cond), m_func(rhs.m_func), m_fixedCond(rhs.m_fixedCond),
																	 m_ifMS(NULL), m_elseMS(NULL)
{
	m_ifMS = rhs.m_ifMS->clone();
	m_elseMS = rhs.m_elseMS->clone();
}

bool ConditionalMating::mate(Population &pop, Population &scratch)
{
	bool res = true;

	if (m_fixedCond != -1)
		res = m_fixedCond == 1;
	else if (m_func.isValid())
	{
		PyObject *args = PyTuple_New(m_func.numArgs());

		DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

		for (size_t i = 0; i < m_func.numArgs(); ++i)
		{
			const string &arg = m_func.arg(i);
			if (arg == "pop")
				PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(&pop)));
			else
			{
				DBG_FAILIF(true, ValueError, "Only parameter 'pop' is acceptable in "
											 "function" +
												 m_func.name());
			}
		}
		res = m_func(PyObj_As_Bool, args);
		Py_XDECREF(args);
	}
	else
	{
		m_cond.setLocalDict(pop.dict());
		res = m_cond.valueAsBool();
	}

	if (res)
		return m_ifMS->mate(pop, scratch);
	else
		return m_elseMS->mate(pop, scratch);
	return true;
}

} // namespace simuPOP
