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

offspringGenerator::offspringGenerator(const baseOperator & transmitter,
	UINT numParents, double numOffspring, PyObject * numOffspringFunc,
	UINT numOffspringParam, UINT mode, double sexParam, UINT sexMode) :
	m_numOffspring(numOffspring), m_numOffspringFunc(NULL),
	m_numOffspringParam(numOffspringParam), m_mode(mode),
	m_sexParam(sexParam), m_sexMode(sexMode),
	m_numParents(numParents), m_transmitter(NULL),
	m_formOffGenotype(true), m_initialized(false)
{
	DBG_FAILIF(mode == MATE_PyNumOffspring && numOffspringFunc == NULL, ValueError,
		"Please provide a python function when mode is MATE_PyNumOffspring");
	if (numOffspringFunc != NULL && numOffspringFunc != Py_None) {
		DBG_ASSERT(PyCallable_Check(numOffspringFunc), ValueError,
			"Passed variable is not a callable python function.");

		Py_INCREF(numOffspringFunc);
		m_numOffspringFunc = numOffspringFunc;
		m_mode = MATE_PyNumOffspring;
	}
	DBG_FAILIF(mode == MATE_BinomialDistribution && numOffspringParam < 2,
		ValueError, "If mode is MATE_BinomialDistribution, numOffspringParam should be > 1");
	DBG_FAILIF(mode == MATE_UniformDistribution && numOffspringParam < static_cast<UINT>(numOffspring),
		ValueError, "If mode is MATE_UniformDistribution, numOffspringParam should be greater than numOffspring");
	DBG_FAILIF(mode == MATE_GeometricDistribution && (fcmp_lt(m_numOffspring, 0) || fcmp_gt(m_numOffspring, 1.)),
		ValueError, "P for a geometric distribution should be within [0,1], given " + toStr(m_numOffspring));
	DBG_FAILIF(m_mode == MATE_BinomialDistribution && (fcmp_lt(m_numOffspring, 0) || fcmp_gt(m_numOffspring, 1.)),
		ValueError, "P for a Bionomial distribution should be within [0,1], given " + toStr(m_numOffspring));
	DBG_FAILIF(m_mode == MATE_BinomialDistribution && m_numOffspringParam < 1,
		ValueError, "Max number of offspring should be greater than 1. Given "
		+ toStr(m_numOffspringParam));
	DBG_FAILIF(m_sexMode == MATE_ProbOfMale && (fcmp_lt(m_sexParam, 0) || fcmp_gt(m_sexParam, 1)),
		ValueError, "Probability of male has to be between 0 and 1");
	// the genotype transmitter that will be used when no during mating
	// operator is used to transmit genotype from parents to offspring.
	m_transmitter = transmitter.clone();
}


offspringGenerator::offspringGenerator(const offspringGenerator & rhs)
	: m_numOffspring(rhs.m_numOffspring),
	m_numOffspringFunc(rhs.m_numOffspringFunc),
	m_numOffspringParam(rhs.m_numOffspringParam),
	m_mode(rhs.m_mode),
	m_sexParam(rhs.m_sexParam),
	m_sexMode(rhs.m_sexMode),
	m_numParents(rhs.m_numParents),
	m_transmitter(rhs.m_transmitter),
	m_formOffGenotype(rhs.m_formOffGenotype),
	m_initialized(rhs.m_initialized)
{
	if (m_numOffspringFunc != NULL)
		Py_INCREF(m_numOffspringFunc);

	if (rhs.m_transmitter != NULL)
		m_transmitter = rhs.m_transmitter->clone();
}


ULONG offspringGenerator::numOffspring(int gen)
{
	switch (m_mode) {
	case MATE_NumOffspring:
		return static_cast<UINT>(m_numOffspring);
	case MATE_PyNumOffspring: {
		int numOff;
		PyCallFunc(m_numOffspringFunc, "(i)", gen, numOff, PyObj_As_Int);
		DBG_FAILIF(numOff < 1, ValueError, "Need at least one offspring.");
		return numOff;
	}
	case MATE_GeometricDistribution:
		return rng().randGeometric(m_numOffspring);
	case MATE_PoissonDistribution:
		return rng().randPoisson(m_numOffspring) + 1;
	case MATE_BinomialDistribution:
		return rng().randBinomial(m_numOffspringParam - 1, m_numOffspring) + 1;
	case MATE_UniformDistribution:
		// max: 5
		// num: 2
		// randint(4)  ==> 0, 1, 2, 3
		// + 2 ==> 2, 3, 4, 5
		return rng().randInt(m_numOffspringParam - static_cast<unsigned long>(m_numOffspring) + 1)
		       + static_cast<UINT>(m_numOffspring);
	default:
		throw ValueError("Wrong mating numoffspring mode. Should be one of \n"
			             "MATE_NumOffspring, MATE_NumOffspringEachFamily and MATE_GEometricDistribution");
	}
	//
	DBG_ASSERT(false, SystemError, "This line should never be reached.");
	return 0;
}


Sex offspringGenerator::getSex(int count)
{
	if (m_sexMode == MATE_RandomSex)
		return rng().randInt(2) == 0 ? Male : Female;
	else if (m_sexMode == MATE_ProbOfMale)
		return rng().randUniform01() < m_sexParam ? Male : Female;
	else if (m_sexMode == MATE_NumOfMale)
		return count < static_cast<int>(m_sexParam) ? Male : Female;
	else if (m_sexMode == MATE_NumOfFemale)
		return count < static_cast<int>(m_sexParam) ? Female : Male;
	DBG_ASSERT(false, SystemError, "This line should not be reached.");
	return Male;
}


void offspringGenerator::initialize(const population & pop, vector<baseOperator *> const & ops)
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
	while (count < numOff && it != itEnd) {
		// set sex, during mating operator will try to
		// follow the offspring sex (e.g. pass X or Y chromosome)
		it->setSex(getSex(count));
		// use the default offspring genotype transmitter.
		if (m_formOffGenotype)
			m_transmitter->applyDuringMating(pop, it, dad, mom);

		accept = true;
		// apply during mating operators
		vector<baseOperator *>::iterator iop = ops.begin();
		vector<baseOperator *>::iterator iopEnd = ops.end();
		for (; iop != iopEnd;  ++iop) {
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


void sequentialParentChooser::initialize(population & pop, SubPopID sp)
{
	m_begin = pop.indBegin(sp);
	m_end = pop.indEnd(sp);
	m_ind = m_begin;
	m_initialized = true;
}


individual * sequentialParentChooser::chooseParent(RawIndIterator)
{
	if (m_ind == m_end)
		m_ind = m_begin;
	return & * m_ind++;
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


// void pedigreeParentsChooser::initialize(population & pop, SubPopID subPop)
// {
//  m_gen = pop.gen();
//  m_subPop = subPop;
//  m_begin = pop.rawIndBegin();
//  m_index = 0;
//  m_initialized = true;
// }
//
//
// parentChooser::individualPair pedigreeParentsChooser::chooseParents(RawIndIterator)
// {
//  individual * dad = NULL;
//  individual * mom = NULL;
//
//  dad = & * (m_begin + m_pedigree.father(m_gen, m_subPop, m_index));
//  if (m_pedigree.numParents() == 2)
//      mom = & * (m_begin + m_pedigree.mother(m_gen, m_subPop, m_index));
//  DBG_FAILIF(m_index >= m_pedigree.subPopSize(m_gen, m_subPop), IndexError,
//      "Trying to retrieve more indiviudals (index=" +
//      toStr(m_index) + " than what are available from the pedigree ("
//      + toStr(m_pedigree.subPopSize(m_gen, m_subPop)) + ")");
//  m_index++;
//  return std::make_pair(dad, mom);
// }


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
        if (m_size == 0) // if m_index is not used (no VSP)
            m_size = pop.subPopSize(sp);
    }

	if (!m_replacement)
		std::random_shuffle(m_index.begin(), m_index.end());

	DBG_FAILIF(!m_replacement && m_selection, ValueError,
		"Selection is not allowed in random sample without replacement");

	m_initialized = true;
}


individual * randomParentChooser::chooseParent(RawIndIterator basePtr)
{
	DBG_ASSERT(initialized(), SystemError,
		"Please initialize this parent chooser before using it");
	// choose a parent
	if (!m_replacement) {
		DBG_FAILIF(m_index.empty(), ValueError,
			"All parents have been chosen.");
		individual * ind = & * m_index.back();
		m_index.pop_back();
		return ind;
	}
    if (m_index.empty()) {
    	if (m_selection)
	    	return & * (basePtr + m_sampler.get());
    	else
	    	return & * (basePtr + rng().randInt(m_size));
    } else  {
    	if (m_selection)
	    	return & * (m_index[m_sampler.get()]);
    	else
	    	return & * (m_index[rng().randInt(m_size)]);
    }
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
			"All females has been chosen.");
		mom = & * m_femaleIndex.back();
		m_femaleIndex.pop_back();
		
        DBG_FAILIF(m_maleIndex.empty(), ValueError,
		    "All males has been chosen.");
		dad = & * m_maleIndex.back();
		m_maleIndex.pop_back();
		return std::make_pair(dad, mom);
	}

	// using weidhted sampler.
    if (m_selection) {                                    // with selection
        if (m_numMale != 0)
            dad = & * (m_maleIndex[m_malesampler.get()]);
        else
            dad = & * (m_femaleIndex[m_femalesampler.get()]);
    } else {
        // using random sample.
        if (m_numMale != 0)
            dad = & * (m_maleIndex[rng().randInt(m_numMale)]);
        else
            dad = & * (m_femaleIndex[rng().randInt(m_numFemale)]);
    }

    if (m_selection) {                                    // with selection
        if (m_numFemale != 0)
            mom = & * (m_femaleIndex[m_femalesampler.get()]);
        else
            mom = & * (m_maleIndex[m_malesampler.get()]);
    } else {         // no selection
        if (m_numFemale != 0)
            mom = & * (m_femaleIndex[rng().randInt(m_numFemale)]);
        else
            mom = & * (m_maleIndex[rng().randInt(m_numMale)]);
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
		if (m_selection) {                                    // with selection
			if (m_numMale != 0)
				dad = & * (m_maleIndex[m_malesampler.get()]);
			else
				dad = & * (m_femaleIndex[m_femalesampler.get()]);
		} else {
			// using random sample.
			if (m_numMale != 0)
				dad = & * (m_maleIndex[rng().randInt(m_numMale)]);
			else
				dad = & * (m_femaleIndex[rng().randInt(m_numFemale)]);
		}
		if (m_polySex == Male && m_polyNum > 1) {
			m_polyCount = m_polyNum - 1;
			m_lastParent = dad;
		}
	}

	if (mom == NULL) {
		if (m_selection) {                                    // with selection
			if (m_numFemale != 0)
				mom = & * (m_femaleIndex[m_femalesampler.get()]);
			else
				mom = & * (m_maleIndex[m_malesampler.get()]);
		} else {         // no selection
			if (m_numFemale != 0)
				mom = & * (m_femaleIndex[rng().randInt(m_numFemale)]);
			else
				mom = & * (m_maleIndex[rng().randInt(m_numMale)]);
		}
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
		if (m_selection) { // fix me, without replacement!
			// using weighted sampler.
			m_newAlphaIndex.push_back(
				m_alphaSex == Male ? m_maleIndex[m_malesampler.get()]
				: m_femaleIndex[m_femalesampler.get()]);
			m_newAlphaFitness.push_back(m_newAlphaIndex.back()->info(fit_id));
		} else  // fix me, without replacement!
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

	// using weidhted sampler.
    if (m_selection) {                                    // with selection
        if (m_numMale != 0)
            dad = & * (m_maleIndex[m_malesampler.get()]);
        else
            dad = & * (m_femaleIndex[m_femalesampler.get()]);
    } else {
        // using random sample.
        if (m_numMale != 0)
            dad = & * (m_maleIndex[rng().randInt(m_numMale)]);
        else
            dad = & * (m_femaleIndex[rng().randInt(m_numFemale)]);
    }

    if (m_selection) {                                    // with selection
        if (m_numFemale != 0)
            mom = & * (m_femaleIndex[m_femalesampler.get()]);
        else
            mom = & * (m_maleIndex[m_malesampler.get()]);
    } else {         // no selection
        if (m_numFemale != 0)
            mom = & * (m_femaleIndex[rng().randInt(m_numFemale)]);
        else
            mom = & * (m_maleIndex[rng().randInt(m_numMale)]);
    }
	return std::make_pair(dad, mom);
}


void infoParentsChooser::initialize(population & pop, SubPopID sp)
{
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
	individual * par1 = chooseParent(basePtr);
	Sex sex1 = par1->sex();
	// there is no valid information field value
	if (m_degenerate) {
		int attempt = 0;
		while (++attempt < 1000) {
			individual * par2 = chooseParent(basePtr);
			if (par2->sex() != sex1)
				return sex1 == Male ? std::make_pair(par1, par2) : std::make_pair(par2, par1);
		}
		throw ValueError("Can not locate any individual of opposite sex");
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
	DBG_FAILIF(validInds.empty(), SystemError,
		"No valid relative is found");
	individual * par2 = validInds[rng().randInt(validInds.size())];
	DBG_DO(DBG_DEVEL, cout << "infoParentsChooser: par1: " << par1 - & * basePtr
		                   << " par2: " << par2 - & * basePtr << endl);
	return sex1 == Male ? std::make_pair(par1, par2) : std::make_pair(par2, par1);
}


pyParentsChooser::pyParentsChooser(PyObject * pc)
	: parentChooser(0), m_func(pc), m_generator(NULL), m_parIterator(NULL)
{
	Py_INCREF(m_func);
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
	PyObject * arglist = Py_BuildValue("(Oi)", popObj, sp);
	m_generator = PyEval_CallObject(m_func, arglist);
	Py_XDECREF(arglist);

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
	DBG_FAILIF(item == NULL, ValueError,
		"User-defined function yield invalid value. This may happen \n"
		"if you function does not provide enough parents for the mating \n"
		"scheme (a 'while True' statement is recommended).");

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


mating::mating(vectorlu newSubPopSize, string newSubPopSizeExpr, PyObject * newSubPopSizeFunc,
	vspID subPop, double weight)
	: m_subPopSize(newSubPopSize),
	m_subPopSizeExpr(newSubPopSizeExpr, ""), m_subPopSizeFunc(NULL),
	m_subPop(subPop), m_weight(weight)
{
	DBG_FAILIF(!m_subPopSizeExpr.empty() && newSubPopSizeFunc != NULL,
		ValueError, "Please only specify one of newSubPopSizeExpr and newSubPopSizeFunc.");


	if (newSubPopSizeFunc != NULL && newSubPopSizeFunc != Py_None) {
		if (!PyCallable_Check(newSubPopSizeFunc))
			throw ValueError("Passed variable is not a callable python function.");

		Py_INCREF(newSubPopSizeFunc);
		m_subPopSizeFunc = newSubPopSizeFunc;
	}
}


mating::mating(const mating & rhs)
	: m_subPopSize(rhs.m_subPopSize),
	m_subPopSizeExpr(rhs.m_subPopSizeExpr),
	m_subPopSizeFunc(rhs.m_subPopSizeFunc),
	m_subPop(rhs.m_subPop),
	m_weight(rhs.m_weight)
{
	if (m_subPopSizeFunc != NULL)
		Py_INCREF(m_subPopSizeFunc);
#ifndef OPTIMIZED
	m_famSize.clear();
#endif
}


void mating::prepareScratchPop(population & pop, population & scratch)
{
	// use population structure of pop
	if (m_subPopSize.empty() && m_subPopSizeExpr.empty() && m_subPopSizeFunc == NULL)
		scratch.fitSubPopStru(pop.subPopSizes(), pop.subPopNames());
	else if (!m_subPopSize.empty())  // set subPoplation size
		scratch.fitSubPopStru(m_subPopSize, pop.subPopNames());
	// evaluate from an expression
	else if (!m_subPopSizeExpr.empty()) {
		m_subPopSizeExpr.setLocalDict(pop.dict());
		vectorf sizef = m_subPopSizeExpr.valueAsArray();
		vectorlu sz(sizef.size());

		for (size_t i = 0, iEnd = sizef.size(); i < iEnd; i++)
			sz[i] = static_cast<ULONG>(sizef[i]);

		scratch.fitSubPopStru(sz, pop.subPopNames());
	} else {                                                                            // use m_subPopSizeFunc
		// get generation number
		int gen = pop.gen();
		// convert current pop size to a tuple
		PyObject * curSize = PyTuple_New(pop.numSubPop());

		DBG_ASSERT(curSize != NULL, SystemError, "Can not convert current pop size to a list");

		for (size_t i = 0; i < pop.numSubPop(); ++i)
			PyTuple_SetItem(curSize, i, PyInt_FromLong(pop.subPopSize(i)));

		vectorf res;
		PyCallFunc2(m_subPopSizeFunc, "(iO)", gen, curSize, res, PyObj_As_Array);
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
}


bool mating::mate(population & pop, population & scratch,
                  vector<baseOperator * > & ops, bool submit)
{
	// scrtach will have the right structure.
	prepareScratchPop(pop, scratch);
	preparePopulation(pop);

	DBG_DO(DBG_MATING, m_famSize.clear());

	for (SubPopID sp = 0; sp < static_cast<SubPopID>(pop.numSubPop()); ++sp)
		if (!mateSubPop(pop, sp, scratch.rawIndBegin(sp),
				scratch.rawIndEnd(sp), ops))
			return false;
	if (submit)
		submitScratch(pop, scratch);
	return true;
}


void mating::submitScratch(population & pop, population & scratch)
{
	pop.turnOffSelection();
	// use scratch population,
	pop.push(scratch);
	scratch.validate("after push and discard");
	DBG_DO(DBG_MATING, pop.setIntVectorVar("famSizes", m_famSize));
}


void countAlleles(population & pop, int subpop, const vectori & loci, const vectori & alleles,
                  vectorlu & alleleNum)
{
	alleleNum = vectorlu(loci.size(), 0L);
	for (size_t l = 0; l < loci.size(); ++l) {
		int loc = loci[l];
		Allele ale = alleles[l];
		// go through all alleles
		if (subpop == -1) {
			IndAlleleIterator a = pop.alleleBegin(loc);
			IndAlleleIterator aEnd = pop.alleleEnd(loc);
			for (; a != aEnd; ++a)
				if (AlleleUnsigned(*a) == ale)
					alleleNum[l]++;
		} else {
			IndAlleleIterator a = pop.alleleBegin(loc, subpop);
			IndAlleleIterator aEnd = pop.alleleEnd(loc, subpop);
			for (; a != aEnd; ++a)
				if (AlleleUnsigned(*a) == ale)
					alleleNum[l]++;
		}
	}
}


// give expected frequency for the whole population, or all subpopulations
// return expected number of alleles at each subpopulations.
void getExpectedAlleles(population & pop, vectorf & expFreq, const vectori & loci,
                        const vectori & alleles, vectoru & expAlleles)
{
	// determine expected number of alleles of each allele
	// at each subpopulation.
	UINT nLoci = loci.size();
	UINT numSP = pop.numSubPop();
	size_t pldy = pop.ploidy();
	size_t i, p;

	for (i = 0; i < nLoci; ++i) {
		if (expFreq[i] < 0.)
			expFreq[i] = 0.;
		if (expFreq[i] > 1.)
			expFreq[i] = 1.;
	}
	//
	// in the order of Loc1: sp1, sp2, sp3, sloc2: p1, sp2, sp3
	expAlleles.resize(nLoci * numSP);
	// if exp frequencies in subpopulation is not specified.
	//  I need to find the current allele frequencies
	// and use them as proportions for the next generation.
	if (numSP > 1 && expFreq.size() == nLoci) {
		vectorlu totalAlleles(nLoci);
		for (i = 0; i < nLoci; ++i) {
			int locus = loci[i];
			Allele allele = static_cast<Allele>(alleles[i]);

			totalAlleles[i] = static_cast<ULONG>(expFreq[i] * pop.popSize() * pop.ploidy());
			// make sure one allele (our seed :-) exists
			if (expFreq[i] > 0. && totalAlleles[i] == 0)
				totalAlleles[i] = 1;

			// determine the number alleles at each subpopulation.
			vectorf curFreq(numSP);
			ULONG numOfAlleles = 0;
			for (size_t sp = 0; sp < numSP; ++sp) {
				ULONG n = 0;
				if (!pop.indOrdered()) {
					for (IndIterator it = pop.indBegin(sp); it.valid(); ++it)
						for (p = 0; p < pldy; ++p)
							if (it->allele(locus, p) == allele)
								n++;
				} else {
					IndAlleleIterator a = pop.alleleBegin(locus, sp);
					IndAlleleIterator aEnd = pop.alleleEnd(locus, sp);
					for (; a != aEnd; ++a)
						if (AlleleUnsigned(*a) == allele)
							n++;
				}
				numOfAlleles += n;
				curFreq[sp] = double(n) / (pop.subPopSize(sp) * pldy);
			}

			DBG_DO(DBG_MATING, cout << "Cur freq " << curFreq << endl);

			// if there is no alleles
			if (numOfAlleles == 0 && expFreq[i] > 0.)
				throw ValueError("No disease allele exists, but exp allele frequency is greater than 0.\n"
					             " Generation " + toStr(pop.gen()) );

			// calculate exp number of affected offspring in the next generation.
			//
			// step 1: totalsize*expFreq is the total number of disease alleles
			// step 2: assign these alleles to each subpopulation according to a multi-nomial
			// distribution with p_i beging allele frequency at each subpopulation.
			// assign these numbers to each subpopulation
			rng().randMultinomial(static_cast<unsigned int>(pop.popSize() * expFreq[i] * pldy),
				curFreq, expAlleles.begin() + numSP * i);
		}
	}
	// simpler case, one subpopulation, or with gieven allele frequency
	else if (expFreq.size() == numSP * nLoci) {
		for (i = 0; i < nLoci; ++i) {
			for (size_t sp = 0; sp < numSP; ++sp) {
#ifndef OPTIMIZED
				int locus = loci[i];
				Allele allele = alleles[i];
				ULONG n = 0;
				// go through all alleles
				if (!pop.indOrdered()) {
					for (IndIterator it = pop.indBegin(sp); it.valid(); ++it)
						for (p = 0; p < pldy; ++p)
							if (it->allele(locus, p) == allele)
								n++;
				} else {
					IndAlleleIterator a = pop.alleleBegin(locus, sp);
					IndAlleleIterator aEnd = pop.alleleEnd(locus, sp);
					for (; a != aEnd; ++a) {
						if (AlleleUnsigned(*a) == allele)
							n++;
					}
				}

				// if there is no alleles
				if (n == 0 && expFreq[numSP * i + sp] > 0.)
					throw ValueError("No disease allele exists, but exp allele frequency is greater than 0.\n"
						             " Generation " + toStr(pop.gen()) );
#endif
				expAlleles[numSP * i + sp] = static_cast<UINT>(pop.subPopSize(sp) * pldy * expFreq[numSP * i + sp]);
				if (expFreq[numSP * i + sp] > 0. && expAlleles[numSP * i + sp] == 0)
					expAlleles[numSP * i + sp] = 1;
			}
		}                                                                                 // each locus
	} else {
		throw ValueError("Returned expected frequency has wrong length");
	}
}


// //
// // This mating scheme is very complicated, it is similar to randomMating, but it tries to control
// // the mating process so that the total disease allele frequency controlled to pre-specified values.
// //
// bool controlledRandomMating::mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit)
// {
//  // scrtach will have the right structure.
//  prepareScratchPop(pop, scratch);
//
//  size_t pldy = pop.ploidy(), nLoci = m_loci.size();
//  size_t i;
//
//  // expected frequency at each locus
//  vectorf expFreq;
//  PyCallFunc(m_freqFunc, "(i)", pop.gen(), expFreq, PyObj_As_Array);
//  DBG_DO(DBG_MATING, cout << "expected freq " << expFreq << endl);
//
//  // determine expected number of alleles of each allele
//  // at each subpopulation.
//  UINT numSP = pop.numSubPop();
//  vectoru expAlleles(nLoci * numSP);
//  getExpectedAlleles(scratch, expFreq, m_loci, m_alleles, expAlleles);
//  DBG_DO(DBG_MATING, cout << "expected alleles " << expAlleles << endl);
//
//  // whether or not use stack.
//  if (!m_offspringGenerator.initialized())
//      m_offspringGenerator.initialize(pop, ops);
//  bool useStack = m_offspringGenerator.mode() == MATE_NumOffspring;
//  // use to go through offspring generation to count alleles
//  UINT totNumLoci = pop.totNumLoci();
//
//  // controlled random mating happens within each subpopulation
//  for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
//      ULONG spSize = scratch.subPopSize(sp);
//      if (spSize == 0)
//          continue;
//
//      // reset stack each time
//      if (useStack)
//          m_stack = stack<RawIndIterator>();
//
//      // total allowed disease alleles.
//      vectoru totAllele(nLoci);
//      // currently available disease allele. (in the offspring generation)
//      vectoru curAllele(nLoci, 0);
//      // We control allele 1 if expected allele frequency is less than 0.5
//      vector<bool> flip(nLoci, false);
//      for (i = 0; i < nLoci; ++i) {
//          totAllele[i] = expAlleles[sp + numSP * i];
//          if (totAllele[i] > spSize * pldy) {
//              cout << "Warning: number of planned affected alleles exceed population size.";
//              totAllele[i] = spSize * pldy;
//          }
//          if (2 * totAllele[i] > spSize * pldy) {
//              flip[i] = true;
//              totAllele[i] = spSize * pldy - totAllele[i];
//          }
//      }
//
//      randomParentsChooser pc;
//      pc.initialize(pop, sp);
//      if (pc.numMale() == 0 || pc.numFemale() == 0) {
//          if (m_contWhenUniSex)
//              cout << "Warning: the subpopulation is uni-sex. Mating will continue with same-sex mate" << endl;
//          else
//              throw ValueError("Subpopulation becomes uni-sex. Can not continue. \n"
//                               "You can use ignoreParentsSex (do not check parents' sex) or \ncontWhenUnixSex "
//                               "(same sex mating if have to) options to get around this problem.");
//      }
//
//      // it is possible that no disease allele is required
//      // so everyone is accepted
//      bool freqRequMet = true;
//      for (i = 0; i < nLoci; ++i) {
//          if (totAllele[i] > 0) {
//              freqRequMet = false;
//              break;
//          }
//      }
//      // if a family has disease allele.
//      bool stackStage = false;
//      // start!
//      RawIndIterator it = scratch.rawIndBegin(sp);
//      RawIndIterator it_end = scratch.rawIndEnd(sp);
//      RawIndIterator itBegin = it;
//      UINT numOff = 0;
//      // if mor ethan noAAattempt times, no pure homo found,
//      // accept non-homo cases.
//      int AAattempt = 200;
//      while (true) {
//          // the logic is complicated here and I will try to be explicit
//          // if already in the stackStage, it is taken from the stack.
//          if (useStack && stackStage) {
//              if (m_stack.empty()) {
//                  DBG_ASSERT(!freqRequMet, SystemError, "Empty stack should only happen when freq requirement is not met.");
//                  cout << "Warning: frequency requirement is not met, for subpop " << sp << " at generation "
//                       << pop.gen() << ".\nThis is usually caused by multiple high disease allele frequency." << endl;
//                  break;
//              }
//              it = m_stack.top();
//          }
//
//          // randomly choose parents
//          parentChooser::individualPair parents = pc.chooseParents(pop.rawIndBegin());
//
//          // generate m_numOffspring offspring per mating
//          // record family size (this may be wrong for the last family)
//
//          itBegin = it;
//          // generate numOffspring offspring per mating
//          // it moves forward
//          numOff = m_offspringGenerator.generateOffspring(pop, parents.first, parents.second, it, it_end, ops);
//
//          // count alleles in this family
//          // count number of alleles in the family.
//          vectori na(nLoci, 0);
//          bool hasAff = false;
//          // we know that scratch population has ordered linear genotype
//          for (i = 0; i < nLoci; ++i) {
//              GenoIterator ptr = itBegin->genoBegin() + m_loci[i];
//              for (size_t j = 0; j < numOff * pldy; ++j, ptr += totNumLoci) {
//                  if (flip[i] ? (*ptr != static_cast<Allele>(m_alleles[i]))
//                      : (*ptr == static_cast<Allele>(m_alleles[i]))) {
//                      na[i]++;
//                      hasAff = true;
//                  }
//              }
//          }
//
//          if (useStack) {
//              // now check if this family is usable.
//              bool accept = false;
//              // all disease alleles have been satidfied.
//              // only accept unaffected families
//              // otherwise, accept any family that can help.
//              if (freqRequMet) {
//                  // and not it == itEnd, that is to say, we only accept unaffected
//                  if (!hasAff) {
//                      // has AA, so we have less reason to compromise
//                      // but it might still be very difficult to find one,
//                      // so we still allow accepting of "wrong" individuals.
//                      AAattempt = 10000;
//                      accept = true;
//                  }
//                  // tried 100 times, no AA is found.
//                  else if (AAattempt == 0) {
//                      AAattempt = 100;
//                      accept = true;
//                  }
//                  AAattempt--;
//              } else {
//                  // we accept affected helpful ones, and purely unsffected ones,
//                  if (hasAff) {
//                      // has the right kind of mutant?
//                      for (i = 0; i < nLoci; ++i) {
//                          // accept the whole family, if we need this allele
//                          if (curAllele[i] < totAllele[i] && na[i] > 0) {
//                              accept = true;
//                              break;
//                          }
//                      }
//                  } else {
//                      // in the stack stage, we only accept affected.
//                      if (!stackStage)
//                          accept = true;
//                  }
//              }
//
//              // reject this family
//              if (!accept) {
//                  // it relocate to its begin point
//                  // DBG_DO(DBG_MATING, cout << "Reject " << na << endl);
//                  it = itBegin;
//                  continue;
//              }
//              DBG_DO(DBG_DEVEL, cout << "Accept " << na << " CUR " << curAllele << " TOT  " << totAllele << endl);
//
//              // accept this family
//              for (i = 0; i < nLoci; ++i)
//                  curAllele[i] += na[i];
//              // accpet this family, see if all done.
//              if (!freqRequMet) {
//                  freqRequMet = true;
//                  for (i = 0; i < nLoci; ++i) {
//                      if (curAllele[i] < totAllele[i])
//                          freqRequMet = false;
//                  }
//              }
//              if (!stackStage && !freqRequMet && !hasAff) {
//                  // this family is in stack, might be
//                  m_stack.push(itBegin);
//                  DBG_DO(DBG_DEVEL, cout << "Push in stack " << m_stack.size() << endl);
//              }
//              // accepted,
//              if (stackStage) {
//                  m_stack.pop();
//                  DBG_DO(DBG_DEVEL, cout << "Pop index " << m_stack.size() << endl);
//              }
//              // see if break
//              if (it == it_end) {
//                  stackStage = true;
//                  DBG_DO(DBG_MATING, cout << "Stack stage " << m_stack.size() << endl);
//              }
//              if (freqRequMet && stackStage) {
// #ifndef OPTIMIZED
//                  if (debug(DBG_MATING)) {
//                      cout << "Finish generating offspring. subpopulation size: " << spSize << endl;
//                      cout << "Expected:";
//                      for (size_t ii = 0; ii < nLoci; ++ii)
//                          cout << " " << expAlleles[sp + numSP * ii] << " ("
//                               << (expAlleles[sp + numSP * ii] * 1.0 / (spSize * pldy)) << ")";
//                      cout << endl << "Simulated:";
//                      for (size_t ii = 0; ii < nLoci; ++ii)
//                          cout << " " << (flip[ii] ? spSize * pldy - curAllele[ii] : curAllele[ii]) << " ("
//                               << (flip[ii] ? 1. - curAllele[ii] / (1.0 * spSize * pldy) : curAllele[ii] / (1.0 * spSize * pldy)) << ")";
//                      cout << endl;
//                  }
// #endif
//                  break;
//              }
//          } else {        // if do not use stack.
//              // now check if this family is usable.
//              bool accept = false;
//              // all disease alleles have been satidfied.
//              // only accept unaffected families
//              // otherwise, accept any family that can help.
//              if (freqRequMet) {
//                  if (!hasAff) {
//                      // has AA, so no need to compromise
//                      AAattempt = 10000;
//                      accept = true;
//                  }
//                  // tried 100 times, no AA is found.
//                  else if (AAattempt == 0) {
//                      AAattempt = 100;
//                      accept = true;
//                  }
//                  AAattempt--;
//              } else {                                                       // do not use stack
//                  for (i = 0; i < nLoci; ++i) {
//                      // accept the whole family, if we need this allele
//                      if (curAllele[i] < totAllele[i] && na[i] > 0) {
//                          accept = true;
//                          break;
//                      }
//                  }
//              }
//
//              // reject this family
//              if (!accept) {
//                  // it relocate to its begin point
//                  // DBG_DO(DBG_MATING, cout << "Reject " << na << endl);
//                  it = itBegin;
//                  continue;
//              }
//              DBG_DO(DBG_MATING, if (it - scratch.rawIndBegin(sp) < 100) cout << "Accept " << na << endl;);
//
//              // accpet this family, see if all done.
//              for (i = 0; i < nLoci; ++i)
//                  curAllele[i] += na[i];
//              if (!freqRequMet) {
//                  freqRequMet = true;
//                  for (i = 0; i < nLoci; ++i) {
//                      if (curAllele[i] < totAllele[i])
//                          freqRequMet = false;
//                  }
//              }
//              // see if break
//              if (it == it_end) {
//                  if (!freqRequMet)
//                      cout << "Warning: frequency requirement is not met, for subpop " << sp << " at generation " << pop.gen() << endl;
// #ifndef OPTIMIZED
//                  if (!freqRequMet || debug(DBG_MATING)) {
//                      cout << "Subpopulation size " << spSize << "\nExpected:";
//                      for (size_t ii = 0; ii < nLoci; ++ii)
//                          cout << " " << expAlleles[sp + numSP * ii];
//                      cout << endl << "Simulated:";
//                      for (size_t ii = 0; ii < nLoci; ++ii)
//                          cout << " " << (flip[ii] ? spSize * pldy - curAllele[ii] : curAllele[ii]);
//                  }
// #endif
//                  break;
//              }
//          }
//      }                                                                                           // nostack scheme
//  }                                                                                               // each subPop
//
//  if (submit)
//      submitScratch(pop, scratch);
//  return true;
// }


pyMating::pyMating(parentChooser & chooser,
	offspringGenerator & generator,
	vectorlu newSubPopSize,
	string newSubPopSizeExpr,
	PyObject * newSubPopSizeFunc,
	vspID subPop,
	double weight)
	: mating(newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc, subPop, weight)
{
	m_parentChooser = chooser.clone();
	m_offspringGenerator = generator.clone();
	DBG_FAILIF(m_parentChooser->numParents() != 0 && m_offspringGenerator->numParents() != 0 
		&& m_parentChooser->numParents() != m_offspringGenerator->numParents(),
		ValueError, "Imcompatible parent chooser and offspring generator");
}


bool pyMating::mateSubPop(population & pop, SubPopID subPop,
                          RawIndIterator offBegin, RawIndIterator offEnd,
                          vector<baseOperator * > & ops)
{
	// nothing to do.
	if (offBegin == offEnd)
		return true;

	if (!m_parentChooser->initialized())
		m_parentChooser->initialize(pop, subPop);

	if (!m_offspringGenerator->initialized())
		m_offspringGenerator->initialize(pop, ops);

	// generate scratch.subPopSize(sp) individuals.
	RawIndIterator it = offBegin;
	while (it != offEnd) {
		individual * dad = NULL;
		individual * mom = NULL;
		if (m_parentChooser->numParents() == 1)
			dad = m_parentChooser->chooseParent(pop.rawIndBegin());
		// 0 or 2, that is to say, dad or mom can be NULL
		else {
			parentChooser::individualPair const parents = m_parentChooser->chooseParents(pop.rawIndBegin());
			dad = parents.first;
			mom = parents.second;
		}

		//
		UINT numOff = m_offspringGenerator->generateOffspring(pop, dad, mom, it, offEnd, ops);
		(void)numOff;             // silent warning about unused variable.
		// record family size (this may be wrong for the last family)
		DBG_DO(DBG_MATING, m_famSize.push_back(numOff));
	}
	m_parentChooser->finalize(pop, subPop);
	m_offspringGenerator->finalize(pop);
	return true;
}


heteroMating::heteroMating(const vectormating & matingSchemes,
	vectorlu newSubPopSize,
	string newSubPopSizeExpr,
	PyObject * newSubPopSizeFunc,
	bool shuffleOffspring,
	vspID subPop,
	double weight)
	: mating(newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc, subPop, weight),
	m_shuffleOffspring(shuffleOffspring)
{
	DBG_WARNING(subPop.subPop() != InvalidSubPopID || subPop.virtualSubPop() != InvalidSubPopID,
		"Parameter subPop or virtualSubPop is specified, but is ignored.");

	vectormating::const_iterator it = matingSchemes.begin();
	vectormating::const_iterator it_end = matingSchemes.end();

	for (; it != it_end; ++it)
		m_matingSchemes.push_back((*it)->clone());
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
		m_matingSchemes.push_back((*it)->clone());
}


bool heteroMating::mate(population & pop, population & scratch,
                        vector<baseOperator * > & ops, bool submit)
{
	// scrtach will have the right structure.
	prepareScratchPop(pop, scratch);
	vectormating::const_iterator it = m_matingSchemes.begin();
	vectormating::const_iterator it_end = m_matingSchemes.end();
	for (; it != it_end; ++it)
		(*it)->preparePopulation(pop);

	DBG_DO(DBG_MATING, m_famSize.clear());

	for (SubPopID sp = 0; sp < static_cast<SubPopID>(pop.numSubPop()); ++sp) {
		vectormating m;
		vectorf w_pos;                          // positive weights
		vectorf w_neg;                          // negative weights
		vectormating::iterator it = m_matingSchemes.begin();
		vectormating::iterator it_end = m_matingSchemes.end();
		for (; it != it_end; ++it) {
			// if it is used for this subpop,
			// for use for all subPops ...
			if ((*it)->subPop() == sp ||
			    (*it)->subPop() == InvalidSubPopID) {
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
		(void)overall_neg;             // silent warning about unused variable.
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

			DBG_DO(DBG_MATING, (*it)->m_famSize.clear(););
			// real mating
			try {
				if (!(*it)->mateSubPop(pop, sp, ind, ind + *itSize, ops))
					return false;
			} catch (...) {
				cout << "Mating in subpopulation " + toStr(sp) + " failed" << endl;
				throw;
			}
			DBG_DO(DBG_MATING,
				m_famSize.insert(m_famSize.end(), (*it)->m_famSize.begin(),
					(*it)->m_famSize.end()););
			ind += *itSize;
		}
		DBG_ASSERT(ind == scratch.rawIndEnd(sp), SystemError,
			"Mating scheme somehow does not fill the whole offspring population.");
		pop.deactivateVirtualSubPop(sp);
		// if more than two mating schemes working on the same subpopulation,
		// it is better to shuffle offspring afterwards,
		if (m.size() > 1 && m_shuffleOffspring) {
			std::random_shuffle(pop.rawIndBegin(sp), pop.rawIndEnd(sp));
			pop.setIndOrdered(false);
		}
	}                        // each subpopulation.
	if (submit)
		submitScratch(pop, scratch);
	return true;
}


}
