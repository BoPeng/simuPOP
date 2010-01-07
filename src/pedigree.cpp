/**
 *  $File: Pedigree.cpp $
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

#include "pedigree.h"
// #include <fstream>
// using std::ifstream;
// using std::ofstream;
//
// #include <sstream>
// using std::istringstream;
// using std::ws;
//
// #include <numeric>
// using std::accumulate;

namespace simuPOP {

Pedigree::Pedigree(const Population & pop, const uintList & loci,
	const stringList & infoFields, const uintList & ancGens, const string & idField,
	const string & fatherField, const string & motherField)
	: m_idField(idField), m_fatherField(fatherField), m_motherField(motherField),
	m_idIdx(-1), m_fatherIdx(-1), m_motherIdx(-1)
{
	vectorstr extractFields = infoFields.elems();

	if (!m_idField.empty() && find(extractFields.begin(), extractFields.end(), idField) == extractFields.end())
		extractFields.push_back(idField);
	if (!m_fatherField.empty() && find(extractFields.begin(), extractFields.end(), fatherField) == extractFields.end())
		extractFields.push_back(fatherField);
	if (!m_motherField.empty() && find(extractFields.begin(), extractFields.end(), motherField) == extractFields.end())
		extractFields.push_back(motherField);
	//
	Population & ped = pop.extract(loci, extractFields, subPopList(), ancGens);

	swap(ped);

	DBG_FAILIF(m_idField.empty() || !pop.hasInfoField(m_idField), ValueError,
		"A valid ID information field is needed to create a pedigree object.");
	DBG_FAILIF(!m_fatherField.empty() && !pop.hasInfoField(m_fatherField), ValueError,
		"Invalid father information field " + m_fatherField);
	DBG_FAILIF(!m_motherField.empty() && !pop.hasInfoField(m_motherField), ValueError,
		"Invalid mother information field " + m_motherField);

	m_idIdx = static_cast<int>(infoIdx(idField));
	m_fatherIdx = m_fatherField.empty() ? -1 : static_cast<int>(infoIdx(fatherField));
	m_motherIdx = m_motherField.empty() ? -1 : static_cast<int>(infoIdx(motherField));

	// build an ID map
	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		for (IndIterator it = indIterator(); it.valid(); ++it) {
			ULONG id = static_cast<ULONG>(it->info(m_idIdx) + 0.5);
			DBG_WARNING(m_idMap.find(id) != m_idMap.end() && *m_idMap[id] != *it,
				"Different individuals share the same ID " + toStr(id) +
				"so only the latest Individual will be used. If this is an "
				"age-structured population, you may want to remove parental generations.");
			m_idMap[id] = &*it;
		}
	}
}


Pedigree::Pedigree(const Pedigree & rhs) :
	Population(rhs),
	m_idField(rhs.m_idField), m_fatherField(rhs.m_fatherField), m_motherField(rhs.m_motherField),
	m_idIdx(rhs.m_idIdx), m_fatherIdx(rhs.m_fatherIdx), m_motherIdx(rhs.m_motherIdx)
{
}


Pedigree * Pedigree::clone() const
{
	return new Pedigree(*this);
}


UINT Pedigree::numParents()
{
	return static_cast<UINT>(m_fatherIdx != -1) + static_cast<UINT>(m_motherIdx != -1);
}


Individual & Pedigree::indByID(double fid)
{
	// essentially m_idMap(static_cast<ULONG>(fid))
	//
	ULONG id = static_cast<ULONG>(fid + 0.5);

	DBG_FAILIF(fabs(fid - id) > 1e-8, ValueError,
		"individual ID has to be integer (or a double round to full iteger).");

	std::map<ULONG, Individual *>::iterator it = m_idMap.find(id);
	// if still cannot be found, raise an IndexError.
	if (it == m_idMap.end())
		throw IndexError("No individual with ID " + toStr(id) + " could be found.");

	return *it->second;
}


bool Pedigree::acceptableSex(Sex mySex, Sex relSex, SexChoice choice)
{
	return choice == ANY_SEX ||
	       (choice == MALE_ONLY && relSex == MALE) ||
	       (choice == FEMALE_ONLY && relSex == FEMALE) ||
	       (choice == SAME_SEX && relSex == mySex) ||
	       (choice == OPPOSITE_SEX && relSex != mySex);
}


bool Pedigree::acceptableAffectionStatus(bool affected, AffectionStatus choice)
{
	return choice == ANY_AFFECTION_STATUS ||
	       (choice == AFFECTED && affected) ||
	       (choice == UNAFFECTED && !affected);
}


void Pedigree::locateRelatives(RelativeType relType, const vectorstr & resultFields,
                               SexChoice sexChoice, AffectionStatus affectionChoice, int ancGen)
{
	DBG_ASSERT(sexChoice == ANY_SEX || sexChoice == MALE_ONLY || sexChoice == FEMALE_ONLY
		|| sexChoice == SAME_SEX || sexChoice == OPPOSITE_SEX, ValueError,
		"Relative sex can only be one of ANY_SEX, MALE_ONLY, FEMALE_ONLY, SAME_SEX or OPPOSITE_SEX.");

	DBG_ASSERT(affectionChoice == AFFECTED || affectionChoice == UNAFFECTED
		|| affectionChoice == ANY_AFFECTION_STATUS, ValueError,
		"Relative affection status can only be one of AFFECTED, UNAFFECTED and ANY_AFFECTION_STATUS.");

	UINT oldGen = curAncestralGen();
	switch (relType) {
	case SPOUSE:
		locateSpouse(sexChoice, affectionChoice, resultFields, ancGen, false);
		break;
	case OUTBRED_SPOUSE:
		locateSpouse(sexChoice, affectionChoice, resultFields, ancGen, true);
		break;
	case OFFSPRING:
		locateOffspring(sexChoice, affectionChoice, resultFields, ancGen);
		break;
	case SIBLING:
		locateSibling(sexChoice, affectionChoice, resultFields, ancGen);
		break;
	case FULLSIBLING:
		locateFullSibling(sexChoice, affectionChoice, resultFields, ancGen);
		break;
	case COMMON_OFFSPRING:
		locateCommonOffspring(sexChoice, affectionChoice, resultFields, ancGen);
		break;
	default:
		throw ValueError("Unrecognized relative type");
	}
	useAncestralGen(oldGen);
}


void Pedigree::locateSpouse(SexChoice sexChoice, AffectionStatus affectionChoice, const vectorstr & resultFields, int ancGen, bool excludeOutbred)
{
	DBG_ASSERT(numParents() == 2, ValueError,
		"This relative only exists when there are two parents for each indidivual");

	DBG_ASSERT(resultFields.size() >= 1, ValueError,
		"Please provide at least one information field to store SPOUSE Individuals");

	DBG_FAILIF(sexChoice == SAME_SEX, ValueError, "Can not locate spouses with the same sex");

	// if we are using individual ID, the algorithm is slower, but it can handle
	// overlapping generations.
	//
	// What we essentially do is getting all the couples and assign them...
	UINT maxSpouse = resultFields.size();

	vectori spouseIdx(maxSpouse);

	// clear all fields
	UINT topGen = ancGen == -1 ? ancestralGens() : std::min(ancestralGens(), static_cast<UINT>(ancGen));
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (size_t i = 0; i < maxSpouse; ++i) {
			spouseIdx[i] = infoIdx(resultFields[i]);
			// clear these fields for the last generation
			IndInfoIterator ptr = infoBegin(spouseIdx[i]);
			IndInfoIterator ptrEnd = infoEnd(spouseIdx[i]);
			for ( ; ptr != ptrEnd; ++ptr)
				*ptr = static_cast<InfoType>(-1);
		}
	}
	// find all the couples
	typedef std::pair<ULONG, ULONG> couple;
	vector<couple> couples;
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (size_t i = 0; i < popSize(); ++i) {
			double f = individual(i).info(m_fatherIdx);
			double m = individual(i).info(m_motherIdx);
			if (f >= 0 && m >= 0) {
				if (excludeOutbred) {
					// if they share a father or a mother.
					try {
						double f1 = indByID(m).info(m_fatherIdx);
						double m1 = indByID(m).info(m_motherIdx);
						double f2 = indByID(f).info(m_fatherIdx);
						double m2 = indByID(f).info(m_motherIdx);
						if ((fcmp_ne(f1, -1) && fcmp_eq(f1, f2)) || (fcmp_ne(m1, -1) && fcmp_eq(m1, m2)))
							continue;
					} catch (IndexError &) {
						// if parent not found, does not matter.
						// pass
					}
				}
				couples.push_back(couple(static_cast<ULONG>(f + 0.5), static_cast<ULONG>(m + 0.5)));
			}
		}
	}
	// now, look for each pair and assign spouse
	std::map<ULONG, UINT> numSpouse;
	for (size_t i = 0; i < couples.size(); ++i) {
		// now p and m are spouse
		ULONG p = couples[i].first;
		ULONG m = couples[i].second;
		try {
			// but these guys might not be in the population...
			Individual & fa = indByID(p);
			Individual & ma = indByID(m);
			DBG_ASSERT(ma.sex() == FEMALE && fa.sex() == MALE, RuntimeError,
				"Sex of parents appear to be wrong.");
			if (numSpouse[p] < maxSpouse && sexChoice != MALE_ONLY
			    && acceptableAffectionStatus(ma.affected(), affectionChoice)) {
				bool valid = true;
				// duplicate spouse
				for (size_t s = 0; s < numSpouse[p]; ++s)
					if (fa.info(spouseIdx[s]) == m) {
						valid = false;
						break;
					}
				if (valid) {
					fa.setInfo(m, spouseIdx[numSpouse[p]]);
					++numSpouse[p];
				}
			}
			if (numSpouse[m] < maxSpouse && sexChoice != FEMALE_ONLY
			    && acceptableAffectionStatus(fa.affected(), affectionChoice)) {
				bool valid = true;
				// duplicate spouse
				for (size_t s = 0; s < numSpouse[m]; ++s)
					if (individual(m).info(spouseIdx[s]) == p) {
						valid = false;
						break;
					}
				if (valid) {
					ma.setInfo(p, spouseIdx[numSpouse[m]]);
					++numSpouse[m];
				}
			}                                                                                               // idx
		} catch (...) {
			// if does not found, ignore this couple.
		}
	}
}


void Pedigree::locateOffspring(SexChoice sexChoice, AffectionStatus affectionChoice, const vectorstr & resultFields, int ancGen)
{
	DBG_ASSERT(resultFields.size() >= 1, ValueError,
		"Please provide at least one information field to store offspring");

	UINT maxOffspring = resultFields.size();

	vectori offspringIdx(maxOffspring);

	// clear offspring fields
	UINT topGen = ancGen == -1 ? ancestralGens() : std::min(ancestralGens(), static_cast<UINT>(ancGen));
	for (unsigned ans = topGen; ans >= 0; --ans) {
		useAncestralGen(ans);
		for (size_t i = 0; i < maxOffspring; ++i) {
			offspringIdx[i] = infoIdx(resultFields[i]);
			// clear these fields for the last generation
			for (IndInfoIterator ptr = infoBegin(offspringIdx[i]);
			     ptr != infoEnd(offspringIdx[i]); ++ptr)
				*ptr = static_cast<InfoType>(-1);
		}
	}

	// use individual ID
	std::map<ULONG, UINT> numOffspring;
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (size_t i = 0; i < popSize(); ++i) {
			// everyone has one or two parents.
			double p = individual(i).info(m_fatherIdx);
			double m = individual(i).info(m_motherIdx);
			if (p < 0 || m < 0)
				continue;
			try {
				// but the parents might not exist
				ULONG pp = static_cast<ULONG>(p);
				ULONG mm = static_cast<ULONG>(m);
				Individual & fa = indByID(pp);
				Individual & ma = indByID(mm);
				// add child as father's offspring
				if (numOffspring[pp] < maxOffspring &&
				    acceptableSex(MALE, individual(i).sex(), sexChoice) &&
				    acceptableAffectionStatus(individual(i).affected(), affectionChoice)) {
					fa.setInfo(individual(i).info(m_idIdx), offspringIdx[numOffspring[pp]]);
					++numOffspring[pp];
				}
				// add child as mother's offspring
				if (numOffspring[mm] < maxOffspring &&
				    acceptableSex(FEMALE, individual(i).sex(), sexChoice) &&
				    acceptableAffectionStatus(individual(i).affected(), affectionChoice)) {
					ma.setInfo(individual(i).info(m_idIdx), offspringIdx[numOffspring[mm]]);
					++numOffspring[mm];
				}
			} catch (...) {
				// does not care if a parent cannot be found
			}
		}
	}

}


void Pedigree::locateSibling(SexChoice sexChoice, AffectionStatus affectionChoice, const vectorstr & resultFields, int ancGen)
{
	DBG_ASSERT(resultFields.size() >= 1, ValueError,
		"Please provide at least one information field to store offspring");

	UINT maxSibling = resultFields.size();
	vectori siblingIdx(maxSibling);

	for (size_t i = 0; i < maxSibling; ++i)
		siblingIdx[i] = infoIdx(resultFields[i]);

	// clear all fields
	UINT topGen = ancGen == -1 ? ancestralGens() : std::min(ancestralGens(), static_cast<UINT>(ancGen));
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		for (IndIterator it = indIterator(); it.valid(); ++it)
			for (size_t i = 0; i < maxSibling; ++i)
				it->setInfo(-1, siblingIdx[i]);
	}

	// find all single families
	map<ULONG, vectoru> families;
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (size_t i = 0; i < popSize(); ++i) {
			double f = individual(i).info(m_fatherIdx);
			double m = individual(i).info(m_motherIdx);
			if (f >= 0)
				families[static_cast<ULONG>(f)].push_back(static_cast<ULONG>(individual(i).info(m_idIdx)));
			if (m >= 0)
				families[static_cast<ULONG>(m)].push_back(static_cast<ULONG>(individual(i).info(m_idIdx)));
		}
	}
	// look in each single-parent family
	std::map<ULONG, UINT> numSibling;
	map<ULONG, vectoru>::iterator it = families.begin();
	map<ULONG, vectoru>::iterator itEnd = families.end();
	for (; it != itEnd; ++it) {
		// these guys share at least one parent!
		const vectoru & sibs = it->second;
		if (sibs.size() <= 1)
			continue;
		for (size_t i = 0; i < sibs.size(); ++i) {
			try {
				Individual & child = indByID(sibs[i]);
				for (size_t j = 0; j < sibs.size(); ++j) {
					if (i == j)
						continue;
					try {
						Individual & sibling = indByID(sibs[j]);

						if (numSibling[sibs[i]] < maxSibling &&
						    acceptableSex(child.sex(), sibling.sex(), sexChoice) &&
						    acceptableAffectionStatus(sibling.affected(), affectionChoice)) {
							bool valid = true;
							// duplicate sibling
							for (size_t s = 0; s < numSibling[sibs[i]]; ++s)
								if (child.info(siblingIdx[s]) == sibs[j]) {
									valid = false;
									break;
								}
							if (valid) {
								child.setInfo(sibs[j], siblingIdx[numSibling[sibs[i]]]);
								++numSibling[sibs[i]];
							}
						}
						if (numSibling[sibs[j]] < maxSibling &&
						    acceptableSex(sibling.sex(), child.sex(), sexChoice) &&
						    acceptableAffectionStatus(child.affected(), affectionChoice)) {
							bool valid = true;
							// duplicate sibling
							for (size_t s = 0; s < numSibling[sibs[j]]; ++s)
								if (sibling.info(siblingIdx[s]) == sibs[i]) {
									valid = false;
									break;
								}
							if (valid) {
								sibling.setInfo(sibs[i], siblingIdx[numSibling[sibs[j]]]);
								++numSibling[sibs[j]];
							}
						}
						// idx
					} catch (...) {
						continue;
					}
				}
			} catch (...) {
				// if does not found, ignore this couple.
			}
		}
	}
}


void Pedigree::locateFullSibling(SexChoice sexChoice, AffectionStatus affectionChoice, const vectorstr & resultFields, int ancGen)
{
	DBG_ASSERT(resultFields.size() >= 1, ValueError,
		"Please provide at least one information field to store offspring");

	DBG_FAILIF(numParents() != 2, ValueError,
		"Please provide two parental information fields");

	UINT maxSibling = resultFields.size();
	vectori siblingIdx(maxSibling);

	for (size_t i = 0; i < maxSibling; ++i)
		siblingIdx[i] = infoIdx(resultFields[i]);

	UINT topGen = ancGen == -1 ? ancestralGens() : std::min(ancestralGens(), static_cast<UINT>(ancGen));
	// clear all fields
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		for (IndIterator it = indIterator(); it.valid(); ++it)
			for (size_t i = 0; i < maxSibling; ++i)
				it->setInfo(-1, siblingIdx[i]);
	}

	// find all full families
	typedef std::pair<ULONG, ULONG> couple;
	map<couple, vectoru> families;
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (size_t i = 0; i < popSize(); ++i) {
			double f = individual(i).info(m_fatherIdx);
			double m = individual(i).info(m_motherIdx);
			if (f >= 0 && m >= 0)
				families[couple(static_cast<ULONG>(f), static_cast<ULONG>(m))].push_back(static_cast<ULONG>(individual(i).info(m_idIdx)));
		}
	}
	// look in each single-parent family
	std::map<ULONG, UINT> numSibling;
	map<couple, vectoru>::iterator it = families.begin();
	map<couple, vectoru>::iterator itEnd = families.end();
	for (; it != itEnd; ++it) {
		// these guys share two parents....
		const vectoru & sibs = it->second;
		if (sibs.size() <= 1)
			continue;
		for (size_t i = 0; i < sibs.size(); ++i) {
			try {
				Individual & child = indByID(sibs[i]);
				for (size_t j = 0; j < sibs.size(); ++j) {
					if (i == j)
						continue;
					Individual & sibling = indByID(sibs[j]);

					if (numSibling[sibs[i]] < maxSibling &&
					    acceptableSex(child.sex(), sibling.sex(), sexChoice) &&
					    acceptableAffectionStatus(sibling.affected(), affectionChoice)) {
						bool valid = true;
						// duplicat sibling
						for (size_t s = 0; s < numSibling[sibs[i]]; ++s)
							if (child.info(siblingIdx[s]) == sibs[j]) {
								valid = false;
								break;
							}
						if (valid) {
							child.setInfo(sibs[j], siblingIdx[numSibling[sibs[i]]]);
							++numSibling[sibs[i]];
						}
					}
					if (numSibling[sibs[j]] < maxSibling &&
					    acceptableSex(sibling.sex(), child.sex(), sexChoice) &&
					    acceptableAffectionStatus(child.affected(), affectionChoice)) {
						bool valid = true;
						// duplicate sibling
						for (size_t s = 0; s < numSibling[sibs[j]]; ++s)
							if (sibling.info(siblingIdx[s]) == sibs[i]) {
								valid = false;
								break;
							}
						if (valid) {
							sibling.setInfo(sibs[i], siblingIdx[numSibling[sibs[j]]]);
							++numSibling[sibs[j]];
						}
					}                                                                                       // idx
				}
			} catch (...) {
				// if does not found, ignore this couple.
			}
		}
	}
}


void Pedigree::locateCommonOffspring(SexChoice sexChoice, AffectionStatus affectionChoice, const vectorstr & resultFields, int ancGen)
{
	DBG_ASSERT(numParents() == 2, ValueError,
		"This relative only exists when there are two parents for each indidivual");

	DBG_ASSERT(resultFields.size() >= 2, ValueError,
		"Please provide at least one information field for spouse and one for offspring.");

	UINT topGen = ancGen == -1 ? ancestralGens() : std::min(ancestralGens(), static_cast<UINT>(ancGen));
	UINT maxOffspring = resultFields.size() - 1;

	vectori offspringIdx(maxOffspring);
	int spouseIdx = infoIdx(resultFields[0]);

	// clear all fields
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (size_t i = 0; i < maxOffspring; ++i) {
			offspringIdx[i] = infoIdx(resultFields[i + 1]);
			// clear these fields for the last generation
			IndInfoIterator ptr = infoBegin(offspringIdx[i]);
			IndInfoIterator ptrEnd = infoEnd(offspringIdx[i]);
			for (; ptr != ptrEnd; ++ptr)
				*ptr = static_cast<InfoType>(-1);
		}
	}

	// find all full families
	typedef std::pair<ULONG, ULONG> couple;
	map<couple, vectoru> families;
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (size_t i = 0; i < popSize(); ++i) {
			double f = individual(i).info(m_fatherIdx);
			double m = individual(i).info(m_motherIdx);
			if (f >= 0 && m >= 0)
				families[couple(static_cast<ULONG>(f), static_cast<ULONG>(m))].push_back(static_cast<ULONG>(individual(i).info(m_idIdx)));
		}
	}
	// look in each family
	map<couple, vectoru>::iterator it = families.begin();
	map<couple, vectoru>::iterator itEnd = families.end();
	std::map<ULONG, UINT> numOffspring;
	for (; it != itEnd; ++it) {
		ULONG p = it->first.first;
		ULONG m = it->first.second;
		// these guys share two parents
		const vectoru & offspring = it->second;
		try {
			Individual & fa = indByID(p);
			Individual & ma = indByID(m);
			// spouse
			double fa_spouse = fa.info(spouseIdx);
			double ma_spouse = ma.info(spouseIdx);
			bool valid_fa = fa_spouse != -1 && static_cast<ULONG>(fa_spouse) == m;
			bool valid_ma = ma_spouse != -1 && static_cast<ULONG>(ma_spouse) == p;
			if (!valid_fa && !valid_ma)
				continue;
			// offspring
			for (size_t j = 0; j < offspring.size(); ++j) {
				try {
					Individual & child = indByID(offspring[j]);
					bool valid = acceptableSex(MALE, child.sex(), sexChoice) &&
					             acceptableAffectionStatus(child.affected(), affectionChoice);
					// duplicate child
					if (valid) {
						for (size_t s = 0; s < numOffspring[p]; ++s)
							if (fa.info(offspringIdx[s]) == offspring[j]) {
								valid = false;
								break;
							}
					}
					if (valid) {
						if (valid_fa && numOffspring[p] < maxOffspring) {
							fa.setInfo(offspring[j], offspringIdx[numOffspring[p]]);
							++numOffspring[p];
						}
						if (valid_ma && numOffspring[m] < maxOffspring) {
							ma.setInfo(offspring[j], offspringIdx[numOffspring[m]]);
							++numOffspring[m];
						}
					}
				} catch (...) {
					// if does not found, ignore this offspring
				}
			}
		} catch (...) {
			// if does not found, ignore this couple.
		}
	}
}


bool Pedigree::traceRelatives(const stringMatrix & fieldPath,
                              const uintList & sexChoiceList,  const uintList & affectionChoiceList,
                              const stringList & resultFieldList, int ancGen)
{
	const matrixstr & pathFields = fieldPath.elems();
	const vectoru & sexChoice = sexChoiceList.elems();
	const vectoru & affectionChoice = affectionChoiceList.elems();
	const vectorstr & resultFields = resultFieldList.elems();

	DBG_FAILIF(!sexChoice.empty() && sexChoice.size() != pathFields.size(),
		ValueError,
		"Parameter sexChoice, if given, should have the same length of pathFields");

	DBG_FAILIF(!affectionChoice.empty() && affectionChoice.size() != pathFields.size(),
		ValueError,
		"Parameter affectionChoice, if given, should have the same length of pathFields");

	UINT topGen = ancGen == -1 ? ancestralGens() : std::min(ancestralGens(), static_cast<UINT>(ancGen));

	vectori resultIdx(resultFields.size());
	for (size_t i = 0; i < resultIdx.size(); ++i)
		resultIdx[i] = infoIdx(resultFields[i]);
	vectoru numResult(popSize(), 0);
	UINT maxResult = resultIdx.size();
	// clear values
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (IndIterator ind = indIterator(); ind.valid(); ++ind)
			for (size_t i = 0; i < maxResult; ++i)
				ind->setInfo(-1, resultIdx[i]);
	}
	// convert pathFields to pathIdx
	intMatrix pathIdx(pathFields.size());
	for (size_t i = 0; i < pathFields.size(); ++i) {
		pathIdx[i] = vectori(pathFields[i].size());
		for (size_t j = 0; j < pathFields[i].size(); ++j)
			pathIdx[i][j] = infoIdx(pathFields[i][j]);
	}
	// convert sexChoice to type SexChoices
	vector<SexChoice> sexes(pathIdx.size(), ANY_SEX);
	for (size_t i = 0; i < sexChoice.size(); ++i) {
		DBG_ASSERT(sexChoice[i] == ANY_SEX || sexChoice[i] == MALE_ONLY || sexChoice[i] == FEMALE_ONLY
			|| sexChoice[i] == SAME_SEX || sexChoice[i] == OPPOSITE_SEX, ValueError,
			"Relative sex can only be one of ANY_SEX, MALE_ONLY, FEMALE_ONLY, SAME_SEX or OPPOSITE_SEX.");
		sexes[i] = static_cast<SexChoice>(sexChoice[i]);
	}
	// convert affectionChoice to type SexChoices
	vector<AffectionStatus> affections(pathIdx.size(), ANY_AFFECTION_STATUS);
	for (size_t i = 0; i < affectionChoice.size(); ++i) {
		DBG_ASSERT(affectionChoice[i] == AFFECTED || affectionChoice[i] == UNAFFECTED
			|| affectionChoice[i] == ANY_AFFECTION_STATUS, ValueError,
			"Relative affection status can only be one of AFFECTED, UNAFFECTED and ANY_AFFECTION_STATUS.");
		affections[i] = static_cast<AffectionStatus>(affectionChoice[i]);
	}

	ULONG idx = 0;
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (IndIterator ind = indIterator(); ind.valid(); ++ind, ++idx) {
			Sex mySex = ind->sex();
			vectoru inds = vectoru(1, static_cast<ULONG>(ind->info(m_idIdx)));
			// go through the path
			for (size_t path = 0; path < pathFields.size(); ++path) {
				DBG_DO(DBG_POPULATION, cerr << "Start of path " << path
					                        << " : " << inds << endl);
				const vectori & fields = pathIdx[path];
				SexChoice sex = sexes[path];
				AffectionStatus affection = affections[path];

				vectoru newInds;
				// for all individuals
				for (size_t i = 0; i < inds.size(); ++i) {
					// for all fields
					for (size_t s = 0; s < fields.size(); ++s) {
						double sID = indByID(inds[i]).info(fields[s]);
						if (sID < 0)
							continue;
						Individual & sind = indByID(sID);
						if (!acceptableSex(mySex, sind.sex(), sex))
							continue;
						if (!acceptableAffectionStatus(sind.affected(), affection))
							continue;
						newInds.push_back(static_cast<ULONG>(sID + 0.5));
					}
				}
				inds.swap(newInds);
				if (inds.empty())
					break;
			}
			DBG_DO(DBG_POPULATION, cerr << "Ind " << idx << " has relatives " << inds << endl);
			// ind has the results
			for (size_t i = 0; i < maxResult; ++i)
				if (i < inds.size())
					ind->setInfo(inds[i], resultIdx[i]);
		}
	}
	useAncestralGen(0);
	return true;
}


vectoru Pedigree::individualsWithRelatives(const stringList & infoFieldList, const uintList & sexChoiceList,
                                           const uintList & affectionChoiceList,
                                           const subPopList & subPops, int ancGen)
{
	const vectoru & sexChoice = sexChoiceList.elems();
	const vectoru & affectionChoice = affectionChoiceList.elems();
	const vectorstr & infoFields = infoFieldList.elems();

	DBG_FAILIF(!sexChoice.empty() && sexChoice.size() != infoFields.size(),
		ValueError,
		"Parameter sexChoice, if given, should have the same length of infoFields");

	DBG_FAILIF(!affectionChoice.empty() && affectionChoice.size() != infoFields.size(),
		ValueError,
		"Parameter affectionChoice, if given, should have the same length of infoFields");

	vectoru fieldIdx(infoFields.size());
	for (size_t i = 0; i < fieldIdx.size(); ++i)
		fieldIdx[i] = infoIdx(infoFields[i]);

	// convert sexChoice to type SexChoices
	vector<SexChoice> sexes(infoFields.size(), ANY_SEX);
	for (size_t i = 0; i < sexChoice.size(); ++i) {
		DBG_ASSERT(sexChoice[i] == ANY_SEX || sexChoice[i] == MALE_ONLY || sexChoice[i] == FEMALE_ONLY
			|| sexChoice[i] == SAME_SEX || sexChoice[i] == OPPOSITE_SEX, ValueError,
			"Relative sex can only be one of ANY_SEX, MALE_ONLY, FEMALE_ONLY, SAME_SEX or OPPOSITE_SEX.");
		sexes[i] = static_cast<SexChoice>(sexChoice[i]);
	}

	// convert affectionChoice to type SexChoices
	vector<AffectionStatus> affections(infoFields.size(), ANY_AFFECTION_STATUS);
	for (size_t i = 0; i < affectionChoice.size(); ++i) {
		DBG_ASSERT(affectionChoice[i] == AFFECTED || affectionChoice[i] == UNAFFECTED
			|| affectionChoice[i] == ANY_AFFECTION_STATUS, ValueError,
			"Relative affection status can only be one of AFFECTED, UNAFFECTED and ANY_AFFECTION_STATUS.");
		affections[i] = static_cast<AffectionStatus>(affectionChoice[i]);
	}

	// mark eligible Individuals
	UINT topGen = ancGen == -1 ? ancestralGens() : std::min(ancestralGens(), static_cast<UINT>(ancGen));
	for (unsigned ans = 0; ans <= ancestralGens(); ++ans) {
		useAncestralGen(ans);
		if (ans > topGen) {
			markIndividuals(vspID(), false);
			continue;
		}
		if (subPops.allAvail())
			markIndividuals(vspID(), true);
		else {
			markIndividuals(vspID(), false);
			subPopList::const_iterator it = subPops.begin();
			subPopList::const_iterator itEnd = subPops.end();
			for (; it != itEnd; ++it)
				markIndividuals(*it, true);
		}
	}

	vectoru IDs;
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (IndIterator ind = indIterator(); ind.valid(); ++ind) {
			if (!ind->marked())
				continue;
			bool valid = true;
			for (size_t i = 0; i < fieldIdx.size(); ++i) {
				double rel = ind->info(fieldIdx[i]);
				if (rel < 0) {
					valid = false;
					break;
				}
				try {
					// valid?
					Individual & rind = indByID(rel);
					if (!rind.marked() ||
					    !acceptableSex(ind->sex(), rind.sex(), sexes[i]) ||
					    !acceptableAffectionStatus(ind->affected(), affections[i])) {
						valid = false;
						break;
					}
				} catch (IndexError &) {
					valid = false;
					break;
				}
			}
			if (valid)
				IDs.push_back(static_cast<ULONG>(ind->info(m_idIdx) + 0.5));
		}
	}
	useAncestralGen(0);
	return IDs;
}


// const unsigned long UnusedIndividual = std::numeric_limits<unsigned long>::max();
//
// Pedigree::Pedigree(int numParents, const string & pedfile)
//  : m_numParents(numParents)
// {
//  DBG_ASSERT(numParents == 1 || numParents == 2, ValueError,
//      "individuals in a pedigree can have one or two parents");
//
//  if (!pedfile.empty())
//      load(pedfile);
// }
////
// vectoru Pedigree::subPopSizes(ULONG gen)
// {
//  CHECK_GEN(gen);
//  return m_pedSize[gen];
// }
//
//
// ULONG Pedigree::subPopSize(ULONG gen, SubPopID subPop)
// {
//  CHECK_GEN(gen);
//  CHECK_SUBPOP(gen, subPop);
//
//  return m_pedSize[gen][subPop];
// }
//
//
// void Pedigree::addGen(const vectoru & sizes)
// {
//  ULONG popSize = accumulate(sizes.begin(), sizes.end(), 0UL);
//
//  m_paternal.push_back(vectoru(popSize));
//  if (m_numParents == 2)
//      m_maternal.push_back(vectoru(popSize));
//  if (m_info.size() > 0) {
//      UINT infoSize = m_infoNames.size();
//      m_info.push_back(vector<vectorf>(popSize));
//      for (size_t i = 0; i < popSize; ++i)
//          m_info.back()[i].resize(infoSize);
//  }
//  m_pedSize.push_back(sizes);
// }
//
//
// void Pedigree::load(const string & filename)
// {
//  m_paternal.clear();
//  m_maternal.clear();
//  m_pedSize.clear();
//
//  ifstream ifs(filename.c_str());
//
//  DBG_FAILIF(!ifs, SystemError, "Can not open pedigree file " + filename + " to read");
//
//  string line;
//  while (getline(ifs, line)) {
//      istringstream input(line);
//      long int idx;
//      vectoru values;
//      vectoru sizes;
//      while (!input.eof()) {
//          input >> idx;
//          values.push_back(idx == -1 ? UnusedIndividual : static_cast<ULONG>(idx));
//          input >> ws;
//          if (input.peek() == '#') {
//              // ignore '#'
//              input.ignore();
//              // start reading population size.
//              while (!input.eof()) {
//                  input >> idx;
//                  sizes.push_back(idx);
//                  input >> ws;
//              }
//              // exit the outer loop
//              break;
//          }
//      }
//      // check if things are all right
//      // determine number of parents...
//      ULONG popSize = accumulate(sizes.begin(), sizes.end(), 0UL);
//      DBG_FAILIF(popSize != 0 && values.empty(),
//          ValueError, "No parent is read");
//      if (popSize == 0 && values.empty())
//          continue;
//      if (m_numParents == 0)
//          m_numParents = values.size() / popSize;
//      DBG_ASSERT(m_numParents * popSize == values.size(), ValueError,
//          "Number of parents does not match subpopulation sizes.\n"
//          "Line: " + toStr(m_paternal.size() + 1) + ", individuals read: "
//          + toStr(values.size()) +
//          ", Pop size: " + toStr(popSize));
//      m_pedSize.push_back(vectoru());
//      m_pedSize.back().swap(sizes);
//      if (m_numParents == 1) {
//          m_paternal.push_back(vectoru());
//          m_paternal.back().swap(values);
//      } else if (m_numParents == 2) {
//          m_paternal.push_back(vectoru(popSize));
//          m_maternal.push_back(vectoru(popSize));
//          for (size_t i = 0; i < popSize; ++i) {
//              m_paternal.back()[i] = values[2 * i];
//              m_maternal.back()[i] = values[2 * i + 1];
//          }
//      } else {
//          DBG_ASSERT(false, SystemError,
//              "Sorry, pedigree does not support more than two parents");
//      }
//  }
//  ifs.close();
// }
//
//
// void Pedigree::loadInfo(const string & filename, const string & name)
// {
//  vectorstr names(1, name);
//
//  loadInfo(filename, names);
// }
//
//
// void Pedigree::loadInfo(const string & filename, const vectorstr & names)
// {
//  ifstream afs(filename.c_str());
//
//  DBG_FAILIF(!afs, SystemError, "Can not open auxiliary information Pedigree" + filename + " to read");
//  size_t gen = 0;
//  size_t numInfo = names.size();
//
//  for (size_t i = 0; i < numInfo; ++i) {
//      DBG_ASSERT(find(m_infoNames.begin(), m_infoNames.end(), names[i]) == m_infoNames.end(),
//          ValueError, "Information " + names[i] + " has already been loaded");
//      m_infoNames.push_back(names[i]);
//  }
//  string line;
//  while (getline(afs, line)) {
//      istringstream input(line);
//      vectorf values;
//      double value;
//      while (!input.eof()) {
//          input >> value;
//          values.push_back(value);
//          input >> ws;
//      }
//      DBG_FAILIF(gen >= m_paternal.size(), ValueError,
//          "Information pedigree is larger than parental Pedigree");
//      ULONG size = m_paternal[gen].size();
//
//      DBG_FAILIF(numInfo * size != values.size(), ValueError,
//          "At generation " + toStr(gen) + ", number of information read is "
//          + toStr(values.size()) + ", which is not a multiple of number of individuals "
//          + toStr(m_paternal[gen].size()));
//      //
//      if (m_info.size() <= gen)
//          m_info.push_back(vector<vectorf>(size));
//      DBG_FAILIF(m_info.size() < gen + 1, ValueError,
//          "Error loading information Pedigree");
//      size_t idx = 0;
//      for (size_t i = 0; i < size; ++i)
//          for (size_t j = 0; j < numInfo; ++j, ++idx)
//              m_info[gen][i].push_back(values[idx]);
//      //
//      gen++;
//  }
//  afs.close();
// }
//
//
// void Pedigree::addInfo(const string & name, double init)
// {
//  DBG_ASSERT(find(m_infoNames.begin(), m_infoNames.end(), name) == m_infoNames.end(),
//      ValueError, "Information " + name + " has already been loaded");
//  m_infoNames.push_back(name);
//
//  for (size_t gen = 0; gen < m_paternal.size(); ++gen) {
//      ULONG size = m_paternal[gen].size();
//      if (m_info.size() <= gen)
//          m_info.push_back(vector<vectorf>(size));
//      DBG_FAILIF(m_info.size() < gen + 1, ValueError,
//          "Error loading information Pedigree");
//      for (size_t i = 0; i < size; ++i)
//          m_info[gen][i].push_back(init);
//  }
// }
//
//
// void Pedigree::save(const string & filename)
// {
//  ofstream ofs(filename.c_str());
//
//  DBG_FAILIF(!ofs, SystemError, "Can not open pedigree file " + filename + " to write.");
//
//  for (size_t gen = 0; gen < m_paternal.size(); ++gen) {
//      size_t sz = m_paternal[gen].size();
//      for (size_t idx = 0; idx < sz; ++idx) {
//          if (m_paternal[gen][idx] == UnusedIndividual)
//              ofs << -1;
//          else
//              ofs << m_paternal[gen][idx];
//          if (m_numParents == 2) {
//              if (m_maternal[gen][idx] == UnusedIndividual)
//                  ofs << "\t-1";
//              else
//                  ofs << '\t' << m_maternal[gen][idx];
//          }
//          ofs << '\t';
//      }
//      ofs << "#\t";
//      for (size_t idx = 0; idx < m_pedSize[gen].size(); ++idx)
//          ofs << m_pedSize[gen][idx] << '\t';
//      ofs << '\n';
//  }
//  ofs.close();
// }
//
//
// void Pedigree::saveInfo(const string & filename, const string & name)
// {
//  vectorstr names(1, name);
//
//  saveInfo(filename, names);
// }
//
//
// void Pedigree::saveInfo(const string & filename, const vectorstr & names)
// {
//  vectoru idx;
//
//  for (size_t i = 0; i < names.size(); ++i) {
//      size_t j = 0;
//      for (; j < m_infoNames.size(); ++j) {
//          if (m_infoNames[j] == names[i]) {
//              idx.push_back(j);
//              break;
//          }
//      }
//      DBG_FAILIF(j == m_infoNames.size(), ValueError, "Invalid information name: " + names[i]);
//  }
//
//  ofstream afs(filename.c_str());
//  DBG_FAILIF(!afs, SystemError, "Can not open information pedigree file " + filename + " to write.");
//
//  for (size_t gen = 0; gen < m_paternal.size(); ++gen) {
//      size_t sz = m_paternal[gen].size();
//      for (size_t i = 0; i < sz; ++i)
//          for (size_t j = 0; j < idx.size(); ++j)
//              afs << m_info[gen][i][j] << '\t';
//      afs << '\n';
//  }
//  afs.close();
// }
//
//
// void Pedigree::selectIndividuals(const vectoru & inds)
// {
//  DBG_FAILIF(m_paternal.empty(), ValueError,
//      "Can not select individuals from an empty Pedigree");
//
//  size_t size = m_paternal.back().size();
//  vector<bool> used(size, false);
//  vectoru::const_iterator it = inds.begin();
//  vectoru::const_iterator it_end = inds.end();
//  for (; it != it_end; ++it) {
//      DBG_FAILIF(*it >= size, IndexError,
//          "Index exceeded the size of the last generation");
//      used[*it] = true;
//  }
//  for (size_t idx = 0; idx < size; ++idx)
//      if (!used[idx]) {
//          m_paternal.back()[idx] = UnusedIndividual;
//          if (m_numParents == 2)
//              m_maternal.back()[idx] = UnusedIndividual;
//      }
// }
//
//
// void Pedigree::markUnrelated()
// {
//  if (m_paternal.size() <= 1)
//      return;
//
//  vector<bool> used;
//  // starting from the last generation, gen=0 etc will be replaced.
//  for (size_t gen = m_paternal.size() - 1; gen > 0; --gen) {
//      used.clear();
//      used.resize(m_paternal[gen - 1].size(), false);
//      for (size_t idx = 0; idx < m_paternal[gen].size(); ++idx)
//          if (m_paternal[gen][idx] != UnusedIndividual)
//              used[m_paternal[gen][idx]] = true;
//      if (m_numParents == 2)
//          for (size_t idx = 0; idx < m_maternal[gen].size(); ++idx)
//              if (m_maternal[gen][idx] != UnusedIndividual)
//                  used[m_maternal[gen][idx]] = true;
//      for (size_t idx = 0; idx < m_paternal[gen - 1].size(); ++idx)
//          if (!used[idx]) {
//              m_paternal[gen - 1][idx] = UnusedIndividual;
//              if (m_numParents == 2)
//                  m_maternal[gen - 1][idx] = UnusedIndividual;
//          }
//  }
// }
//
//
// void Pedigree::removeUnrelated(bool shift_index)
// {
//  if (m_paternal.size() <= 1)
//      return;
//
//  // starting from the last generation, gen=0 etc will be replaced.
//  if (shift_index) {
//      for (size_t gen = 0; gen < m_paternal.size() - 1; ++gen) {
//          vectoru & curGen = m_paternal[gen];
//          vectoru & nextGen = m_paternal[gen + 1];
//          size_t shift = 0;
//          for (size_t idx = 0; idx < curGen.size(); ++idx)
//              if (curGen[idx] == UnusedIndividual) {
//                  // the next generation, with value > idx will be shifted by 1.
//                  for (size_t idx1 = 0; idx1 < nextGen.size(); ++idx1)
//                      if (nextGen[idx1] != UnusedIndividual && nextGen[idx1] + shift > idx)
//                          nextGen[idx1]--;
//                  shift++;
//              }
//      }
//      // maternal
//      if (m_numParents == 2) {
//          for (size_t gen = 0; gen < m_maternal.size() - 1; ++gen) {
//              vectoru & curGen = m_maternal[gen];
//              vectoru & nextGen = m_maternal[gen + 1];
//              size_t shift = 0;
//              for (size_t idx = 0; idx < curGen.size(); ++idx)
//                  if (curGen[idx] == UnusedIndividual) {
//                      // the next generation, with value > idx will be shifted by 1.
//                      for (size_t idx1 = 0; idx1 < nextGen.size(); ++idx1)
//                          if (nextGen[idx1] != UnusedIndividual && nextGen[idx1] + shift > idx)
//                              nextGen[idx1]--;
//                      shift++;
//                  }
//          }
//      }
//  }
//  // adjust m_pedSize
//  for (size_t gen = 0; gen < m_paternal.size(); ++gen) {
//      size_t idx = 0;
//      for (UINT sp = 0; sp < m_pedSize[gen].size(); ++sp) {
//          UINT spSize = m_pedSize[gen][sp];
//          for (size_t i = 0; i < spSize; ++i, ++idx)
//              if (m_paternal[gen][idx] == UnusedIndividual)
//                  m_pedSize[gen][sp]--;
//      }
//  }
//  // remove individuals
//  // new pedigree generation entries
//  vectoru l_pat;
//  vectoru l_mat;
//  vector<vectorf> l_info;
//  for (size_t gen = 0; gen < m_paternal.size(); ++gen) {
//      l_pat.clear();
//      l_mat.clear();
//      l_info.clear();
//      for (size_t idx = 0; idx < m_paternal[gen].size(); ++idx)
//          if (m_paternal[gen][idx] != UnusedIndividual) {
//              l_pat.push_back(m_paternal[gen][idx]);
//              if (m_numParents == 2) {
//                  DBG_ASSERT(m_maternal[gen][idx] != UnusedIndividual,
//                      ValueError, "Inconsistent maternal and matermal Pedigree");
//                  l_mat.push_back(m_maternal[gen][idx]);
//              }
//              if (!m_info.empty())
//                  l_info.push_back(m_info[gen][idx]);
//          }
//      m_paternal[gen].swap(l_pat);
//      if (m_numParents == 2)
//          m_maternal[gen].swap(l_mat);
//      if (!m_info.empty())
//          m_info[gen].swap(l_info);
//  }
// }
//
//
}


