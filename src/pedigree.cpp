/**
 *  $File: Pedigree.cpp $
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

#include "pedigree.h"
#include <fstream>
using std::ifstream;
using std::ofstream;

#include <set>

namespace simuPOP {

Pedigree::Pedigree(const Population & pop, const lociList & loci,
	const stringList & infoFields, const uintList & ancGens, const string & idField,
	const string & fatherField, const string & motherField, bool stealPop)
	: m_idField(idField), m_fatherField(fatherField), m_motherField(motherField),
	m_idIdx(-1), m_fatherIdx(-1), m_motherIdx(-1)
{
	vectorstr extractFields = infoFields.allAvail() ? pop.infoFields() : infoFields.elems();

	if (!m_idField.empty() && find(extractFields.begin(), extractFields.end(), idField) == extractFields.end())
		extractFields.push_back(idField);
	if (!m_fatherField.empty() && find(extractFields.begin(), extractFields.end(), fatherField) == extractFields.end())
		extractFields.push_back(fatherField);
	if (!m_motherField.empty() && find(extractFields.begin(), extractFields.end(), motherField) == extractFields.end())
		extractFields.push_back(motherField);
	//
	if (stealPop) {
		// swap and leaving an empty population
		swap(const_cast<Population &>(pop));
		// Remove information fields.
		vectorstr removedFields;
		for (size_t i = 0; i < infoSize(); ++i)
			if (find(extractFields.begin(), extractFields.end(), infoField(i)) == extractFields.end())
				removedFields.push_back(infoField(i));
		if (!removedFields.empty())
			removeInfoFields(removedFields);
		// Remove loci
		if (!loci.allAvail())
			removeLoci(vectoru(false), loci);
		// remove individuals???
		if (!ancGens.allAvail())
			keepAncestralGens(ancGens);
	} else {
		// leaving population intact
		Population & ped = pop.extract(loci, extractFields, subPopList(), ancGens);
		swap(ped);
	}

	PARAM_FAILIF(m_idField.empty() || !hasInfoField(m_idField), ValueError,
		"A valid ID information field is needed to create a pedigree object.");
	PARAM_FAILIF(!m_fatherField.empty() && !hasInfoField(m_fatherField), ValueError,
		"Invalid father information field " + m_fatherField);
	PARAM_FAILIF(!m_motherField.empty() && !hasInfoField(m_motherField), ValueError,
		"Invalid mother information field " + m_motherField);

	m_idIdx = static_cast<int>(infoIdx(idField));
	m_fatherIdx = m_fatherField.empty() ? -1 : static_cast<int>(infoIdx(fatherField));
	m_motherIdx = m_motherField.empty() ? -1 : static_cast<int>(infoIdx(motherField));

	buildIDMap();
}


Pedigree::Pedigree(const Pedigree & rhs) :
	Population(rhs),
	m_idField(rhs.m_idField), m_fatherField(rhs.m_fatherField), m_motherField(rhs.m_motherField),
	m_idIdx(rhs.m_idIdx), m_fatherIdx(rhs.m_fatherIdx), m_motherIdx(rhs.m_motherIdx)
{
	// ID map needs to be rebuilt because individuals are copied.
	buildIDMap();
}


void Pedigree::buildIDMap()
{
	// build an ID map
	m_idMap.clear();
	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		for (IndIterator it = indIterator(); it.valid(); ++it) {
			size_t id = toID(it->info(m_idIdx));
			DBG_WARNIF(m_idMap.find(id) != m_idMap.end() && *m_idMap[id] != *it,
				(boost::format("Different individuals share the same ID %1%"
					           " so only the latest Individual will be used. If this is an "
					           "age-structured population, you may want to remove parental generations.") % id).str());
			m_idMap[id] = &*it;
		}
	}
}


Pedigree * Pedigree::clone() const
{
	return new Pedigree(*this);
}


void Pedigree::save(const string & filename, const stringList & fieldList,
                    const lociList & lociList) const
{
	ofstream file(filename.c_str());

	if (!file)
		throw RuntimeError("Cannot open file " + filename + " for write.");

	vectorstr fields = fieldList.allAvail() ? infoFields() : fieldList.elems();
	vectoru indexes;
	for (size_t i = 0; i < fields.size(); ++i)
		indexes.push_back(infoIdx(fields[i]));

	// out << .... is very slow compared to the sprintf implementation.
	//
	// three numbers (maximum 20 charameters) + M F, the buffer should be long enough
	char buffer[96];

	size_t ply = ploidy();
	const vectoru & loci = lociList.elems(this);

	size_t nParents = numParents();
	size_t curGen = curAncestralGen();
	for (int gen = ancestralGens(); gen >= 0; --gen) {
		const_cast<Pedigree *>(this)->useAncestralGen(gen);
		ConstRawIndIterator it = rawIndBegin();
		ConstRawIndIterator it_end = rawIndEnd();
		for (; it != it_end; ++it) {
			size_t myID = toID(it->info(m_idIdx));
			size_t fatherID = 0;
			size_t motherID = 0;
			if (m_fatherIdx != -1) {
				fatherID = toID(it->info(m_fatherIdx));
				if (fatherID && m_idMap.find(fatherID) == m_idMap.end())
					fatherID = 0;
			}
			if (m_motherIdx != -1) {
				motherID = toID(it->info(m_motherIdx));
				if (motherID && m_idMap.find(motherID) == m_idMap.end())
					motherID = 0;
			}
			char sexChar = it->sex() == MALE ? 'M' : 'F';
			char affChar = it->affected() ? 'A' : 'U';
			if (nParents == 0)
				sprintf(buffer, SIZE_T_FORMAT " %c %c", myID, sexChar, affChar);
			else if (nParents == 1)
				sprintf(buffer, SIZE_T_FORMAT " " SIZE_T_FORMAT " %c %c", myID, fatherID ? fatherID : motherID,
					sexChar, affChar);
			else
				sprintf(buffer, SIZE_T_FORMAT " " SIZE_T_FORMAT " " SIZE_T_FORMAT " %c %c", myID, fatherID, motherID,
					sexChar, affChar);
			file << buffer;
			for (size_t i = 0; i < indexes.size(); ++i)
				file << " " << it->info(indexes[i]);
			for (size_t i = 0; i < loci.size(); ++i)
				for (size_t p = 0; p < ply; ++p)
					file << " " << it->allele(loci[i], p);
			file << '\n';
		}
	}
	const_cast<Pedigree *>(this)->useAncestralGen(curGen);
	file.close();
}


size_t Pedigree::numParents() const
{
	return static_cast<size_t>(m_fatherIdx != -1) + static_cast<size_t>(m_motherIdx != -1);
}


Individual & Pedigree::indByID(double fid) const
{
	// essentially m_idMap(toID(fid))
	//
	size_t id = toID(fid);

	DBG_FAILIF(fabs(fid - id) > 1e-8, ValueError,
		"individual ID has to be integer (or a double round to full iteger).");

	IdMap::iterator it = m_idMap.find(id);
	// if still cannot be found, raise an IndexError.
	if (it == m_idMap.end())
		throw IndexError((boost::format("No individual with ID %1% could be found.") % id).str());

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
                               SexChoice sexChoice, AffectionStatus affectionChoice,
                               const uintList & ancGens)
{
	DBG_ASSERT(sexChoice == ANY_SEX || sexChoice == MALE_ONLY || sexChoice == FEMALE_ONLY
		|| sexChoice == SAME_SEX || sexChoice == OPPOSITE_SEX, ValueError,
		"Relative sex can only be one of ANY_SEX, MALE_ONLY, FEMALE_ONLY, SAME_SEX or OPPOSITE_SEX.");

	DBG_ASSERT(affectionChoice == AFFECTED || affectionChoice == UNAFFECTED
		|| affectionChoice == ANY_AFFECTION_STATUS, ValueError,
		"Relative affection status can only be one of AFFECTED, UNAFFECTED and ANY_AFFECTION_STATUS.");

	size_t oldGen = curAncestralGen();
	vectoru gens = ancGens.elems();
	if (ancGens.allAvail())
		for (int gen = 0; gen <= ancestralGens(); ++gen)
			gens.push_back(gen);
	else if (ancGens.unspecified())
		gens.push_back(curAncestralGen());

	switch (relType) {
	case SPOUSE:
		locateSpouse(sexChoice, affectionChoice, resultFields, gens, false);
		break;
	case OUTBRED_SPOUSE:
		locateSpouse(sexChoice, affectionChoice, resultFields, gens, true);
		break;
	case OFFSPRING:
		locateOffspring(sexChoice, affectionChoice, resultFields, gens);
		break;
	case SIBLING:
		locateSibling(sexChoice, affectionChoice, resultFields, gens);
		break;
	case FULLSIBLING:
		locateFullSibling(sexChoice, affectionChoice, resultFields, gens);
		break;
	case COMMON_OFFSPRING:
		locateCommonOffspring(sexChoice, affectionChoice, resultFields, gens);
		break;
	default:
		throw ValueError("Unrecognized relative type");
	}
	useAncestralGen(oldGen);
}


void Pedigree::locateSpouse(SexChoice sexChoice, AffectionStatus affectionChoice, const vectorstr & resultFields,
                            const vectoru & ancGens, bool excludeOutbred)
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
	size_t maxSpouse = resultFields.size();

	vectoru spouseIdx(maxSpouse);

	// clear all fields
	for (unsigned genIdx = 0; genIdx < ancGens.size(); ++genIdx) {
		useAncestralGen(ancGens[genIdx]);
		for (size_t i = 0; i < maxSpouse; ++i) {
			spouseIdx[i] = infoIdx(resultFields[i]);
			// clear these fields for the last generation
			IndInfoIterator ptr = infoBegin(spouseIdx[i]);
			IndInfoIterator ptrEnd = infoEnd(spouseIdx[i]);
			for ( ; ptr != ptrEnd; ++ptr)
				*ptr = -1;
		}
	}
	// find all the couples
	typedef pairu couple;
	vector<couple> couples;
	for (unsigned genIdx = 0; genIdx < ancGens.size(); ++genIdx) {
		useAncestralGen(ancGens[genIdx]);
		for (size_t i = 0; i < popSize(); ++i) {
			double f = individual(i).info(m_fatherIdx);
			double m = individual(i).info(m_motherIdx);
			if (f > 0 && m > 0) {
				if (excludeOutbred) {
					// if they share a father or a mother.
					try {
						double f1 = indByID(m).info(m_fatherIdx);
						double m1 = indByID(m).info(m_motherIdx);
						double f2 = indByID(f).info(m_fatherIdx);
						double m2 = indByID(f).info(m_motherIdx);
						if ((fcmp_ge(f1, 1) && fcmp_eq(f1, f2)) || (fcmp_ge(m1, 1) && fcmp_eq(m1, m2)))
							continue;
					} catch (IndexError &) {
						// if parent not found, does not matter.
						// pass
					}
				}
				couples.push_back(couple(toID(f), toID(m)));
			}
		}
	}
	// now, look for each pair and assign spouse
	std::map<size_t, size_t> numSpouse;
	for (size_t i = 0; i < couples.size(); ++i) {
		// now p and m are spouse
		size_t p = couples[i].first;
		size_t m = couples[i].second;
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
					fa.setInfo(static_cast<double>(m), spouseIdx[numSpouse[p]]);
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
					ma.setInfo(static_cast<double>(p), spouseIdx[numSpouse[m]]);
					++numSpouse[m];
				}
			}                                                                                               // idx
		} catch (...) {
			// if does not found, ignore this couple.
		}
	}
}


void Pedigree::locateOffspring(SexChoice sexChoice, AffectionStatus affectionChoice, const vectorstr & resultFields, const vectoru & ancGens)
{
	DBG_ASSERT(resultFields.size() >= 1, ValueError,
		"Please provide at least one information field to store offspring");

	size_t maxOffspring = resultFields.size();

	vectoru offspringIdx(maxOffspring);

	// clear offspring fields
	for (unsigned genIdx = 0; genIdx < ancGens.size(); ++genIdx) {
		useAncestralGen(ancGens[genIdx]);
		for (size_t i = 0; i < maxOffspring; ++i) {
			offspringIdx[i] = infoIdx(resultFields[i]);
			// clear these fields for the last generation
			for (IndInfoIterator ptr = infoBegin(offspringIdx[i]);
			     ptr != infoEnd(offspringIdx[i]); ++ptr)
				*ptr = static_cast<double>(-1);
		}
	}

	// use individual ID
	std::map<size_t, size_t> numOffspring;
	for (unsigned genIdx = 0; genIdx < ancGens.size(); ++genIdx) {
		useAncestralGen(ancGens[genIdx]);
		for (size_t i = 0; i < popSize(); ++i) {
			// everyone has one or two parents.
			double p = individual(i).info(m_fatherIdx);
			double m = individual(i).info(m_motherIdx);
			if (p < 0 || m < 0)
				continue;
			try {
				// but the parents might not exist
				size_t pp = toID(p);
				size_t mm = toID(m);
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


void Pedigree::locateSibling(SexChoice sexChoice, AffectionStatus affectionChoice, const vectorstr & resultFields, const vectoru & ancGens)
{
	DBG_ASSERT(resultFields.size() >= 1, ValueError,
		"Please provide at least one information field to store offspring");

	size_t maxSibling = resultFields.size();
	vectoru siblingIdx(maxSibling);

	for (size_t i = 0; i < maxSibling; ++i)
		siblingIdx[i] = infoIdx(resultFields[i]);

	// clear all fields
	for (unsigned genIdx = 0; genIdx < ancGens.size(); ++genIdx) {
		useAncestralGen(ancGens[genIdx]);
		for (IndIterator it = indIterator(); it.valid(); ++it)
			for (size_t i = 0; i < maxSibling; ++i)
				it->setInfo(-1, siblingIdx[i]);
	}

	// find all single families
	map<size_t, vectoru> families;
	for (unsigned genIdx = 0; genIdx < ancGens.size(); ++genIdx) {
		useAncestralGen(ancGens[genIdx]);
		for (size_t i = 0; i < popSize(); ++i) {
			double f = individual(i).info(m_fatherIdx);
			double m = individual(i).info(m_motherIdx);
			if (f >= 1)
				families[toID(f)].push_back(toID(individual(i).info(m_idIdx)));
			if (m >= 1)
				families[toID(m)].push_back(toID(individual(i).info(m_idIdx)));
		}
	}
	// look in each single-parent family
	std::map<size_t, size_t> numSibling;
	map<size_t, vectoru>::iterator it = families.begin();
	map<size_t, vectoru>::iterator itEnd = families.end();
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
								child.setInfo(static_cast<double>(sibs[j]), siblingIdx[numSibling[sibs[i]]]);
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
								sibling.setInfo(static_cast<double>(sibs[i]), siblingIdx[numSibling[sibs[j]]]);
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


void Pedigree::locateFullSibling(SexChoice sexChoice, AffectionStatus affectionChoice, const vectorstr & resultFields, const vectoru & ancGens)
{
	DBG_ASSERT(resultFields.size() >= 1, ValueError,
		"Please provide at least one information field to store offspring");

	DBG_FAILIF(numParents() != 2, ValueError,
		"Please provide two parental information fields");

	size_t maxSibling = resultFields.size();
	vectoru siblingIdx(maxSibling);

	for (size_t i = 0; i < maxSibling; ++i)
		siblingIdx[i] = infoIdx(resultFields[i]);

	for (unsigned genIdx = 0; genIdx < ancGens.size(); ++genIdx) {
		useAncestralGen(ancGens[genIdx]);
		// clear all fields
		for (IndIterator it = indIterator(); it.valid(); ++it)
			for (size_t i = 0; i < maxSibling; ++i)
				it->setInfo(-1, siblingIdx[i]);
	}

	// find all full families
	typedef pairu couple;
	map<couple, vectoru> families;
	for (unsigned genIdx = 0; genIdx < ancGens.size(); ++genIdx) {
		useAncestralGen(ancGens[genIdx]);
		for (size_t i = 0; i < popSize(); ++i) {
			double f = individual(i).info(m_fatherIdx);
			double m = individual(i).info(m_motherIdx);
			if (f >= 1 && m >= 1)
				families[couple(toID(f), toID(m))].push_back(toID(individual(i).info(m_idIdx)));
		}
	}
	// look in each single-parent family
	std::map<size_t, size_t> numSibling;
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
							child.setInfo(static_cast<double>(sibs[j]), siblingIdx[numSibling[sibs[i]]]);
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
							sibling.setInfo(static_cast<double>(sibs[i]), siblingIdx[numSibling[sibs[j]]]);
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


void Pedigree::locateCommonOffspring(SexChoice sexChoice, AffectionStatus affectionChoice,
                                     const vectorstr & resultFields, const vectoru & ancGens)
{
	DBG_ASSERT(numParents() == 2, ValueError,
		"This relative only exists when there are two parents for each indidivual");

	DBG_ASSERT(resultFields.size() >= 2, ValueError,
		"Please provide at least one information field for spouse and one for offspring.");

	size_t maxOffspring = resultFields.size() - 1;

	vectoru offspringIdx(maxOffspring);
	size_t spouseIdx = infoIdx(resultFields[0]);

	// clear all fields
	for (unsigned genIdx = 0; genIdx < ancGens.size(); ++genIdx) {
		useAncestralGen(ancGens[genIdx]);
		for (size_t i = 0; i < maxOffspring; ++i) {
			offspringIdx[i] = infoIdx(resultFields[i + 1]);
			// clear these fields for the last generation
			IndInfoIterator ptr = infoBegin(offspringIdx[i]);
			IndInfoIterator ptrEnd = infoEnd(offspringIdx[i]);
			for (; ptr != ptrEnd; ++ptr)
				*ptr = static_cast<double>(-1);
		}
	}

	// find all full families
	typedef pairu couple;
	map<couple, vectoru> families;
	for (unsigned genIdx = 0; genIdx < ancGens.size(); ++genIdx) {
		useAncestralGen(ancGens[genIdx]);
		for (size_t i = 0; i < popSize(); ++i) {
			double f = individual(i).info(m_fatherIdx);
			double m = individual(i).info(m_motherIdx);
			if (f >= 1 && m >= 1)
				families[couple(toID(f), toID(m))].push_back(toID(individual(i).info(m_idIdx)));
		}
	}
	// look in each family
	map<couple, vectoru>::iterator it = families.begin();
	map<couple, vectoru>::iterator itEnd = families.end();
	std::map<size_t, size_t> numOffspring;
	for (; it != itEnd; ++it) {
		size_t p = it->first.first;
		size_t m = it->first.second;
		// these guys share two parents
		const vectoru & offspring = it->second;
		try {
			Individual & fa = indByID(p);
			Individual & ma = indByID(m);
			// spouse
			double fa_spouse = fa.info(spouseIdx);
			double ma_spouse = ma.info(spouseIdx);
			bool valid_fa = fa_spouse != -1 && toID(fa_spouse) == m;
			bool valid_ma = ma_spouse != -1 && toID(ma_spouse) == p;
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
							fa.setInfo(static_cast<double>(offspring[j]), offspringIdx[numOffspring[p]]);
							++numOffspring[p];
						}
						if (valid_ma && numOffspring[m] < maxOffspring) {
							ma.setInfo(static_cast<double>(offspring[j]), offspringIdx[numOffspring[m]]);
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
                              const stringList & resultFieldList, const uintList & ancGens)
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

	size_t oldGen = curAncestralGen();

	vectoru gens = ancGens.elems();
	if (ancGens.allAvail())
		for (int gen = 0; gen <= ancestralGens(); ++gen)
			gens.push_back(gen);
	else if (ancGens.unspecified())
		gens.push_back(curAncestralGen());

	vectoru resultIdx(resultFields.size());
	for (size_t i = 0; i < resultIdx.size(); ++i)
		resultIdx[i] = infoIdx(resultFields[i]);
	vectoru numResult(popSize(), 0);
	size_t maxResult = resultIdx.size();
	// clear values
	for (unsigned genIdx = 0; genIdx < gens.size(); ++genIdx) {
		useAncestralGen(gens[genIdx]);
		for (IndIterator ind = indIterator(); ind.valid(); ++ind)
			for (size_t i = 0; i < maxResult; ++i)
				ind->setInfo(-1, resultIdx[i]);
	}
	// convert pathFields to pathIdx
	matrixi pathIdx(pathFields.size());
	for (size_t i = 0; i < pathFields.size(); ++i) {
		pathIdx[i] = vectori(pathFields[i].size());
		for (size_t j = 0; j < pathFields[i].size(); ++j)
			pathIdx[i][j] = static_cast<int>(infoIdx(pathFields[i][j]));
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

	size_t idx = 0;
	for (unsigned genIdx = 0; genIdx < gens.size(); ++genIdx) {
		useAncestralGen(gens[genIdx]);
		for (IndIterator ind = indIterator(); ind.valid(); ++ind, ++idx) {
			Sex mySex = ind->sex();
			vectoru inds = vectoru(1, toID(ind->info(m_idIdx)));
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
						newInds.push_back(toID(sID));
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
					ind->setInfo(static_cast<double>(inds[i]), resultIdx[i]);
		}
	}
	useAncestralGen(oldGen);
	return true;
}


vectoru Pedigree::individualsWithRelatives(const stringList & infoFieldList, const uintList & sexChoiceList,
                                           const uintList & affectionChoiceList,
                                           const subPopList & subPops, const uintList & ancGens)
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

	vectoru gens = ancGens.elems();
	if (ancGens.allAvail())
		for (int gen = 0; gen <= ancestralGens(); ++gen)
			gens.push_back(gen);
	else if (ancGens.unspecified())
		gens.push_back(curAncestralGen());

	size_t oldGen = curAncestralGen();
	// mark eligible Individuals
	for (int ans = 0; ans <= ancestralGens(); ++ans) {
		useAncestralGen(ans);
		if (std::find(gens.begin(), gens.end(), ans) == gens.end()) {
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
	for (int ans = 0; ans <= ancestralGens(); ++ans) {
		if (std::find(gens.begin(), gens.end(), ans) == gens.end())
			continue;
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
				IDs.push_back(toID(ind->info(m_idIdx)));
		}
	}
	useAncestralGen(oldGen);
	return IDs;
}


vectoru Pedigree::identifyFamilies(const string & pedField, const subPopList & subPops,
                                   const uintList & ancGens)
{
	// step 1: Mark eligible individuals and collect IDs
	std::map<size_t, ssize_t> famID;

	vectoru gens = ancGens.elems();
	if (ancGens.allAvail())
		for (int gen = 0; gen <= ancestralGens(); ++gen)
			gens.push_back(gen);
	else if (ancGens.unspecified())
		gens.push_back(curAncestralGen());

	size_t oldGen = curAncestralGen();
	// mark eligible Individuals
	for (int ans = 0; ans <= ancestralGens(); ++ans) {
		useAncestralGen(ans);
		if (std::find(gens.begin(), gens.end(), static_cast<size_t>(ans)) == gens.end()) {
			markIndividuals(vspID(), false);
			continue;
		}
		if (subPops.allAvail())
			markIndividuals(vspID(), true);
		else {
			markIndividuals(vspID(), false);
			subPopList::const_iterator sp = subPops.begin();
			subPopList::const_iterator spEnd = subPops.end();
			for (; sp != spEnd; ++sp)
				markIndividuals(*sp, true);
		}
		// collect ID.
		RawIndIterator it = rawIndBegin();
		RawIndIterator itEnd = rawIndEnd();
		for (; it != itEnd; ++it)
			if (it->marked())
				famID[toID(it->info(m_idIdx))] = -1;
	}
	// step 2: decide family ID
	size_t famCount = 0;
	//
	std::map<size_t, ssize_t>::iterator it = famID.begin();
	std::map<size_t, ssize_t>::iterator it_end = famID.end();
	for (; it != it_end; ++it) {
		// CASE ONE: if this guy is someone's parent, and has already been processed
		if (it->second >= 0)
			continue;
		// this guy should exist
		Individual * ind = m_idMap[it->first];
		Individual * dad = NULL;
		Individual * mom = NULL;
		ssize_t dadFam = -2;
		ssize_t momFam = -2;
		// try to identify father and mother....
		if (m_fatherIdx != -1) {
			size_t dad_id = toID(ind->info(m_fatherIdx));
			// ok father
			std::map<size_t, ssize_t>::iterator dad_fam = famID.find(dad_id);
			// because father exists in famID
			if (dad_fam != famID.end()) {
				dadFam = dad_fam->second;
				dad = m_idMap[dad_id];
			}
		}
		if (m_motherIdx != -1) {
			size_t mom_id = toID(ind->info(m_motherIdx));
			// ok mother
			std::map<size_t, ssize_t>::iterator mom_fam = famID.find(mom_id);
			// because father exists in famID
			if (mom_fam != famID.end()) {
				momFam = mom_fam->second;
				mom = m_idMap[mom_id];
			}
		}
		// CASE TWO: no parent
		if (dad == NULL && mom == NULL) {
			it->second = famCount++;
		} else if (dad != NULL && mom == NULL) {
			// CASE THREE: One father
			if (dadFam >= 0)
				// follow dad fam
				it->second = dadFam;
			else {
				it->second = famCount;
				famID[toID(dad->info(m_idIdx))] = famCount;
				++famCount;
			}
		} else if (dad == NULL && mom != NULL) {
			// CASE FOUR: One mother
			if (momFam >= 0)
				// follow mom fam
				it->second = momFam;
			else {
				it->second = famCount;
				famID[toID(mom->info(m_idIdx))] = famCount;
				++famCount;
			}
		} else if (dadFam < 0 && momFam < 0) {
			// CASE FIVE: fresh father and mother
			it->second = famCount;
			famID[toID(mom->info(m_idIdx))] = famCount;
			famID[toID(dad->info(m_idIdx))] = famCount;
			++famCount;
		} else if (dadFam >= 0 && momFam < 0) {
			// CASE SIX: fresh mother
			it->second = dadFam;
			famID[toID(mom->info(m_idIdx))] = dadFam;
		} else if (dadFam < 0 && momFam >= 0) {
			// CASE SEVEN: fresh father
			it->second = momFam;
			famID[toID(dad->info(m_idIdx))] = momFam;
		} else if (dadFam == momFam) {
			// CASE EIGHT: a sibling?
			it->second = momFam;
		} else {
			// CASE NINE: we have a conflict.
			ssize_t oldID = std::max(dadFam, momFam);
			ssize_t newID = std::min(dadFam, momFam);
			// everyone with the large ID need to be converted to small ID.
			it->second = newID;
			std::map<size_t, ssize_t>::iterator iit = famID.begin();
			std::map<size_t, ssize_t>::iterator iit_end = famID.end();
			for (; iit != iit_end; ++iit) {
				if (iit->second == oldID)
					iit->second = newID;
				if (iit->second > oldID)
					--iit->second;
			}
			--famCount;
		}
	}
	int pedIdx = pedField.empty() ? -1 : static_cast<int>(infoIdx(pedField));
	// return result
	vectoru famSize(famCount, 0);
	it = famID.begin();
	for (; it != it_end; ++it) {
		ssize_t famID = it->second;
		++famSize[famID];
		if (pedIdx >= 0)
			m_idMap[it->first]->setInfo(static_cast<double>(famID), static_cast<size_t>(pedIdx));
	}
	useAncestralGen(oldGen);
	return famSize;
}


vectoru Pedigree::identifyAncestors(const uintList & IDs,
                                    const subPopList & subPops,
                                    const uintList & ancGens)
{
	vectoru gens = ancGens.elems();

	if (ancGens.allAvail())
		for (int gen = 0; gen <= ancestralGens(); ++gen)
			gens.push_back(gen);
	else if (ancGens.unspecified())
		gens.push_back(curAncestralGen());

	size_t oldGen = curAncestralGen();
	// mark eligible Individuals
	for (int ans = 0; ans <= ancestralGens(); ++ans) {
		useAncestralGen(ans);
		if (std::find(gens.begin(), gens.end(), static_cast<size_t>(ans)) == gens.end()) {
			markIndividuals(vspID(), false);
			continue;
		}
		if (subPops.allAvail())
			markIndividuals(vspID(), true);
		else {
			markIndividuals(vspID(), false);
			subPopList::const_iterator sp = subPops.begin();
			subPopList::const_iterator spEnd = subPops.end();
			for (; sp != spEnd; ++sp)
				markIndividuals(*sp, true);
		}
	}
	useAncestralGen(oldGen);

	// step 2: source IDs
	vectoru res;
	if (IDs.allAvail()) {
		RawIndIterator it = rawIndBegin();
		RawIndIterator itEnd = rawIndEnd();
		for (; it != itEnd; ++it)
			if (it->marked())
				res.push_back(toID(it->info(m_idIdx)));
	} else {
		const vectoru & inputIDs = IDs.elems();
		res.reserve(inputIDs.size());
		for (size_t i = 0; i < inputIDs.size(); ++i)
			if (m_idMap.find(inputIDs[i]) != m_idMap.end())
				res.push_back(inputIDs[i]);
	}
	// step 3: trace back like a spider
	size_t start = 0;
	while (true) {
		size_t end = res.size();
		if (start == end)
			break;
		for (size_t i = start; i < end; ++i) {
			size_t ID = res[i];
			// true ID starts from 1
			size_t father_ID = 0;
			size_t mother_ID = 0;
			try {
				// if this ID exists
				Individual & ind = indByID(ID);
				if (m_fatherIdx != -1)
					father_ID = toID(ind.info(m_fatherIdx));
				if (m_motherIdx != -1)
					mother_ID = toID(ind.info(m_motherIdx));
			} catch (IndexError &) {
				//
			}
			if (father_ID) {
				try {
					Individual & father = indByID(father_ID);
					if (father.marked()) {
						res.push_back(father_ID);
						// this father is already included
						father.setMarked(false);
					}
				} catch (IndexError &) {
					///
				}
			}
			if (mother_ID) {
				try {
					Individual & mother = indByID(mother_ID);
					if (mother.marked()) {
						res.push_back(mother_ID);
						// this mother is already included
						mother.setMarked(false);
					}
				} catch (IndexError &) {
					///
				}
			}
		}
		// all parents of individuals between start and end has been located
		// start to find the parent of these parents
		start = end;
	}
	return res;
}


vectoru Pedigree::identifyOffspring(const uintList & IDs,
                                    const subPopList & subPops,
                                    const uintList & ancGens)
{
	// record offspring of everyone
	std::map<size_t, vectoru> offspringMap;
	vectoru gens = ancGens.elems();
	if (ancGens.allAvail())
		for (int gen = 0; gen <= ancestralGens(); ++gen)
			gens.push_back(gen);
	else if (ancGens.unspecified())
		gens.push_back(curAncestralGen());

	// mark eligible Individuals
	size_t oldGen = curAncestralGen();
	for (int ans = 0; ans <= ancestralGens(); ++ans) {
		useAncestralGen(ans);
		if (std::find(gens.begin(), gens.end(), static_cast<size_t>(ans)) == gens.end()) {
			markIndividuals(vspID(), false);
			continue;
		}
		if (subPops.allAvail())
			markIndividuals(vspID(), true);
		else {
			markIndividuals(vspID(), false);
			subPopList::const_iterator sp = subPops.begin();
			subPopList::const_iterator spEnd = subPops.end();
			for (; sp != spEnd; ++sp)
				markIndividuals(*sp, true);
		}
		// collect ID.
		RawIndIterator it = rawIndBegin();
		RawIndIterator itEnd = rawIndEnd();
		for (; it != itEnd; ++it) {
			// I am a valid offspring
			if (it->marked()) {
				size_t myID = toID(it->info(m_idIdx));
				size_t fatherID = m_fatherIdx == -1 ? 0 : toID(it->info(m_fatherIdx));
				size_t motherID = m_motherIdx == -1 ? 0 : toID(it->info(m_motherIdx));
				// we do not care if father or mother is valid.
				if (fatherID) {
					if (offspringMap.find(fatherID) == offspringMap.end())
						offspringMap[fatherID] = vectoru(1, myID);
					else
						offspringMap[fatherID].push_back(myID);
				}
				if (motherID) {
					if (offspringMap.find(motherID) == offspringMap.end())
						offspringMap[motherID] = vectoru(1, myID);
					else
						offspringMap[motherID].push_back(myID);
				}
			}
		}
	}
	useAncestralGen(oldGen);

	// step 2: locate all offspring
	vectoru res;
	const vectoru & inputIDs = IDs.elems();
	res.reserve(inputIDs.size());
	for (size_t i = 0; i < inputIDs.size(); ++i)
		if (m_idMap.find(inputIDs[i]) != m_idMap.end())
			res.push_back(inputIDs[i]);
	size_t start = 0;
	while (true) {
		size_t end = res.size();
		if (start == end)
			break;
		for (size_t i = start; i < end; ++i) {
			size_t ID = res[i];
			if (offspringMap.find(ID) == offspringMap.end())
				continue;
			res.insert(res.end(), offspringMap[ID].begin(), offspringMap[ID].end());
		}
		// all offspring of individuals between start and end has been located
		// start to find the offspring of these offspring
		start = end;
	}
	// return a unique list
	std::sort(res.begin(), res.end());
	vectoru::iterator new_end = std::unique(res.begin(), res.end());
	res.erase(new_end, res.end());
	return res;
}


void Pedigree::removeIndividuals(const uintList & indexes,
                                 const floatList & IDs, const string & idField, PyObject * filter)
{
	Population::removeIndividuals(indexes, IDs, idField, filter);
	buildIDMap();
}


void Pedigree::removeSubPops(const subPopList & subPops)
{
	Population::removeSubPops(subPops);
	buildIDMap();
}


void Pedigree::push(Population & pop)
{
	Population::push(pop);
	buildIDMap();
}


void Pedigree::addChrom(const vectorf & lociPos, const vectorstr & lociNames,
                        const string & chromName, const stringMatrix & alleleNames,
                        size_t chromType)
{
	Population::addChrom(lociPos, lociNames, chromName, alleleNames, chromType);
	buildIDMap();
}


void Pedigree::addChromFrom(const Population & pop)
{
	Population::addChromFrom(pop);
	buildIDMap();
}


void Pedigree::addIndFrom(const Population & pop)
{
	Population::addIndFrom(pop);
	buildIDMap();
}


size_t Pedigree::mergeSubPops(const uintList & subPops, const string & name)
{
	size_t res = Population::mergeSubPops(subPops, name);

	buildIDMap();
	return res;
}


void Pedigree::resize(const uintList & sizes, bool propagate)
{
	Population::resize(sizes, propagate);
	buildIDMap();
}


void Pedigree::setSubPopByIndInfo(const string & field)
{
	Population::setSubPopByIndInfo(field);
	buildIDMap();
}


struct IndInfo
{
	Sex sex;
	vectoru parents;
	vectoru offspring;
	bool affectionStatus;
	vectorf fields;
	vectora genotype;
	IndInfo() : sex(MALE), parents(), offspring(), affectionStatus(false), fields(), genotype() {}
	IndInfo(size_t off) : sex(MALE), parents(), offspring(1, off), affectionStatus(false), fields(), genotype() {}
};


Pedigree loadPedigree(const string & file, const string & idField, const string & fatherField,
                      const string & motherField, float ploidy, const uintList & _lociList, const uintList & chromTypes,
                      const floatList & lociPos, const stringList & chromNames, const stringMatrix & alleleNames,
                      const stringList & lociNames, const stringList & subPopNames, const stringList & fieldList)
{
	initClock();
	int pldy = ploidy == HAPLODIPLOID ? 2 : static_cast<int>(ploidy);
	//
	const vectorstr & infoFields = fieldList.elems();

	for (size_t i = 0; i < infoFields.size(); ++i) {
		DBG_FAILIF(infoFields[i] == idField || infoFields[i] == fatherField || infoFields[i] == motherField,
			ValueError, "Parameter infoFields can only specify additional fields other than idField, fatherField and motherField.");
	}
	vectoru loci = _lociList.elems();
	size_t genoCols = accumulate(loci.begin(), loci.end(), size_t(0)) * pldy;

	ifstream input(file.c_str());
	if (!input)
		throw RuntimeError("Cannot open specified pedigree file " + file);
	//
	size_t max_parents = 0;
	string line;
	// individual and their parents
#if TR1_SUPPORT == 0
	typedef std::map<size_t, IndInfo *> IdMap;
#elif TR1_SUPPORT == 1
	typedef std::unordered_map<size_t, IndInfo *> IdMap;
#else
	typedef std::tr1::unordered_map<size_t, IndInfo *> IdMap;
#endif
	IdMap individuals;
	while (getline(input, line)) {
		if (line.empty())
			continue;
		//
		IndInfo * info = NULL;
		size_t myID = 0;
		int part = 0;
		char * p = strtok(const_cast<char *>(line.c_str()), " ");
		// boost::tokenizer is proven to be too slow..... (5.5s vs. 1.5s)
		while (p) {
			char * q = p;
			p = strtok(NULL, " ");
			// collect self ID
			if (part == 0) {
				myID = atoi(q);
				if (individuals.find(myID) != individuals.end())
					throw ValueError((boost::format("Duplicate individual ID %1%") % myID).str());
				info = (individuals.insert(IdMap::value_type(myID, new IndInfo())).first)->second;
				++part;
				continue;
				// parental ID and sex
			} else if (part == 1) {
				if (*q == 'M') {
					info->sex = MALE;
					++part;
					continue;
				} else if (*q == 'F') {
					info->sex = FEMALE;
					++part;
					continue;
				} else {
					size_t id = atoi(q);
					if (id) {
						info->parents.push_back(id);
						IdMap::iterator it = individuals.find(id);
						// this is a parent but we do not know if he or she has parent
						if (it == individuals.end())
							individuals[id] = new IndInfo(myID);
						else
							it->second->offspring.push_back(myID);
					}
					if (info->parents.size() > 2)
						throw ValueError("At most two parental IDs are allowed before sex information");
				}
				// parental affection status, can be ignored
			} else if (part == 2) {
				if (*q == 'A')
					info->affectionStatus = true;
				else if (*q == 'U')
					info->affectionStatus = false;
				else
					++part;
			}

			// information fields, can be ignored
			if (part == 3) {
				if (info->fields.size() == infoFields.size())
					++part;
				else
					info->fields.push_back(atof(q));
			}

			// genotype
			if (part == 4)
				info->genotype.push_back(TO_ALLELE(atoi(q)));
		}
		// if there is no valid input...
		if (part == 0)
			continue;
		if (!info->genotype.empty()) {
			if (loci.empty()) {
				loci.push_back(info->genotype.size() / pldy);
				genoCols = info->genotype.size();
				if (loci.back() * pldy != genoCols)
					throw ValueError("Incorrect number of genotype colmns for a diploid population.");
			} else {
				if (genoCols != info->genotype.size())
					throw ValueError("Inconsistent number of columns of genotypes.");
			}
		}
		//
		if (max_parents < info->parents.size())
			max_parents = info->parents.size();
	}
	input.close();
	elapsedTime("Readfile");
	DBG_DO(DBG_POPULATION, cerr << "Information about " << individuals.size() << " individuals are loaded." << endl);
	// create the top most ancestral generation
	// find parents who do not have parents...
	vectorstr fields(1, idField);
	if (!fatherField.empty() && fields.size() < max_parents + 1)
		fields.push_back(fatherField);
	if (!motherField.empty() && fields.size() < max_parents + 1)
		fields.push_back(motherField);
	// this is the end of parental IDs and start of additional fields
	size_t fieldIndex = fields.size();
	// additional information fields
	fields.insert(fields.end(), infoFields.begin(), infoFields.end());
	DBG_DO(DBG_POPULATION, cerr << "Using information fields " << fields << endl);
	//
	if (individuals.empty()) {
		return Population(0, ploidy, loci, chromTypes, lociPos,
			-1, chromNames, alleleNames, lociNames, subPopNames, fields);
	}
	//
	typedef std::set<size_t> IdSet;
	IdSet parents;
	IdMap::iterator it = individuals.begin();
	IdMap::iterator it_end = individuals.end();
	for (; it != it_end; ++it)
		if (it->second->parents.empty())
			parents.insert(it->first);
	if (parents.empty())
		throw ValueError("No parents in the top-most ancestral generation");
	//
	Population pop(vectoru(1, parents.size()), ploidy, loci, chromTypes, lociPos,
	               -1, chromNames, alleleNames, lociNames, subPopNames, fields);
	// set individual info
	RawIndIterator ind = pop.rawIndBegin();
	RawIndIterator ind_end = pop.rawIndEnd();
	IdSet::iterator pit = parents.begin();
	for (; ind != ind_end; ++ind, ++pit) {
		ind->setInfo(static_cast<double>(*pit), 0);
		const IndInfo * info = individuals.find(*pit)->second;
		ind->setSex(info->sex);
		ind->setAffected(info->affectionStatus);
		for (size_t i = 0; i < infoFields.size() && i < info->fields.size(); ++i)
			ind->setInfo(info->fields[i], i + fieldIndex);
		for (size_t i = 0, k = 0; i < genoCols / pldy; ++i)
			for (int j = 0; j < pldy; ++j, ++k)
				ind->setAllele(info->genotype[k], i, j);
	}
	DBG_DO(DBG_POPULATION, cerr << parents.size() << " individuals are located for the top-most ancestral generation" << endl);
	//
	while (true) {
		if (parents.empty())
			break;
		//
		IdSet offspring;
		pit = parents.begin();
		IdSet::iterator pit_end = parents.end();
		for (; pit != pit_end; ++pit) {
			const IdMap::const_iterator info = individuals.find(*pit);
			const vectoru & off = info->second->offspring;
			for (size_t i = 0; i < off.size(); ++i)
				offspring.insert(off[i]);
		}
		DBG_DO(DBG_POPULATION, cerr << offspring.size() << " individuals are located from "
			                        << individuals.size() << " individuals for an ancestral generation" << endl);

		if (offspring.empty())
			break;

		Population off_pop(vectoru(1, offspring.size()), ploidy, loci, chromTypes, lociPos,
		                   0, chromNames, alleleNames, lociNames, subPopNames, fields);
		// set individual info
		ind = off_pop.rawIndBegin();
		ind_end = off_pop.rawIndEnd();
		pit = offspring.begin();
		for (; ind != ind_end; ++ind, ++pit) {
			ind->setInfo(static_cast<double>(*pit), 0);
			const IndInfo * info = individuals.find(*pit)->second;
			for (size_t i = 0; i < info->parents.size(); ++i)
				ind->setInfo(static_cast<double>(info->parents[i]), 1 + i);
			ind->setSex(info->sex);
			ind->setAffected(info->affectionStatus);
			for (size_t i = 0; i < infoFields.size() && i < info->fields.size(); ++i)
				ind->setInfo(info->fields[i], i + fieldIndex);
			for (size_t i = 0, k = 0; i < genoCols / pldy; ++i)
				for (int j = 0; j < pldy; ++j, ++k)
					ind->setAllele(info->genotype[k], i, j);
		}
		//
		parents.swap(offspring);
		pop.push(off_pop);
	}
	elapsedTime("Generation");
	DBG_DO(DBG_POPULATION, cerr << "A pedigree with " << pop.ancestralGens()
		                        << " ancestral generations are created." << endl);

	it = individuals.begin();
	it_end = individuals.end();
	for (; it != it_end; ++it)
		delete it->second;
	// uintList means ALL_AVAIL
	return Pedigree(pop, lociList(), pop.infoFields(), uintList(),
		idField, max_parents > 0 ? fatherField : string(),
		max_parents > 1 ? motherField : string(), true);
}


}
