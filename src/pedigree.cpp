/**
 *  $File: pedigree.cpp $
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

pedigree::pedigree(const population & pop, const uintList & loci,
	const stringList & infoFields, int ancGen, const string & idField,
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
	population & ped = pop.extract(string(), loci, extractFields, ancGen, NULL);

	swap(ped);

	DBG_FAILIF(!m_idField.empty() && !pop.hasInfoField(m_idField), ValueError,
		"Invalid id information field " + m_idField);
	DBG_FAILIF(!m_fatherField.empty() && !pop.hasInfoField(m_fatherField), ValueError,
		"Invalid father information field " + m_fatherField);
	DBG_FAILIF(!m_motherField.empty() && !pop.hasInfoField(m_motherField), ValueError,
		"Invalid mother information field " + m_motherField);

	m_idIdx = m_idField.empty() ? -1 : static_cast<int>(infoIdx(idField));
	m_fatherIdx = m_fatherField.empty() ? -1 : static_cast<int>(infoIdx(fatherField));
	m_motherIdx = m_motherField.empty() ? -1 : static_cast<int>(infoIdx(motherField));

	if (m_idIdx == -1)
		return;

	// build an ID map
	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		for (IndIterator it = indIterator(); it.valid(); ++it)
			m_idMap[static_cast<ULONG>(it->info(m_idIdx))] = &*it;
	}
}


pedigree::pedigree(const pedigree & rhs) :
	population(rhs),
	m_idField(rhs.m_idField), m_fatherField(rhs.m_fatherField), m_motherField(rhs.m_motherField),
	m_idIdx(rhs.m_idIdx), m_fatherIdx(rhs.m_fatherIdx), m_motherIdx(rhs.m_motherIdx)
{
}


pedigree * pedigree::clone() const
{
	return new pedigree(*this);
}


UINT pedigree::numParents()
{
	return static_cast<UINT>(m_fatherIdx != -1) + static_cast<UINT>(m_motherIdx != -1);
}


int pedigree::father(ULONG idx, SubPopID subPop)
{
	if (m_fatherIdx < 0)
		return -1;

	return ind(idx, subPop).intInfo(m_fatherIdx);
}


int pedigree::mother(ULONG idx, SubPopID subPop)
{
	if (m_motherIdx < 0)
		return -1;

	return ind(idx, subPop).intInfo(m_motherIdx);
}


individual & pedigree::indByID(double fid, int ancGen, const string & idField)
{
	DBG_WARNING(!idField.empty(), "Parameter idField of pedigree.indByID is ignored.");

	// essentially m_idMap(static_cast<ULONG>(fid))
	//
	ULONG id = static_cast<ULONG>(fid + 0.5);

	DBG_FAILIF(fabs(fid - id) > 1e-8, ValueError,
		"Individual ID has to be integer (or a double round to full iteger).");

	std::map<ULONG, individual *>::iterator it = m_idMap.find(id);
	// if still cannot be found, raise an IndexError.
	if (it == m_idMap.end())
		throw IndexError("No individual with ID " + toStr(id) + " could be found.");

	return *it->second;
}


void pedigree::locateRelatives(uintList fullRelType, const vectorstr & relFields, int ancGen)
{
	const vectoru & fullType = fullRelType.elems();

	if (fullType.empty() || relFields.empty())
		return;

	RelativeType relType = static_cast<RelativeType>(fullType[0]);

	DBG_FAILIF(fullType.size() > 2, ValueError, "Unrecognized relative type.");

	SexChoice relSex = fullType.size() == 2 ? static_cast<SexChoice>(fullType[1]) : AnySex;

	DBG_ASSERT(relSex == AnySex || relSex == MaleOnly || relSex == FemaleOnly
		|| relSex == SameSex || relSex == OppositeSex, ValueError,
		"Relative sex can only be one of AnySex, MaleOnly, FemaleOnly, SameSex or OppositeSex.");

	UINT topGen = ancGen == -1 ? ancestralGens() : std::min(ancestralGens(), static_cast<UINT>(ancGen));

	UINT oldGen = curAncestralGen();
	if (relType == Self) {
		DBG_ASSERT(relFields.size() == 1, ValueError,
			"Please provide one information field to store Self individuals");
		if (m_idIdx == -1)
			locateSelfByIdx(relSex, relFields, topGen);
		else
			locateSelfByID(relSex, relFields, topGen);
	} else if (relType == Spouse) {
		DBG_ASSERT(numParents() == 2, ValueError,
			"This relative only exists when there are two parents for each indidivual");

		DBG_ASSERT(relFields.size() >= 1, ValueError,
			"Please provide at least one information field to store Self individuals");

		DBG_FAILIF(relSex == SameSex, ValueError, "Can not locate spouses with the same sex");

		if (m_idIdx == -1)
			locateSpouseByIdx(relSex, relFields, topGen);
		else
			locateSpouseByID(relSex, relFields, topGen);
	} else if (relType == Offspring) {
		DBG_ASSERT(relFields.size() >= 1, ValueError,
			"Please provide at least one information field to store offspring");

		if (m_idIdx == -1)
			locateOffspringByIdx(relSex, relFields, topGen);
		else
			locateOffspringByID(relSex, relFields, topGen);
	} else if (relType == Sibling) {
		DBG_ASSERT(relFields.size() >= 1, ValueError,
			"Please provide at least one information field to store offspring");

		if (m_idIdx == -1)
			locateSiblingByIdx(relSex, relFields, topGen);
		else
			locateSiblingByID(relSex, relFields, topGen);
	} else if (relType == FullSibling) {
		DBG_ASSERT(relFields.size() >= 1, ValueError,
			"Please provide at least one information field to store offspring");

		DBG_FAILIF(numParents() != 2, ValueError,
			"Please provide two parental information fields");

		if (m_idIdx == -1)
			locateFullSiblingByIdx(relSex, relFields, topGen);
		else
			locateFullSiblingByID(relSex, relFields, topGen);
	} else if (relType == SpouseAndOffspring) {
		DBG_ASSERT(numParents() == 2, ValueError,
			"This relative only exists when there are two parents for each indidivual");

		DBG_ASSERT(relFields.size() >= 2, ValueError,
			"Please provide at least one information field for spouse and one for offspring.");

		DBG_WARNING(topGen == 0, "Spouse can not be located because there is no parental generation.");

		if (m_idIdx == -1)
			locateSpouseAndOffspringByIdx(relSex, relFields, topGen);
		else
			locateSpouseAndOffspringByID(relSex, relFields, topGen);
	} else {
		throw ValueError("Unrecognized relative type");
	}
	useAncestralGen(oldGen);
}


void pedigree::locateSelfByIdx(SexChoice relSex, const vectorstr & relFields, UINT topGen)
{
	UINT fieldIdx = infoIdx(relFields[0]);

	for (size_t ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (size_t idx = 0; idx < popSize(); ++idx)
			ind(idx).setInfo(idx, fieldIdx);
	}
}


void pedigree::locateSelfByID(SexChoice relSex, const vectorstr & relFields, UINT topGen)
{
	UINT fieldIdx = infoIdx(relFields[0]);

	// copy individual ID from one field to another...
	for (size_t ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (size_t idx = 0; idx < popSize(); ++idx)
			ind(idx).setInfo(ind(idx).info(m_idIdx), fieldIdx);
	}
}


void pedigree::locateSpouseByIdx(SexChoice relSex, const vectorstr & relFields, UINT topGen)
{
	DBG_WARNING(topGen == 0, "Spouse can not be located because there is no parental generation.");

	// clear all...
	UINT maxSpouse = relFields.size();

	vectori spouseIdx(maxSpouse);
	// clear spouse for the last generation.
	useAncestralGen(0);
	for (size_t i = 0; i < maxSpouse; ++i) {
		spouseIdx[i] = infoIdx(relFields[i]);
		// clear these fields for the last generation
		IndInfoIterator ptr = infoBegin(spouseIdx[i]);
		IndInfoIterator ptrEnd = infoEnd(spouseIdx[i]);
		for ( ; ptr != ptrEnd; ++ptr)
			*ptr = static_cast<InfoType>(-1);
	}

	// start from the parental generation
	for (unsigned ans = 1; ans <= topGen; ++ans) {
		vectoru numSpouse;
		// go to offspring generation
		useAncestralGen(ans - 1);
		vectorf father = indInfo(m_fatherIdx);
		vectorf mother = indInfo(m_motherIdx);
		//
		useAncestralGen(ans);
		if (numSpouse.empty())
			numSpouse.resize(popSize(), 0);
		std::fill(numSpouse.begin(), numSpouse.end(), 0);
		//
		for (size_t idx = 0; idx < father.size(); ++idx) {
			DBG_FAILIF(fcmp_eq(father[idx], -1) || fcmp_eq(mother[idx], -1), ValueError,
				"Invalid parental index (-1)");
			ULONG p = static_cast<ULONG>(father[idx]);
			ULONG m = static_cast<ULONG>(mother[idx]);
			DBG_ASSERT(p < popSize() && m < popSize(), IndexError,
				"Parental index out of range of 0 ~ " + toStr(popSize() - 1));
			if (numSpouse[p] < maxSpouse) {
				bool valid = true;
				// if sex is not interested
				if ((ind(m).sex() == Male && relSex == FemaleOnly) ||
				    (ind(m).sex() == Female && relSex == MaleOnly))
					valid = false;
				// duplicate spouse
				if (valid) {
					for (size_t s = 0; s < numSpouse[p]; ++s)
						if (ind(p).info(spouseIdx[s]) == m) {
							valid = false;
							break;
						}
				}
				if (valid) {
					ind(p).setInfo(m, spouseIdx[numSpouse[p]]);
					++numSpouse[p];
				}
			}
			if (numSpouse[m] < maxSpouse) {
				bool valid = true;
				// if sex is not interested
				if ((ind(p).sex() == Male && relSex == FemaleOnly) ||
				    (ind(p).sex() == Female && relSex == MaleOnly))
					valid = false;
				// duplicate spouse
				if (valid) {
					for (size_t s = 0; s < numSpouse[m]; ++s)
						if (ind(m).info(spouseIdx[s]) == p) {
							valid = false;
							break;
						}
				}
				if (valid) {
					ind(m).setInfo(p, spouseIdx[numSpouse[m]]);
					++numSpouse[m];
				}
			}                                                                                               // idx
		}                                                                                                   // ancestal generations
		// set the rest of the field to -1
		for (size_t idx = 0; idx < popSize(); ++idx) {
			for (size_t ns = numSpouse[idx]; ns < maxSpouse; ++ns)
				ind(idx).setInfo(-1, spouseIdx[ns]);
		}
	}     // each genearation
}


void pedigree::locateSpouseByID(SexChoice relSex, const vectorstr & relFields, UINT topGen)
{
	// if we are using individual ID, the algorithm is slower, but it can handle
	// overlapping generations.
	//
	// What we essentially do is getting all the couples and assign them...
	UINT maxSpouse = relFields.size();

	vectori spouseIdx(maxSpouse);

	// clear all fields
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (size_t i = 0; i < maxSpouse; ++i) {
			spouseIdx[i] = infoIdx(relFields[i]);
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
			double f = ind(i).info(m_fatherIdx);
			double m = ind(i).info(m_motherIdx);
			if (f >= 0 && m >= 0)
				couples.push_back(couple(static_cast<ULONG>(f), static_cast<ULONG>(m)));
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
			individual & fa = indByID(p);
			individual & ma = indByID(m);
			if (numSpouse[p] < maxSpouse) {
				bool valid = true;
				// if sex is not interested
				if ((ma.sex() == Male && relSex == FemaleOnly) ||
				    (ma.sex() == Female && relSex == MaleOnly))
					valid = false;
				// duplicate spouse
				if (valid) {
					for (size_t s = 0; s < numSpouse[p]; ++s)
						if (fa.info(spouseIdx[s]) == m) {
							valid = false;
							break;
						}
				}
				if (valid) {
					fa.setInfo(m, spouseIdx[numSpouse[p]]);
					++numSpouse[p];
				}
			}
			if (numSpouse[m] < maxSpouse) {
				bool valid = true;
				// if sex is not interested
				if ((fa.sex() == Male && relSex == FemaleOnly) ||
				    (fa.sex() == Female && relSex == MaleOnly))
					valid = false;
				// duplicate spouse
				if (valid) {
					for (size_t s = 0; s < numSpouse[m]; ++s)
						if (ind(m).info(spouseIdx[s]) == p) {
							valid = false;
							break;
						}
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


void pedigree::locateOffspringByIdx(SexChoice relSex, const vectorstr & relFields, UINT topGen)
{
	UINT maxOffspring = relFields.size();

	vectori offspringIdx(maxOffspring);

	useAncestralGen(0);
	for (size_t i = 0; i < maxOffspring; ++i) {
		offspringIdx[i] = infoIdx(relFields[i]);
		// clear these fields for the last generation
		for (IndInfoIterator ptr = infoBegin(offspringIdx[i]);
		     ptr != infoEnd(offspringIdx[i]); ++ptr)
			*ptr = static_cast<InfoType>(-1);
	}

	DBG_WARNING(topGen == 0, "Offspring can not be located because there is no parental generation.");
	// start from the parental generation
	for (unsigned ans = 1; ans <= topGen; ++ans) {
		vectoru numOffspring;
		// for each type of parental relationship
		vectoru parentFields;
		if (m_fatherIdx != -1)
			parentFields.push_back(m_fatherIdx);
		if (m_motherIdx != -1)
			parentFields.push_back(m_motherIdx);
		for (vectoru::const_iterator field = parentFields.begin();
		     field != parentFields.end(); ++field) {
			useAncestralGen(ans - 1);
			vectorf parent = indInfo(*field);
			DBG_DO(DBG_POPULATION, cerr << "Parents " << parent << endl);
			//
			useAncestralGen(ans);
			if (numOffspring.empty())
				numOffspring.resize(popSize(), 0);
			std::fill(numOffspring.begin(), numOffspring.end(), 0);
			//
			for (size_t idx = 0; idx < parent.size(); ++idx) {
				DBG_ASSERT(fcmp_ne(parent[idx], -1), ValueError, "Invalid parental index (-1)");
				ULONG p = static_cast<ULONG>(parent[idx]);
				DBG_ASSERT(p < popSize(), IndexError, "Parental index out of range of 0 ~ " + toStr(popSize() - 1));
				if (numOffspring[p] < maxOffspring) {
					Sex offSex = relSex == AnySex ? Male : ancestor(idx, ans - 1).sex();
					if ((relSex == MaleOnly && offSex != Male) ||
					    (relSex == FemaleOnly && offSex != Female) ||
					    (relSex == SameSex && offSex != ind(p).sex()) ||
					    (relSex == OppositeSex && offSex == ind(p).sex()))
						continue;
					ind(p).setInfo(idx, offspringIdx[numOffspring[p]]);
					++numOffspring[p];
				}
			}                                                                                               // idx
		}                                                                                                   // ancestal generations
		// only the last gen is cleared in advance.
		// other generations should be done...
		for (size_t idx = 0; idx < popSize(); ++idx) {
			for (size_t no = numOffspring[idx]; no < maxOffspring; ++no)
				ind(idx).setInfo(-1, offspringIdx[no]);
		}
	}
	useAncestralGen(0);

}


void pedigree::locateOffspringByID(SexChoice relSex, const vectorstr & relFields, UINT topGen)
{
	UINT maxOffspring = relFields.size();

	vectori offspringIdx(maxOffspring);

	// clear offspring fields
	for (unsigned ans = topGen; ans >= 0; --ans) {
		useAncestralGen(ans);
		for (size_t i = 0; i < maxOffspring; ++i) {
			offspringIdx[i] = infoIdx(relFields[i]);
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
			double p = ind(i).info(m_fatherIdx);
			double m = ind(i).info(m_motherIdx);
			if (p < 0 || m < 0)
				continue;
			try {
				// but the parents might not exist
				ULONG pp = static_cast<ULONG>(p);
				ULONG mm = static_cast<ULONG>(m);
				individual & fa = indByID(pp);
				individual & ma = indByID(mm);
				// add child as father's offspring
				if (numOffspring[pp] < maxOffspring) {
					Sex offSex = relSex == AnySex ? Male : ind(i).sex();
					if ((relSex == MaleOnly && offSex != Male) ||
					    (relSex == FemaleOnly && offSex != Female) ||
					    (relSex == SameSex && offSex != fa.sex()) ||
					    (relSex == OppositeSex && offSex == fa.sex()))
						continue;
					fa.setInfo(ind(i).info(m_idIdx), offspringIdx[numOffspring[pp]]);
					++numOffspring[pp];
				}
				// add child as mother's offspring
				if (numOffspring[mm] < maxOffspring) {
					Sex offSex = relSex == AnySex ? Male : ind(i).sex();
					if ((relSex == MaleOnly && offSex != Male) ||
					    (relSex == FemaleOnly && offSex != Female) ||
					    (relSex == SameSex && offSex != fa.sex()) ||
					    (relSex == OppositeSex && offSex == fa.sex()))
						continue;
					ma.setInfo(ind(i).info(m_idIdx), offspringIdx[numOffspring[mm]]);
					++numOffspring[mm];
				}
			} catch (...) {
				// does not care if a parent cannot be found
			}
		}
	}

}


void pedigree::locateSiblingByIdx(SexChoice relSex, const vectorstr & relFields, UINT topGen)
{
	DBG_WARNING(topGen == 0, "Sibling can not be located because there is no parental generation.");

	UINT maxSibling = relFields.size();
	vectori siblingIdx(maxSibling);
	for (size_t i = 0; i < maxSibling; ++i)
		siblingIdx[i] = infoIdx(relFields[i]);

	// start from the parental generation
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		// if top generation, no information about sibling
		if (ans == ancestralGens()) {
			for (IndIterator it = indIterator(); it.valid(); ++it)
				for (size_t i = 0; i < maxSibling; ++i)
					it->setInfo(-1, siblingIdx[i]);
			continue;
		}
		// when parents information are available.
		vectoru numSibling(popSize(), 0);
		// for each type of parental relationship
		map<pair<ULONG, ULONG>, vector<ULONG> > par_map;

		vectoru parentFields;
		if (m_fatherIdx != -1)
			parentFields.push_back(m_fatherIdx);
		if (m_motherIdx != -1)
			parentFields.push_back(m_motherIdx);
		for (vectoru::const_iterator field = parentFields.begin();
		     field != parentFields.end(); ++field) {
			vectorf parent = indInfo(*field);
			for (size_t idx = 0; idx != parent.size(); ++idx) {
				pair<ULONG, ULONG> parents(static_cast<ULONG>(parent[idx]), 0);
				map<pair<ULONG, ULONG>, vector<ULONG> >::iterator item = par_map.find(parents);
				if (item == par_map.end())
					par_map[parents] = vector<ULONG>(1, idx);
				else
					item->second.push_back(idx);
			}
		}
		// now, each par_map has offsprings as siblings
		map<pair<ULONG, ULONG>, vector<ULONG> >::const_iterator it = par_map.begin();
		map<pair<ULONG, ULONG>, vector<ULONG> >::const_iterator end = par_map.end();
		for (; it != end; ++it) {
			const vector<ULONG> & sibs = it->second;
			if (sibs.size() <= 1)
				continue;
			//
			vector<Sex> sexes(sibs.size());
			for (size_t i = 0; i < sibs.size(); ++i)
				sexes[i] = ind(sibs[i]).sex();
			//
			for (size_t i = 0; i < sibs.size(); ++i) {
				for (size_t j = 0; j < sibs.size(); ++j) {
					if (i == j)
						continue;
					if ((relSex == MaleOnly && sexes[j] == Female) ||
					    (relSex == FemaleOnly && sexes[j] == Male) ||
					    (relSex == SameSex && sexes[j] != sexes[j]) ||
					    (relSex == OppositeSex && sexes[i] == sexes[j]))
						continue;
					if (numSibling[sibs[i]] < maxSibling) {
						ind(sibs[i]).setInfo(sibs[j], siblingIdx[numSibling[sibs[i]]]);
						++numSibling[sibs[i]];
					}
				}
			}                                                                                               // idx
		}
		// set the rest of the field to -1
		for (size_t idx = 0; idx < popSize(); ++idx) {
			for (size_t no = numSibling[idx]; no < maxSibling; ++no)
				ind(idx).setInfo(-1, siblingIdx[no]);
		}
	}

}


void pedigree::locateSiblingByID(SexChoice relSex, const vectorstr & relFields, UINT topGen)
{
	UINT maxSibling = relFields.size();
	vectori siblingIdx(maxSibling);

	for (size_t i = 0; i < maxSibling; ++i)
		siblingIdx[i] = infoIdx(relFields[i]);

	// clear all fields
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
			double f = ind(i).info(m_fatherIdx);
			double m = ind(i).info(m_motherIdx);
			if (f >= 0)
				families[static_cast<ULONG>(f)].push_back(static_cast<ULONG>(ind(i).info(m_idIdx)));
			if (m >= 0)
				families[static_cast<ULONG>(m)].push_back(static_cast<ULONG>(ind(i).info(m_idIdx)));
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
				individual & child = indByID(sibs[i]);
				for (size_t j = 0; j < sibs.size(); ++j) {
					if (i == j)
						continue;
					try {
						individual & sibling = indByID(sibs[j]);

						if (numSibling[sibs[i]] < maxSibling) {
							bool valid = true;
							// if sex is not interested
							if ((sibling.sex() == Male && relSex == FemaleOnly) ||
							    (sibling.sex() == Female && relSex == MaleOnly))
								valid = false;
							// duplicate spouse
							if (valid) {
								for (size_t s = 0; s < numSibling[sibs[i]]; ++s)
									if (child.info(siblingIdx[s]) == sibs[j]) {
										valid = false;
										break;
									}
							}
							if (valid) {
								child.setInfo(sibs[j], siblingIdx[numSibling[sibs[i]]]);
								++numSibling[sibs[i]];
							}
						}
						if (numSibling[sibs[j]] < maxSibling) {
							bool valid = true;
							if ((sibling.sex() == Male && relSex == FemaleOnly) ||
							    (sibling.sex() == Female && relSex == MaleOnly))
								// if sex is not interested
								valid = false;
							// duplicate spouse
							if (valid) {
								for (size_t s = 0; s < numSibling[sibs[j]]; ++s)
									if (sibling.info(siblingIdx[s]) == sibs[i]) {
										valid = false;
										break;
									}
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


void pedigree::locateFullSiblingByIdx(SexChoice relSex, const vectorstr & relFields, UINT topGen)
{
	UINT maxSibling = relFields.size();
	vectori siblingIdx(maxSibling);

	for (size_t i = 0; i < maxSibling; ++i)
		siblingIdx[i] = infoIdx(relFields[i]);

	DBG_WARNING(topGen == 0, "Sibling can not be located because there is no parental generation.");
	// start from the parental generation
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		// if top generation, no information about sibling
		if (ans == ancestralGens()) {
			for (IndIterator it = indIterator(); it.valid(); ++it)
				for (size_t i = 0; i < maxSibling; ++i)
					it->setInfo(-1, siblingIdx[i]);
			continue;
		}
		// when parents information are available.
		vectoru numSibling(popSize(), 0);
		// for each type of parental relationship
		map<pair<ULONG, ULONG>, vector<ULONG> > par_map;

		vectorf father = indInfo(m_fatherIdx);
		vectorf mother = indInfo(m_motherIdx);
		for (size_t idx = 0; idx != father.size(); ++idx) {
			pair<ULONG, ULONG> parents(static_cast<ULONG>(father[idx]),
			                           static_cast<ULONG>(mother[idx]));
			map<pair<ULONG, ULONG>, vector<ULONG> >::iterator item = par_map.find(parents);
			if (item == par_map.end())
				par_map[parents] = vector<ULONG>(1, idx);
			else
				item->second.push_back(idx);
		}
		// now, each par_map has offsprings as siblings
		map<pair<ULONG, ULONG>, vector<ULONG> >::const_iterator it = par_map.begin();
		map<pair<ULONG, ULONG>, vector<ULONG> >::const_iterator end = par_map.end();
		for (; it != end; ++it) {
			const vector<ULONG> & sibs = it->second;
			if (sibs.size() <= 1)
				continue;
			//
			vector<Sex> sexes(sibs.size());
			for (size_t i = 0; i < sibs.size(); ++i)
				sexes[i] = ind(sibs[i]).sex();
			//
			for (size_t i = 0; i < sibs.size(); ++i) {
				for (size_t j = 0; j < sibs.size(); ++j) {
					if (i == j)
						continue;
					if ((relSex == MaleOnly && sexes[j] == Female) ||
					    (relSex == FemaleOnly && sexes[j] == Male) ||
					    (relSex == SameSex && sexes[j] != sexes[j]) ||
					    (relSex == OppositeSex && sexes[i] == sexes[j]))
						continue;
					if (numSibling[sibs[i]] < maxSibling) {
						ind(sibs[i]).setInfo(sibs[j], siblingIdx[numSibling[sibs[i]]]);
						++numSibling[sibs[i]];
					}
				}
			}                                                                                               // idx
		}
		// set the rest of the field to -1
		for (size_t idx = 0; idx < popSize(); ++idx) {
			for (size_t no = numSibling[idx]; no < maxSibling; ++no)
				ind(idx).setInfo(-1, siblingIdx[no]);
		}
	}
}


void pedigree::locateFullSiblingByID(SexChoice relSex, const vectorstr & relFields, UINT topGen)
{
	UINT maxSibling = relFields.size();
	vectori siblingIdx(maxSibling);

	for (size_t i = 0; i < maxSibling; ++i)
		siblingIdx[i] = infoIdx(relFields[i]);

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
			double f = ind(i).info(m_fatherIdx);
			double m = ind(i).info(m_motherIdx);
			if (f >= 0 && m >= 0)
				families[couple(static_cast<ULONG>(f), static_cast<ULONG>(m))].push_back(static_cast<ULONG>(ind(i).info(m_idIdx)));
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
				individual & child = indByID(sibs[i]);
				for (size_t j = 0; j < sibs.size(); ++j) {
					if (i == j)
						continue;
					individual & sibling = indByID(sibs[j]);

					if (numSibling[sibs[i]] < maxSibling) {
						bool valid = true;
						// if sex is not interested
						if ((sibling.sex() == Male && relSex == FemaleOnly) ||
						    (sibling.sex() == Female && relSex == MaleOnly))
							valid = false;
						// duplicate spouse
						if (valid) {
							for (size_t s = 0; s < numSibling[sibs[i]]; ++s)
								if (child.info(siblingIdx[s]) == sibs[j]) {
									valid = false;
									break;
								}
						}
						if (valid) {
							child.setInfo(sibs[j], siblingIdx[numSibling[sibs[i]]]);
							++numSibling[sibs[i]];
						}
					}
					if (numSibling[sibs[j]] < maxSibling) {
						bool valid = true;
						if ((sibling.sex() == Male && relSex == FemaleOnly) ||
						    (sibling.sex() == Female && relSex == MaleOnly))
							// if sex is not interested
							valid = false;
						// duplicate spouse
						if (valid) {
							for (size_t s = 0; s < numSibling[sibs[j]]; ++s)
								if (sibling.info(siblingIdx[s]) == sibs[i]) {
									valid = false;
									break;
								}
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


void pedigree::locateSpouseAndOffspringByIdx(SexChoice relSex, const vectorstr & relFields, UINT topGen)
{
	UINT maxOffspring = relFields.size() - 1;

	int spouseIdx = infoIdx(relFields[0]);
	// clear these fields for the last generation
	IndInfoIterator ptr = infoBegin(spouseIdx);
	IndInfoIterator ptrEnd = infoEnd(spouseIdx);

	for (; ptr != ptrEnd; ++ptr)
		*ptr = static_cast<InfoType>(-1);
	vectori offspringIdx(maxOffspring);
	for (size_t i = 0; i < maxOffspring; ++i) {
		offspringIdx[i] = infoIdx(relFields[i + 1]);
		// clear these fields for the last generation
		ptr = infoBegin(offspringIdx[i]);
		ptrEnd = infoEnd(offspringIdx[i]);
		for (; ptr != ptrEnd; ++ptr)
			*ptr = static_cast<InfoType>(-1);
	}

	// start from the parental generation
	for (unsigned ans = 1; ans <= topGen; ++ans) {
		vectoru numOffspring;
		// go to offspring generation
		useAncestralGen(ans - 1);
		vectorf father = indInfo(m_fatherIdx);
		vectorf mother = indInfo(m_motherIdx);
		//
		useAncestralGen(ans);
		if (numOffspring.empty())
			numOffspring.resize(popSize(), 0);
		std::fill(numOffspring.begin(), numOffspring.end(), 0);
		setIndInfo(-1, spouseIdx);
		//
		for (size_t idx = 0; idx < father.size(); ++idx) {
			DBG_FAILIF(fcmp_eq(father[idx], -1) || fcmp_eq(mother[idx], -1), ValueError,
				"Invalid parental index (-1)");
			ULONG p = static_cast<ULONG>(father[idx]);
			ULONG m = static_cast<ULONG>(mother[idx]);
			DBG_ASSERT(p < popSize() && m < popSize(), IndexError,
				"Parental index out of range of 0 ~ " + toStr(popSize() - 1));
			// if new parents
			if (numOffspring[p] == 0 && numOffspring[m] == 0) {
				if (relSex != AnySex) {
					Sex offSex = ancestor(idx, ans - 1).sex();
					if ((offSex == Male && relSex == FemaleOnly) ||
					    (offSex == Female && relSex == MaleOnly))
						continue;
				}
				ind(p).setInfo(m, spouseIdx);
				ind(m).setInfo(p, spouseIdx);
				ind(p).setInfo(idx, offspringIdx[0]);
				ind(m).setInfo(idx, offspringIdx[0]);
				++numOffspring[p];
				++numOffspring[m];
			} else if (numOffspring[p] == numOffspring[m]) {
				// not the original spouse
				if (ind(p).intInfo(spouseIdx) != static_cast<int>(m))
					continue;
				// no room for another offspring
				if (numOffspring[p] >= maxOffspring)
					continue;
				// sex does not match
				if (relSex != AnySex) {
					Sex offSex = ancestor(idx, ans - 1).sex();
					if ((offSex == Male && relSex == FemaleOnly) ||
					    (offSex == Female && relSex == MaleOnly))
						continue;
				}
				// great! The same parents
				ind(p).setInfo(idx, offspringIdx[numOffspring[p]]);
				ind(m).setInfo(idx, offspringIdx[numOffspring[m]]);
				++numOffspring[p];
				++numOffspring[m];
			} else
				continue;
		}                                                                                                   // ancestal generations
		// set the rest of the field to -1
		for (size_t idx = 0; idx < popSize(); ++idx) {
			for (size_t ns = numOffspring[idx]; ns < maxOffspring; ++ns)
				ind(idx).setInfo(-1, offspringIdx[ns]);
		}
	}

}


void pedigree::locateSpouseAndOffspringByID(SexChoice relSex, const vectorstr & relFields, UINT topGen)
{
	UINT maxOffspring = relFields.size() - 1;

	int spouseIdx = infoIdx(relFields[0]);
	// clear these fields for the last generation
	IndInfoIterator ptr = infoBegin(spouseIdx);
	IndInfoIterator ptrEnd = infoEnd(spouseIdx);

	vectori offspringIdx(maxOffspring);

	// clear all fields
	for (unsigned ans = 0; ans <= topGen; ++ans) {
		useAncestralGen(ans);
		for (; ptr != ptrEnd; ++ptr)
			*ptr = static_cast<InfoType>(-1);
		for (size_t i = 0; i < maxOffspring; ++i) {
			offspringIdx[i] = infoIdx(relFields[i + 1]);
			// clear these fields for the last generation
			ptr = infoBegin(offspringIdx[i]);
			ptrEnd = infoEnd(offspringIdx[i]);
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
			double f = ind(i).info(m_fatherIdx);
			double m = ind(i).info(m_motherIdx);
			if (f >= 0 && m >= 0)
				families[couple(static_cast<ULONG>(f), static_cast<ULONG>(m))].push_back(static_cast<ULONG>(ind(i).info(m_idIdx)));
		}
	}
	// look in each family
	std::map<ULONG, UINT> numSibling;
	map<couple, vectoru>::iterator it = families.begin();
	map<couple, vectoru>::iterator itEnd = families.end();
	std::map<ULONG, UINT> numSpouse;
	std::map<ULONG, UINT> numOffspring;
	for (; it != itEnd; ++it) {
		ULONG p = it->first.first;
		ULONG m = it->first.second;
		// these guys share two parents
		const vectoru & offspring = it->second;
		try {
			// someone has kids with someone else
			if (numSpouse[p] == 1 || numSpouse[m] == 1)
				continue;
			individual & fa = indByID(p);
			individual & ma = indByID(m);
			// spouse
			fa.setInfo(m, spouseIdx);
			ma.setInfo(p, spouseIdx);
			++numSpouse[p];
			++numSpouse[m];
			// offspring
			for (size_t j = 0; j < offspring.size(); ++j) {
				if (numOffspring[p] >= maxOffspring || numOffspring[m] >= maxOffspring)
					continue;
				try {
					individual & child = indByID(offspring[j]);
					bool valid = true;
					// if sex is not interested
					if ((child.sex() == Male && relSex == FemaleOnly) ||
					    (child.sex() == Female && relSex == MaleOnly))
						valid = false;
					// duplicate child
					if (valid) {
						for (size_t s = 0; s < numOffspring[p]; ++s)
							if (fa.info(offspringIdx[s]) == offspring[j]) {
								valid = false;
								break;
							}
					}
					if (valid) {
						fa.setInfo(offspring[j], offspringIdx[numOffspring[p]]);
						ma.setInfo(offspring[j], offspringIdx[numOffspring[m]]);
						++numOffspring[p];
						++numOffspring[m];
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


bool pedigree::traceRelatives(const vectoru & pathGen,
                              const stringMatrix & pathFieldsMatrix, const vectori & pathSex,
                              const vectorstr & resultFields)
{
	if (pathGen.empty())
		return true;

	const matrixstr & pathFields = pathFieldsMatrix.elems();

	DBG_ASSERT(m_idIdx == -1 && pathGen.size() == pathFields.size() + 1, ValueError,
		"Parameter pathGen should be one element longer than pathFields");
	DBG_FAILIF(!pathSex.empty() && pathSex.size() != pathFields.size(),
		ValueError,
		"Parameter pathSex, if given, should have the same length of pathFields");

	vectori resultIdx(resultFields.size());
	for (size_t i = 0; i < resultIdx.size(); ++i)
		resultIdx[i] = infoIdx(resultFields[i]);
	// start generation, and at which result will be saved.
	useAncestralGen(pathGen[0]);
	vectoru numResult(popSize(), 0);
	UINT maxResult = resultIdx.size();
	// clear values
	for (IndIterator ind = indIterator(); ind.valid(); ++ind)
		for (size_t i = 0; i < maxResult; ++i)
			ind->setInfo(-1, resultIdx[i]);
	// convert pathFields to pathIdx
	intMatrix pathIdx(pathFields.size());
	for (size_t i = 0; i < pathFields.size(); ++i) {
		pathIdx[i] = vectori(pathFields[i].size());
		for (size_t j = 0; j < pathFields[i].size(); ++j)
			pathIdx[i][j] = infoIdx(pathFields[i][j]);
	}
	// convert pathSex to type SexChoices
	vector<SexChoice> sexes(pathIdx.size(), AnySex);
	for (size_t i = 0; i < pathSex.size(); ++i)
		sexes[i] = static_cast<SexChoice>(pathSex[i]);

	ULONG idx = 0;
	for (IndIterator ind = indIterator(); ind.valid(); ++ind, ++idx) {
		// start from one individual from pathGen[0]
		Sex mySex = ind->sex();
		vectoru inds = vectoru(1, m_idIdx == -1 ? idx : static_cast<ULONG>(ind->info(m_idIdx)));
		// go through the path
		for (size_t path = 0; path < pathFields.size(); ++path) {
			DBG_DO(DBG_POPULATION, cerr << "Start of path " << path
				                        << " : " << inds << endl);
			UINT fromGen = pathGen[path];
			UINT toGen = pathGen[path + 1];

			if (m_idIdx == -1 && (fromGen > ancestralGens() || toGen > ancestralGens())) {
				DBG_WARNING(true, "Insufficient ancestral generations to trace relatives.");
				return false;
			}

			const vectori & fields = pathIdx[path];
			SexChoice sex = sexes[path];

			vectoru newInds;
			// for all individuals
			for (size_t i = 0; i < inds.size(); ++i) {
				// for all fields
				for (size_t s = 0; s < fields.size(); ++s) {
					InfoType sIdx = m_idIdx == -1 ? ancestor(inds[i], fromGen).info(fields[s])
									: indByID(inds[i], fromGen).info(fields[s]);
					if (sIdx < 0)
						continue;
					Sex indSex = m_idIdx == -1 ? ancestor(static_cast<ULONG>(sIdx), toGen).sex() :
					             indByID(static_cast<ULONG>(sIdx), toGen).sex();
					if ((sex == MaleOnly && indSex == Female) ||
					    (sex == FemaleOnly && indSex == Male) ||
					    (sex == SameSex && indSex != mySex) ||
					    (sex == OppositeSex && indSex == mySex))
						continue;
					newInds.push_back(static_cast<ULONG>(sIdx));
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
	return true;
}


// const unsigned long UnusedIndividual = std::numeric_limits<unsigned long>::max();
//
// pedigree::pedigree(int numParents, const string & pedfile)
//  : m_numParents(numParents)
// {
//  DBG_ASSERT(numParents == 1 || numParents == 2, ValueError,
//      "Individuals in a pedigree can have one or two parents");
//
//  if (!pedfile.empty())
//      load(pedfile);
// }
////
// vectoru pedigree::subPopSizes(ULONG gen)
// {
//  CHECK_GEN(gen);
//  return m_pedSize[gen];
// }
//
//
// ULONG pedigree::subPopSize(ULONG gen, SubPopID subPop)
// {
//  CHECK_GEN(gen);
//  CHECK_SUBPOP(gen, subPop);
//
//  return m_pedSize[gen][subPop];
// }
//
//
// void pedigree::addGen(const vectoru & sizes)
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
// void pedigree::load(const string & filename)
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
//          "Line: " + toStr(m_paternal.size() + 1) + ", Individuals read: "
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
// void pedigree::loadInfo(const string & filename, const string & name)
// {
//  vectorstr names(1, name);
//
//  loadInfo(filename, names);
// }
//
//
// void pedigree::loadInfo(const string & filename, const vectorstr & names)
// {
//  ifstream afs(filename.c_str());
//
//  DBG_FAILIF(!afs, SystemError, "Can not open auxiliary information pedigree" + filename + " to read");
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
//          "Information pedigree is larger than parental pedigree");
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
//          "Error loading information pedigree");
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
// void pedigree::addInfo(const string & name, double init)
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
//          "Error loading information pedigree");
//      for (size_t i = 0; i < size; ++i)
//          m_info[gen][i].push_back(init);
//  }
// }
//
//
// void pedigree::save(const string & filename)
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
// void pedigree::saveInfo(const string & filename, const string & name)
// {
//  vectorstr names(1, name);
//
//  saveInfo(filename, names);
// }
//
//
// void pedigree::saveInfo(const string & filename, const vectorstr & names)
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
// void pedigree::selectIndividuals(const vectoru & inds)
// {
//  DBG_FAILIF(m_paternal.empty(), ValueError,
//      "Can not select individuals from an empty pedigree");
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
// void pedigree::markUnrelated()
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
// void pedigree::removeUnrelated(bool shift_index)
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
//                      ValueError, "Inconsistent maternal and matermal pedigree");
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


