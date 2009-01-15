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

pedigree::pedigree(const population & pop, const vectoru & loci,
	const vectorstr & infoFields, int ancGen,
	const string & fatherField, const string & motherField)
	: m_fatherField(fatherField), m_motherField(motherField),
	m_fatherIdx(-1), m_motherIdx(-1)
{
	vectorstr extractFields = infoFields;

	if (!m_fatherField.empty() && find(extractFields.begin(), extractFields.end(), fatherField) == extractFields.end())
		extractFields.push_back(fatherField);
	if (!m_motherField.empty() && find(extractFields.begin(), extractFields.end(), motherField) == extractFields.end())
		extractFields.push_back(motherField);
	//
	population & ped = pop.extract(false, string(), loci.size() != pop.totNumLoci(),
		loci, extractFields.size() != pop.infoSize(), extractFields, ancGen, NULL);

	swap(ped);

	DBG_FAILIF(!m_fatherField.empty() && !pop.hasInfoField(m_fatherField), ValueError,
		"Invalid father information field " + m_fatherField);
	DBG_FAILIF(!m_motherField.empty() && !pop.hasInfoField(m_motherField), ValueError,
		"Invalid mother information field " + m_motherField);

	m_fatherIdx = m_fatherField.empty() ? -1 : infoIdx(fatherField);
	m_motherIdx = m_motherField.empty() ? -1 : infoIdx(motherField);
}


pedigree::pedigree(const pedigree & rhs) :
	population(rhs),
	m_fatherField(rhs.m_fatherField), m_motherField(rhs.m_motherField),
	m_fatherIdx(rhs.m_fatherIdx), m_motherIdx(rhs.m_motherIdx)
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


void pedigree::locateRelatives(uintList fullRelType, const vectorstr & relFields, int ancGen)
{
	const vectorlu & fullType = fullRelType.elems();

	if (fullType.empty() || relFields.empty())
		return;

	RelativeType relType = static_cast<RelativeType>(fullType[0]);

	DBG_FAILIF(fullType.size() > 2, ValueError, "Unrecognized relative type.");

	SexChoice relSex = fullType.size() == 2 ? static_cast<SexChoice>(fullType[1]) : AnySex;

	DBG_ASSERT(relSex == AnySex || relSex == MaleOnly || relSex == FemaleOnly
		|| relSex == SameSex || relSex == OppositeSex, ValueError,
		"Relative sex can only be one of AnySex, MaleOnly, FemaleOnly, SameSex or OppositeSex.");

	UINT topGen = ancGen == -1 ? ancestralGens() : std::min(ancestralGens(), static_cast<UINT>(ancGen));

	if (relType == Self) {
		DBG_ASSERT(relFields.size() == 1, ValueError,
			"Please provide one information field to store Self individuals");
		UINT fieldIdx = infoIdx(relFields[0]);

		for (size_t ans = 0; ans <= topGen; ++ans) {
			useAncestralGen(ans);
			for (size_t idx = 0; idx < popSize(); ++idx)
				ind(idx).setInfo(idx, fieldIdx);
		}
		useAncestralGen(0);
	} else if (relType == Spouse) {
		DBG_ASSERT(numParents() == 2, ValueError,
			"This relative only exists when there are two parents for each indidivual");

		DBG_ASSERT(relFields.size() >= 1, ValueError,
			"Please provide at least one information field to store Self individuals");

		DBG_FAILIF(relSex == SameSex, ValueError, "Can not locate spouses with the same sex");

		UINT maxSpouse = relFields.size();

		vectori spouseIdx(maxSpouse);
		for (size_t i = 0; i < maxSpouse; ++i) {
			spouseIdx[i] = infoIdx(relFields[i]);
			// clear these fields for the last generation
			for (IndInfoIterator ptr = infoBegin(spouseIdx[i]);
			     ptr != infoEnd(spouseIdx[i]); ++ptr)
				*ptr = static_cast<InfoType>(-1);
		}

		DBG_WARNING(topGen == 0, "Spouse can not be located because there is no parental generation.");
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
				}                                                                                           // idx
			}                                                                                               // ancestal generations
			// set the rest of the field to -1
			for (size_t idx = 0; idx < popSize(); ++idx) {
				for (size_t ns = numSpouse[idx]; ns < maxSpouse; ++ns)
					ind(idx).setInfo(-1, spouseIdx[ns]);
			}
		}
		useAncestralGen(0);
	} else if (relType == Offspring) {
		DBG_ASSERT(relFields.size() >= 1, ValueError,
			"Please provide at least one information field to store offspring");

		UINT maxOffspring = relFields.size();

		vectori offspringIdx(maxOffspring);
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
				DBG_DO(DBG_POPULATION, cout << "Parents " << parent << endl);
				//
				useAncestralGen(ans);
				if (numOffspring.empty())
					numOffspring.resize(popSize(), 0);
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
				}                                                                                           // idx
			}                                                                                               // ancestal generations
			// only the last gen is cleared in advance.
			// other generations should be done...
			for (size_t idx = 0; idx < popSize(); ++idx) {
				for (size_t no = numOffspring[idx]; no < maxOffspring; ++no)
					ind(idx).setInfo(-1, offspringIdx[no]);
			}
		}
		useAncestralGen(0);
	} else if (relType == Sibling || relType == FullSibling) {
		DBG_ASSERT(relFields.size() >= 1, ValueError,
			"Please provide at least one information field to store offspring");

		DBG_FAILIF(relType == FullSibling && numParents() != 2, ValueError,
			"Please provide two parental information fields");
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
				for (IndIterator it = indBegin(); it.valid(); ++it)
					for (size_t i = 0; i < maxSibling; ++i)
						it->setInfo(-1, siblingIdx[i]);
				continue;
			}
			// when parents information are available.
			vectoru numSibling(popSize(), 0);
			// for each type of parental relationship
			map<pair<ULONG, ULONG>, vector<ULONG> > par_map;

			if (relType == Sibling) { // one or two parents
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
			} else {
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
				}                                                                                           // idx
			}
			// set the rest of the field to -1
			for (size_t idx = 0; idx < popSize(); ++idx) {
				for (size_t no = numSibling[idx]; no < maxSibling; ++no)
					ind(idx).setInfo(-1, siblingIdx[no]);
			}
		}
		useAncestralGen(0);
	} else {
		throw ValueError("Unrecognized relative type");
	}
}


bool pedigree::traceRelatives(const vectoru & pathGen,
                              const stringMatrix & pathFields, const vectori & pathSex,
                              const vectorstr & resultFields)
{
	if (pathGen.empty())
		return true;

	DBG_ASSERT(pathGen.size() == pathFields.size() + 1, ValueError,
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
	for (IndIterator ind = indBegin(); ind.valid(); ++ind)
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
	for (IndIterator ind = indBegin(); ind.valid(); ++ind, ++idx) {
		// start from one individual from pathGen[0]
		Sex mySex = ind->sex();
		vectorlu inds = vectorlu(1, idx);
		// go through the path
		for (size_t path = 0; path < pathFields.size(); ++path) {
			DBG_DO(DBG_POPULATION, cout << "Start of path " << path
				                        << " : " << inds << endl);
			UINT fromGen = pathGen[path];
			UINT toGen = pathGen[path + 1];

			if (fromGen > ancestralGens() || toGen > ancestralGens()) {
				DBG_WARNING(true, "Insufficient ancestral generations to trace relatives.");
				return false;
			}

			const vectori & fields = pathIdx[path];
			SexChoice sex = sexes[path];

			vectorlu newInds;
			// for all individuals
			for (size_t i = 0; i < inds.size(); ++i) {
				// for all fields
				for (size_t s = 0; s < fields.size(); ++s) {
					InfoType sIdx = ancestor(inds[i], fromGen).info(fields[s]);
					if (sIdx < 0)
						continue;
					Sex indSex = ancestor(static_cast<ULONG>(sIdx), toGen).sex();
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
		DBG_DO(DBG_POPULATION, cout << "Ind " << idx << " has relatives " << inds << endl);
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
// vectorlu pedigree::subPopSizes(ULONG gen)
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
// void pedigree::addGen(const vectorlu & sizes)
// {
//  ULONG popSize = accumulate(sizes.begin(), sizes.end(), 0UL);
//
//  m_paternal.push_back(vectorlu(popSize));
//  if (m_numParents == 2)
//      m_maternal.push_back(vectorlu(popSize));
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
//      vectorlu values;
//      vectorlu sizes;
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
//      m_pedSize.push_back(vectorlu());
//      m_pedSize.back().swap(sizes);
//      if (m_numParents == 1) {
//          m_paternal.push_back(vectorlu());
//          m_paternal.back().swap(values);
//      } else if (m_numParents == 2) {
//          m_paternal.push_back(vectorlu(popSize));
//          m_maternal.push_back(vectorlu(popSize));
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
// void pedigree::selectIndividuals(const vectorlu & inds)
// {
//  DBG_FAILIF(m_paternal.empty(), ValueError,
//      "Can not select individuals from an empty pedigree");
//
//  size_t size = m_paternal.back().size();
//  vector<bool> used(size, false);
//  vectorlu::const_iterator it = inds.begin();
//  vectorlu::const_iterator it_end = inds.end();
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
//          vectorlu & curGen = m_paternal[gen];
//          vectorlu & nextGen = m_paternal[gen + 1];
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
//              vectorlu & curGen = m_maternal[gen];
//              vectorlu & nextGen = m_maternal[gen + 1];
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
//  vectorlu l_pat;
//  vectorlu l_mat;
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
