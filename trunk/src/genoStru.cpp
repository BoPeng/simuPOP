/**
 *  $File: genoStru.cpp $
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

#include "genoStru.h"

#include "boost_pch.hpp"

using namespace boost::lambda;

#include <set>
using std::set;

namespace simuPOP {

GenoStructure::GenoStructure(UINT ploidy, const vectoru & loci, const vectoru & chromTypes, bool haplodiploid,
	const vectorf & lociPos, const vectorstr & chromNames, const matrixstr & alleleNames,
	const vectorstr & lociNames, const vectorstr & infoFields)
	: m_ploidy(ploidy), m_numLoci(loci), m_chromTypes(),
	m_haplodiploid(haplodiploid), m_lociPos(lociPos), m_chromIndex(loci.size() + 1),
	m_chromNames(chromNames), m_alleleNames(alleleNames), m_lociNames(lociNames),
	m_infoFields(infoFields), m_lociPosMap(), m_refCount(0)
{
	DBG_ASSERT(ploidy >= 1, ValueError,
		(boost::format("Ploidy must be >= 1. Given %1%") % ploidy).str());

	// default: one chromosome, one locus
	// otherwise, Loci copies from loci
	if (loci.empty()) {
		m_numLoci.resize(1);
		m_numLoci[0] = 0;
		m_chromIndex.resize(2);
	}

	// chromosome type
	setChromTypes(chromTypes);

	// build chromosome index
	m_chromIndex[0] = 0;
	for (size_t i = 1; i <= m_numLoci.size(); ++i)
		m_chromIndex[i] = m_chromIndex[i - 1] + m_numLoci[i - 1];

	m_totNumLoci = m_chromIndex[m_numLoci.size()];

	DBG_ASSERT(m_lociNames.empty() || m_lociNames.size() == m_totNumLoci, ValueError,
		"Loci names, if specified, should be given to every loci");
	DBG_FAILIF(m_alleleNames.size() > 1 && m_alleleNames.size() != m_totNumLoci,
		ValueError, "Please specify a common set allele names or names for all loci.");

	// if lociPos not specified, use 1,2,3.. 1,2,3. .. on each chromosome.
	if (m_lociPos.empty()) {
		m_lociPos.resize(m_totNumLoci);
		for (size_t i = 0; i < m_numLoci.size(); ++i)
			for (size_t j = 0; j < m_numLoci[i]; j++)
				m_lociPos[m_chromIndex[i] + j] = static_cast<double>(j + 1);
	} else {                                                                            // check loci distance
		// loci distance, if specified, chould have length of chromosome.
		DBG_FAILIF(m_lociPos.size() != m_totNumLoci, ValueError,
			(boost::format("You should specify loci distance for every locus (%1%)") % m_totNumLoci).str());

		for (size_t ch = 0; ch < m_numLoci.size(); ++ch) {
			if (m_numLoci[ch] <= 1)
				continue;

			size_t begin = m_chromIndex[ch];
			size_t end = m_chromIndex[ch + 1];

			bool ordered = true;
			for (size_t j = begin + 1; j < end; ++j) {
				if (m_lociPos[j] <= m_lociPos[j - 1]) {
					ordered = false;
					break;
				}
			}
			if (ordered)
				continue;

			DBG_DO(DBG_POPULATION, cerr << "Loci on chromosome " << ch << " are unordered.");

			vectorf lociPos(m_lociPos.begin() + begin, m_lociPos.begin() + end);

			// rank
			vector<size_t> rank(m_numLoci[ch]);
			for (size_t i = 0; i < rank.size(); ++i)
				rank[i] = i;
			// sort according to value of pos
			std::sort(rank.begin(), rank.end(),
				var(lociPos)[_1] < var(lociPos)[_2]);
			// apply sorted loci positions
			for (size_t i = 0; i < rank.size(); ++i)
				m_lociPos[begin + i] = lociPos[rank[i]];
			// check again for loci duplication etc
			for (size_t j = begin + 1; j < end; ++j) {
				DBG_FAILIF(m_lociPos[j] <= m_lociPos[j - 1], ValueError,
					"Loci on the same chromosome should have different positions.");
			}
			// lociNames?
			if (!m_lociNames.empty()) {
				vectorstr lociNames(m_lociNames.begin() + begin, m_lociNames.begin() + end);
				for (size_t i = 0; i < rank.size(); ++i)
					m_lociNames[begin + i] = lociNames[rank[i]];
			}
			// allele names
			if (m_alleleNames.size() > 1) {
				matrixstr alleleNames(m_alleleNames.begin() + begin, m_alleleNames.begin() + end);
				for (size_t i = 0; i < rank.size(); ++i)
					m_alleleNames[begin + i] = alleleNames[rank[i]];
			}
		}
	}

	// set up a map for loci names and check uniqueness of the names
	if (!m_lociNames.empty()) {
		for (size_t i = 0; i < m_totNumLoci; ++i) {
			if (m_lociNames[i].empty())
				continue;
			if (m_lociNameMap.find(m_lociNames[i]) != m_lociNameMap.end())
				throw ValueError("Given loci names should be unique");
			m_lociNameMap[m_lociNames[i]] = i;
		}
	}
	DBG_ASSERT(m_chromNames.empty() || m_chromNames.size() == m_numLoci.size(), ValueError,
		"Chromosome names, if specified, should be given to every chromosomes");

	if (m_chromNames.empty()) {
		m_chromNames.resize(m_numLoci.size(), string());
	}
#ifndef OPTIMIZED
	else {
		map<string, int> nameMap;
		// check uniqueness of the names
		for (size_t i = 0; i < m_numLoci.size(); ++i) {
			if (!m_chromNames[i].empty() && nameMap.find(m_chromNames[i]) != nameMap.end())
				throw ValueError("Given chromosome names should be unique");
			else
				nameMap[m_chromNames[i]] = 0;
		}
	}

#endif
	// remove duplicated and empty information field names
	set<string> infoMap;
	for (size_t i = 0; i < m_infoFields.size(); ++i) {
		DBG_WARNIF(m_infoFields[i].empty(),
			"Empty information field name is removed.");
		DBG_WARNIF(infoMap.find(m_infoFields[i]) != infoMap.end(),
			"Duplicated information field name will be removed.");
		if ((!m_infoFields[i].empty()) && infoMap.find(m_infoFields[i]) == infoMap.end())
			infoMap.insert(m_infoFields[i]);
	}
	// if there are empty or duplicated info fields, re-create m_infoFields
	if (m_infoFields.size() != infoMap.size()) {
		m_infoFields.clear();
		std::copy(infoMap.begin(), infoMap.end(), std::back_inserter(m_infoFields));
	}
	// shrink allele names
	for (size_t i = 0; i < m_alleleNames.size(); ++i)
		if (!m_alleleNames[i].empty() && static_cast<UINT>(m_alleleNames[i].size() - 1) > ModuleMaxAllele)
			m_alleleNames[i].resize(ModuleMaxAllele + 1);
	//
	bool allEmpty = true;
	for (size_t i = 0; i < m_alleleNames.size(); ++i) {
		if (!m_alleleNames[i].empty()) {
			allEmpty = false;
			break;
		}
	}
	if (allEmpty)
		m_alleleNames.clear();
	//
	if (m_alleleNames.size() > 1) {
		bool allEqual = true;
		for (size_t i = 1; i < m_alleleNames.size(); ++i) {
			if (m_alleleNames[i] != m_alleleNames[0]) {
				allEqual = false;
				break;
			}
		}
		if (allEqual)
			m_alleleNames.resize(1);
	}
}


bool GenoStructure::operator==(const GenoStructure & rhs)
{
	// compare pointer directly will be fastest
	if (this == &rhs || (
	                     (m_ploidy == rhs.m_ploidy) &&
	                     (m_numLoci == rhs.m_numLoci) &&
	                     (m_chromTypes == rhs.m_chromTypes) &&
	                     (m_haplodiploid == rhs.m_haplodiploid) &&
	                     (m_lociPos == rhs.m_lociPos) &&
	                     (m_chromNames == rhs.m_chromNames) &&
	                     (m_alleleNames == rhs.m_alleleNames) &&
	                     (m_lociNames == rhs.m_lociNames) &&
	                     (m_infoFields == rhs.m_infoFields)
	                     ))
		return true;
	else
		return false;
}


void GenoStructure::setChromTypes(const vectoru & chromTypes)
{
	DBG_ASSERT(chromTypes.empty() || chromTypes.size() == m_numLoci.size(),
		ValueError, "If chromosome type is given, it should be given to all chromosomes");

	if (chromTypes.empty())
		m_chromTypes.resize(m_numLoci.size(), AUTOSOME);
	else
		m_chromTypes = chromTypes;
	// has only one chromX?
	m_chromX = -1;
	// check if the type is valid.
#ifndef OPTIMIZED
	for (size_t i = 0; i < m_chromTypes.size(); ++i) {
		size_t type = m_chromTypes[i];
		DBG_ASSERT(type == AUTOSOME || type == CHROMOSOME_X || type == CHROMOSOME_Y || type == CUSTOMIZED || type == MITOCHONDRIAL,
			ValueError, "Chromsome type can only be one of AUTOSOME, CHROMOSOME_X, CHROMOSOME_Y and CUSTOMIZED.");
	}
#endif
	for (int i = 0; i < static_cast<int>(m_chromTypes.size()); ++i) {
		if (m_chromTypes[i] == CHROMOSOME_X) {
			DBG_ASSERT(m_chromX == -1, ValueError,
				"Only one chromosome X can be specified");
			DBG_ASSERT(m_ploidy == 2, ValueError,
				"Sex chromosome can only be specified in a diploid or haplodiploid population.");
			m_chromX = i;
		}
	}
	m_chromY = -1;
	for (int i = 0; i < static_cast<int>(m_chromTypes.size()); ++i) {
		if (m_chromTypes[i] == CHROMOSOME_Y) {
			DBG_ASSERT(m_chromY == -1, ValueError,
				"Only one chromosome Y can be specified");
			DBG_ASSERT(m_ploidy == 2, ValueError,
				"Sex chromosome can only be specified in a diploid or haplodiploid population.");
			m_chromY = i;
		}
	}
	DBG_WARNIF(m_chromX * m_chromY < 0,
		"Chromosome X and Y should be both present for sexual transmission to work.");
	//
	m_mitochondrial = -1;
	for (int i = 0; i < static_cast<int>(m_chromTypes.size()); ++i) {
		if (m_chromTypes[i] == MITOCHONDRIAL) {
			DBG_ASSERT(m_mitochondrial == -1, ValueError,
				"Only one mitochondrial chromosome can be specified");
			m_mitochondrial = i;
			for (int j = i + 1; j < static_cast<int>(m_chromTypes.size()); ++j) {
				DBG_ASSERT(m_chromTypes[j] == CUSTOMIZED, ValueError,
					"Mitochondrial DNA must be specified after autosome and sex chromosomes.");
			}
		}
	}
	//
	m_customized.clear();
	for (size_t i = 0; i < m_chromTypes.size(); ++i) {
		if (m_chromTypes[i] == CUSTOMIZED) {
			DBG_FAILIF(!m_customized.empty() && m_customized.back() != i - 1,
				ValueError,
				"There can be several customized chromosmes, but they need to be adjacent to each other.");
			m_customized.push_back(static_cast<ULONG>(i));
		}
	}
}


void GenoStructure::buildLociPosMap() const
{
	if (!m_lociPosMap.empty())
		return;
	for (size_t ch = 0, loc = 0; ch < m_numLoci.size(); ++ch)
		for (size_t i = 0; i < m_numLoci[ch]; ++i, ++loc)
			m_lociPosMap[genomic_pos(m_chromNames[ch], PRECISION(m_lociPos[loc]))] = loc;
}


// initialize static variable s)genoStruRepository.
vector<GenoStructure> GenoStruTrait::s_genoStruRepository = vector<GenoStructure>();
vectorstr g_alleleName;


double GenoStruTrait::lociDist(size_t loc1, size_t loc2) const
{
	// right now, it is assumed that locus is not the first
	// on a chromosome
	DBG_FAILIF(chromLocusPair(loc1).first != chromLocusPair(loc2).first,
		ValueError, "locusDist assumes that two loci are on the same chromosome");
	return locusPos(loc2) - locusPos(loc1);
}


size_t GenoStruTrait::lociLeft(size_t loc) const
{
	CHECKRANGEABSLOCUS(loc);

	for (size_t i = 1, iEnd = numChrom(); i <= iEnd; ++i)
		if (s_genoStruRepository[m_genoStruIdx].m_chromIndex[i] > loc)
			return s_genoStruRepository[m_genoStruIdx].m_chromIndex[i] - loc;
	DBG_ASSERT(false, SystemError, "This should not be reached.");
	return 0;
}


double GenoStruTrait::distLeft(size_t loc) const
{
	CHECKRANGEABSLOCUS(loc);

	for (size_t i = 1, iEnd = numChrom(); i <= iEnd; ++i)
		if (s_genoStruRepository[m_genoStruIdx].m_chromIndex[i] > loc)
			return locusPos(s_genoStruRepository[m_genoStruIdx].m_chromIndex[i] - 1) - locusPos(loc);
	DBG_ASSERT(false, SystemError, "This should not be reached.");
	return 0;
}


size_t GenoStruTrait::lociCovered(size_t loc, double dist) const
{
	DBG_FAILIF(fcmp_lt(dist, 0.), ValueError,
		"Distance has to be positive in function lociCovered");

	const vectorf & lociPos = s_genoStruRepository[m_genoStruIdx].m_lociPos;

	size_t chrom = chromLocusPair(loc).first;
	size_t endLoc = chromEnd(chrom);
	double beginPos = lociPos[loc];

	for (size_t i = loc + 1; i < endLoc; ++i)
		if (lociPos[i] - beginPos > dist)
			return i - loc;
	return endLoc - loc;
}


void GenoStruTrait::setGenoStructure(UINT ploidy, const vectoru & loci, const vectoru & chromTypes, bool haplodiploid,
                                     const vectorf & lociPos, const vectorstr & chromNames, const matrixstr & alleleNames,
                                     const vectorstr & lociNames, const vectorstr & infoFields)
{
	GenoStructure tmp = GenoStructure(ploidy, loci, chromTypes, haplodiploid,
		lociPos, chromNames, alleleNames, lociNames, infoFields);

	setGenoStructure(tmp);
}


void GenoStruTrait::setGenoStructure(const GenoStructure & rhs)
{
	// only allow for MaxTraitIndex-1 different genotype structures
	// As a matter of fact, most simuPOP scripts have only one
	// Population type.
	if (s_genoStruRepository.size() + 1 == MaxTraitIndex) {
		throw SystemError((boost::format("This simuPOP library only allows %1%"
							             " different genotype structures. \n"
							             "If you do need more structures, modify individual.h/TraitMaxType and "
							             "recompile simuPOP.") % (MaxTraitIndex - 1)).str());
	}

	if (m_genoStruIdx != MaxTraitIndex)
		decGenoStruRef();

	for (TraitIndexType it = 0; it < s_genoStruRepository.size(); ++it) {
		// object comparison
		if (s_genoStruRepository[it] == rhs) {
			m_genoStruIdx = it;
			incGenoStruRef();
			return;
		}
	}
	// if not found, replace zero-referenced structure if necessary
	for (TraitIndexType it = 0; it < s_genoStruRepository.size(); ++it) {
		if (s_genoStruRepository[it].m_refCount == 0) {
			DBG_DO(DBG_POPULATION, cerr << "Replacing an existing geno structure." << endl);
			DBG_ASSERT(rhs.m_refCount == 0, SystemError,
				"Invalid ref count for new genotypic structure.");
			s_genoStruRepository[it] = rhs;
			m_genoStruIdx = it;
			incGenoStruRef();
			return;
		}
	}
	// no zero-referenced structure
	s_genoStruRepository.push_back(rhs);
	DBG_DO(DBG_POPULATION, cerr << "Adding an geno structure. (tot size: "
		                        << s_genoStruRepository.size() << ")" << endl);
	// the last one.
	m_genoStruIdx = static_cast<TraitIndexType>(s_genoStruRepository.size() - 1);
	// increase reference count
	incGenoStruRef();
}


const GenoStructure GenoStruTrait::gsAddChromFromStru(size_t idx) const
{
	GenoStructure & gs1 = s_genoStruRepository[m_genoStruIdx];
	GenoStructure & gs2 = s_genoStruRepository[idx];

	// identify another
	DBG_FAILIF(gs1.m_ploidy != gs2.m_ploidy || gs1.m_haplodiploid != gs2.m_haplodiploid,
		ValueError, "Added chromosome should have the same ploidy");
	//
	vectoru numLoci = gs1.m_numLoci;
	numLoci.insert(numLoci.end(), gs2.m_numLoci.begin(), gs2.m_numLoci.end());
	//
	vectorf lociPos = gs1.m_lociPos;
	lociPos.insert(lociPos.end(), gs2.m_lociPos.begin(), gs2.m_lociPos.end());
	//
	vectorstr chromNames = gs1.m_chromNames;
	chromNames.insert(chromNames.end(), gs2.m_chromNames.begin(), gs2.m_chromNames.end());
	//
	matrixstr alleleNames;
	if (!gs1.m_alleleNames.empty() || !gs2.m_alleleNames.empty()) {
		if (gs1.m_alleleNames.empty())
			alleleNames.resize(gs1.m_totNumLoci, vectorstr());
		else if (gs1.m_alleleNames.size() == 1)
			alleleNames.resize(gs1.m_totNumLoci, gs1.m_alleleNames[0]);
		else
			alleleNames = gs1.m_alleleNames;
		//
		if (gs2.m_alleleNames.empty())
			alleleNames.resize(gs1.m_totNumLoci + gs2.m_totNumLoci, vectorstr());
		else if (gs2.m_alleleNames.size() == 1)
			alleleNames.resize(gs1.m_totNumLoci + gs2.m_totNumLoci, gs2.m_alleleNames[0]);
		else
			alleleNames.insert(alleleNames.end(), gs2.m_alleleNames.begin(), gs2.m_alleleNames.end());
	}
	//
	vectorstr lociNames = gs1.m_lociNames;
	if (gs1.m_lociNames.empty()) {
		if (!gs2.m_lociNames.empty()) {
			lociNames.resize(gs1.m_totNumLoci, string());
			lociNames.insert(lociNames.end(), gs2.m_lociNames.begin(), gs2.m_lociNames.end());
		}
		// if both empty, do nothing.
	} else {
		if (gs2.m_lociNames.empty())
			lociNames.resize(gs1.m_totNumLoci + gs2.m_totNumLoci, string());
		else
			lociNames.insert(lociNames.end(), gs2.m_lociNames.begin(), gs2.m_lociNames.end());
	}
	//
	vectoru chromTypes = gs1.m_chromTypes;
	chromTypes.insert(chromTypes.end(), gs2.m_chromTypes.begin(), gs2.m_chromTypes.end());
	//
	return GenoStructure(gs1.m_ploidy, numLoci, chromTypes, gs1.m_haplodiploid, lociPos,
		chromNames, alleleNames, lociNames, gs1.m_infoFields);
}


const GenoStructure GenoStruTrait::gsAddLociFromStru(size_t idx, vectoru & index1, vectoru & index2) const
{
	GenoStructure & gs1 = s_genoStruRepository[m_genoStruIdx];
	GenoStructure & gs2 = s_genoStruRepository[idx];

	// identify another
	DBG_FAILIF(gs1.m_ploidy != gs2.m_ploidy || gs1.m_haplodiploid != gs2.m_haplodiploid,
		ValueError, "Merged population should have the same ploidy");

	// which pop has more chromosomes?
	vectoru loci(std::max(gs1.m_numLoci.size(), gs2.m_numLoci.size()));
	vectorstr chromNames;
	vectoru chromTypes;
	vectorf lociPos;
	vectorstr lociNames;
	matrixstr alleleNames;

	for (size_t ch = 0; ch < loci.size(); ++ch) {
		DBG_FAILIF(ch < gs1.m_numLoci.size() && ch < gs2.m_numLoci.size()
			&& gs1.m_chromTypes[ch] != gs2.m_chromTypes[ch],
			ValueError, "Chromosomes of different types can not be merged.");
		if (ch < gs1.m_numLoci.size() && ch < gs2.m_numLoci.size()) {
			loci[ch] = gs1.m_numLoci[ch] + gs2.m_numLoci[ch];
			chromNames.push_back(gs1.m_chromNames[ch]);
			DBG_WARNIF(gs1.m_chromNames[ch] != gs2.m_chromNames[ch],
				"Chromosome '" + gs2.m_chromNames[ch] + "' is merged to chromosome '" +
				gs1.m_chromNames[ch] + "'.");
			chromTypes.push_back(gs1.m_chromTypes[ch]);
			DBG_FAILIF(gs1.m_chromTypes[ch] != gs2.m_chromTypes[ch], ValueError,
				"Cannot merge chromosome due to conflicting chromosome type.")
			lociPos.insert(lociPos.end(), gs1.m_lociPos.begin() + gs1.m_chromIndex[ch],
				gs1.m_lociPos.begin() + gs1.m_chromIndex[ch + 1]);
			lociPos.insert(lociPos.end(), gs2.m_lociPos.begin() + gs2.m_chromIndex[ch],
				gs2.m_lociPos.begin() + gs2.m_chromIndex[ch + 1]);
			if (gs1.m_lociNames.empty()) {
				if (!gs2.m_lociNames.empty()) {
					lociNames.resize(lociNames.size() + gs1.m_numLoci[ch], string());
					lociNames.insert(lociNames.end(), gs2.m_lociNames.begin() + gs2.m_chromIndex[ch],
						gs2.m_lociNames.begin() + gs2.m_chromIndex[ch + 1]);
				}
			} else {
				lociNames.insert(lociNames.end(), gs1.m_lociNames.begin() + gs1.m_chromIndex[ch],
					gs1.m_lociNames.begin() + gs1.m_chromIndex[ch + 1]);
				if (gs2.m_lociNames.empty())
					lociNames.resize(lociNames.size() + gs2.m_numLoci[ch], string());
				else
					lociNames.insert(lociNames.end(), gs2.m_lociNames.begin() + gs2.m_chromIndex[ch],
						gs2.m_lociNames.begin() + gs2.m_chromIndex[ch + 1]);
			}
			// allele names
			if (gs1.m_alleleNames.empty())
				alleleNames.resize(alleleNames.size() + gs1.m_numLoci[ch], vectorstr());
			else if (gs1.m_alleleNames.size() == 1)
				alleleNames.resize(alleleNames.size() + gs1.m_numLoci[ch], gs1.m_alleleNames[0]);
			else
				alleleNames.insert(alleleNames.end(), gs1.m_alleleNames.begin() + gs1.m_chromIndex[ch],
					gs1.m_alleleNames.begin() + gs1.m_chromIndex[ch + 1]);
			// add from gs2
			if (gs2.m_alleleNames.empty())
				alleleNames.resize(alleleNames.size() + gs2.m_numLoci[ch], vectorstr());
			else if (gs2.m_alleleNames.size() == 1)
				alleleNames.resize(alleleNames.size() + gs2.m_numLoci[ch], gs2.m_alleleNames[0]);
			else
				alleleNames.insert(alleleNames.end(), gs2.m_alleleNames.begin() + gs2.m_chromIndex[ch],
					gs2.m_alleleNames.begin() + gs2.m_chromIndex[ch + 1]);
		} else if (ch < gs1.m_numLoci.size() && ch >= gs2.m_numLoci.size()) {     // gs1 > gs2
			loci[ch] = gs1.m_numLoci[ch];
			chromNames.push_back(gs1.m_chromNames[ch]);
			chromTypes.push_back(gs1.m_chromTypes[ch]);
			lociPos.insert(lociPos.end(), gs1.m_lociPos.begin() + gs1.m_chromIndex[ch],
				gs1.m_lociPos.begin() + gs1.m_chromIndex[ch + 1]);
			if (gs1.m_lociNames.empty()) {
				if (!lociNames.empty())
					lociNames.resize(lociNames.size() + gs1.m_numLoci[ch], string());
			} else {
				lociNames.insert(lociNames.end(), gs1.m_lociNames.begin() + gs1.m_chromIndex[ch],
					gs1.m_lociNames.begin() + gs1.m_chromIndex[ch + 1]);
			}
			// allele names
			if (gs1.m_alleleNames.empty())
				alleleNames.resize(alleleNames.size() + gs1.m_numLoci[ch], vectorstr());
			else if (gs1.m_alleleNames.size() == 1)
				alleleNames.resize(alleleNames.size() + gs1.m_numLoci[ch], gs1.m_alleleNames[0]);
			else
				alleleNames.insert(alleleNames.end(), gs1.m_alleleNames.begin() + gs1.m_chromIndex[ch],
					gs1.m_alleleNames.begin() + gs1.m_chromIndex[ch + 1]);

		} else {     // gs2 > gs1
			loci[ch] = gs2.m_numLoci[ch];
			chromNames.push_back(gs2.m_chromNames[ch]);
			chromTypes.push_back(gs2.m_chromTypes[ch]);
			lociPos.insert(lociPos.end(), gs2.m_lociPos.begin() + gs2.m_chromIndex[ch],
				gs2.m_lociPos.begin() + gs2.m_chromIndex[ch + 1]);
			if (gs2.m_lociNames.empty()) {
				if (!lociNames.empty())
					lociNames.resize(lociNames.size() + gs2.m_numLoci[ch], string());
			} else {
				lociNames.insert(lociNames.end(), gs2.m_lociNames.begin() + gs2.m_chromIndex[ch],
					gs2.m_lociNames.begin() + gs2.m_chromIndex[ch + 1]);
			}
			if (gs2.m_alleleNames.empty())
				alleleNames.resize(alleleNames.size() + gs2.m_numLoci[ch], vectorstr());
			else if (gs2.m_alleleNames.size() == 1)
				alleleNames.resize(alleleNames.size() + gs2.m_numLoci[ch], gs2.m_alleleNames[0]);
			else
				alleleNames.insert(alleleNames.end(), gs2.m_alleleNames.begin() + gs2.m_chromIndex[ch],
					gs2.m_alleleNames.begin() + gs2.m_chromIndex[ch + 1]);
		}
	}

	//
	GenoStructure ret = GenoStructure(gs1.m_ploidy, loci, chromTypes, gs1.m_haplodiploid, lociPos,
		chromNames, alleleNames, lociNames, gs1.m_infoFields);
	index1.clear();
	UINT locIdx = 0;
	for (size_t ch = 0; ch < gs1.m_numLoci.size(); ++ch) {
		for (size_t loc = 0; loc < gs1.m_numLoci[ch]; ++loc, ++locIdx) {
			double pos = gs1.m_lociPos[locIdx];
			index1.push_back(find(ret.m_lociPos.begin() + ret.m_chromIndex[ch],
					ret.m_lociPos.begin() + ret.m_chromIndex[ch + 1], pos) - ret.m_lociPos.begin());
		}
	}
	index2.clear();
	locIdx = 0;
	for (size_t ch = 0; ch < gs2.m_numLoci.size(); ++ch) {
		for (size_t loc = 0; loc < gs2.m_numLoci[ch]; ++loc, ++locIdx) {
			double pos = gs2.m_lociPos[locIdx];
			index2.push_back(find(ret.m_lociPos.begin() + ret.m_chromIndex[ch],
					ret.m_lociPos.begin() + ret.m_chromIndex[ch + 1], pos) - ret.m_lociPos.begin());
		}
	}
	return ret;
}


const GenoStructure GenoStruTrait::gsAddLociByNameFromStru(size_t idx, vectoru & index1, vectoru & index2) const
{
	GenoStructure & gs1 = s_genoStruRepository[m_genoStruIdx];
	GenoStructure & gs2 = s_genoStruRepository[idx];

	// identify another
	DBG_FAILIF(gs1.m_ploidy != gs2.m_ploidy || gs1.m_haplodiploid != gs2.m_haplodiploid,
		ValueError, "Merged population should have the same ploidy");

	// which pop has more chromosomes?
	vectoru loci(gs1.m_numLoci.size());
	vectorstr chromNames;
	vectoru chromTypes;
	vectorf lociPos;
	vectorstr lociNames;
	matrixstr alleleNames;

	// fail if there are multiple unnamed chromosomes in gs2 and gs1.
	int nUnnamed2 = 0;
	for (size_t ch2 = 0; ch2 < gs2.m_numLoci.size(); ++ch2)
		if (gs2.m_chromNames[ch2] == "")
			++nUnnamed2;
	if (nUnnamed2 > 1)
		throw ValueError("Cannot merge chromosomes by name because there are multiple ununamed chromosomes.");
	else if (nUnnamed2 == 1) {
		int nUnnamed1 = 0;
		for (size_t ch1 = 0; ch1 < gs1.m_numLoci.size(); ++ch1)
			if (gs1.m_chromNames[ch1] == "")
				++nUnnamed1;
		if (nUnnamed1 > 1)
			throw ValueError("Cannot merge chromosomes by name because an unnamed chromosome can be merged to multiple ununamed chromosomes.");
	}
	// look at chromosomes in gs1
	for (size_t ch1 = 0; ch1 < loci.size(); ++ch1) {
		// is there a matching chromosome?
		vectorstr::iterator pch2 = find(gs2.m_chromNames.begin(), gs2.m_chromNames.end(), gs1.m_chromNames[ch1]);
		// if there is no matching chromosome, use chromosome in gs1
		if (pch2 == gs2.m_chromNames.end()) {
			loci[ch1] = gs1.m_numLoci[ch1];
			chromNames.push_back(gs1.m_chromNames[ch1]);
			chromTypes.push_back(gs1.m_chromTypes[ch1]);
			lociPos.insert(lociPos.end(), gs1.m_lociPos.begin() + gs1.m_chromIndex[ch1],
				gs1.m_lociPos.begin() + gs1.m_chromIndex[ch1 + 1]);
			if (gs1.m_lociNames.empty()) {
				if (!lociNames.empty())
					lociNames.resize(lociNames.size() + gs1.m_numLoci[ch1], string());
			} else {
				lociNames.insert(lociNames.end(), gs1.m_lociNames.begin() + gs1.m_chromIndex[ch1],
					gs1.m_lociNames.begin() + gs1.m_chromIndex[ch1 + 1]);
			}
			// allele names
			if (gs1.m_alleleNames.empty())
				alleleNames.resize(alleleNames.size() + gs1.m_numLoci[ch1], vectorstr());
			else if (gs1.m_alleleNames.size() == 1)
				alleleNames.resize(alleleNames.size() + gs1.m_numLoci[ch1], gs1.m_alleleNames[0]);
			else
				alleleNames.insert(alleleNames.end(), gs1.m_alleleNames.begin() + gs1.m_chromIndex[ch1],
					gs1.m_alleleNames.begin() + gs1.m_chromIndex[ch1 + 1]);
		} else {
			// if there is matching chromosome, merge chromosomes
			size_t ch2 = pch2 - gs2.m_chromNames.begin();
			loci[ch1] = gs1.m_numLoci[ch1] + gs2.m_numLoci[ch2];
			chromNames.push_back(gs1.m_chromNames[ch1]);
			DBG_FAILIF(gs1.m_chromTypes[ch1] != gs2.m_chromTypes[ch2], ValueError,
				"Cannot merge chromosome '" + gs2.m_chromNames[ch2] + "' due to conflicting chromosome type.");
			chromTypes.push_back(gs1.m_chromTypes[ch1]);
			lociPos.insert(lociPos.end(), gs1.m_lociPos.begin() + gs1.m_chromIndex[ch1],
				gs1.m_lociPos.begin() + gs1.m_chromIndex[ch1 + 1]);
			lociPos.insert(lociPos.end(), gs2.m_lociPos.begin() + gs2.m_chromIndex[ch2],
				gs2.m_lociPos.begin() + gs2.m_chromIndex[ch2 + 1]);
			if (gs1.m_lociNames.empty()) {
				if (!gs2.m_lociNames.empty()) {
					lociNames.resize(lociNames.size() + gs1.m_numLoci[ch1], string());
					lociNames.insert(lociNames.end(), gs2.m_lociNames.begin() + gs2.m_chromIndex[ch2],
						gs2.m_lociNames.begin() + gs2.m_chromIndex[ch2 + 1]);
				}
			} else {
				lociNames.insert(lociNames.end(), gs1.m_lociNames.begin() + gs1.m_chromIndex[ch1],
					gs1.m_lociNames.begin() + gs1.m_chromIndex[ch1 + 1]);
				if (gs2.m_lociNames.empty())
					lociNames.resize(lociNames.size() + gs2.m_numLoci[ch2], string());
				else
					lociNames.insert(lociNames.end(), gs2.m_lociNames.begin() + gs2.m_chromIndex[ch2],
						gs2.m_lociNames.begin() + gs2.m_chromIndex[ch2 + 1]);
			}
			// allele names
			if (gs1.m_alleleNames.empty())
				alleleNames.resize(alleleNames.size() + gs1.m_numLoci[ch1], vectorstr());
			else if (gs1.m_alleleNames.size() == 1)
				alleleNames.resize(alleleNames.size() + gs1.m_numLoci[ch1], gs1.m_alleleNames[0]);
			else
				alleleNames.insert(alleleNames.end(), gs1.m_alleleNames.begin() + gs1.m_chromIndex[ch1],
					gs1.m_alleleNames.begin() + gs1.m_chromIndex[ch1 + 1]);
			// add from gs2
			if (gs2.m_alleleNames.empty())
				alleleNames.resize(alleleNames.size() + gs2.m_numLoci[ch2], vectorstr());
			else if (gs2.m_alleleNames.size() == 1)
				alleleNames.resize(alleleNames.size() + gs2.m_numLoci[ch2], gs2.m_alleleNames[0]);
			else
				alleleNames.insert(alleleNames.end(), gs2.m_alleleNames.begin() + gs2.m_chromIndex[ch2],
					gs2.m_alleleNames.begin() + gs2.m_chromIndex[ch2 + 1]);
		}
	}
	// now for unmapped chromosomes
	for (size_t ch2 = 0; ch2 < gs2.m_numLoci.size(); ++ch2) {
		// is there a matching chromosome?
		vectorstr::iterator pch1 = find(gs1.m_chromNames.begin(), gs1.m_chromNames.end(), gs2.m_chromNames[ch2]);
		// if there is no matching chromosome, add the chromosome from gs2
		if (pch1 == gs1.m_chromNames.end()) {
			loci.push_back(gs2.m_numLoci[ch2]);
			chromNames.push_back(gs2.m_chromNames[ch2]);
			chromTypes.push_back(gs2.m_chromTypes[ch2]);
			lociPos.insert(lociPos.end(), gs2.m_lociPos.begin() + gs2.m_chromIndex[ch2],
				gs2.m_lociPos.begin() + gs2.m_chromIndex[ch2 + 1]);
			if (gs2.m_lociNames.empty()) {
				if (!lociNames.empty())
					lociNames.resize(lociNames.size() + gs2.m_numLoci[ch2], string());
			} else {
				lociNames.insert(lociNames.end(), gs2.m_lociNames.begin() + gs2.m_chromIndex[ch2],
					gs2.m_lociNames.begin() + gs2.m_chromIndex[ch2 + 1]);
			}
			if (gs2.m_alleleNames.empty())
				alleleNames.resize(alleleNames.size() + gs2.m_numLoci[ch2], vectorstr());
			else if (gs2.m_alleleNames.size() == 1)
				alleleNames.resize(alleleNames.size() + gs2.m_numLoci[ch2], gs2.m_alleleNames[0]);
			else
				alleleNames.insert(alleleNames.end(), gs2.m_alleleNames.begin() + gs2.m_chromIndex[ch2],
					gs2.m_alleleNames.begin() + gs2.m_chromIndex[ch2 + 1]);
		}
	}

	//
	GenoStructure ret = GenoStructure(gs1.m_ploidy, loci, chromTypes, gs1.m_haplodiploid, lociPos,
		chromNames, alleleNames, lociNames, gs1.m_infoFields);
	//
	// for index 1, chromosomes have kept their original order
	index1.clear();
	UINT locIdx = 0;
	for (size_t ch = 0; ch < gs1.m_numLoci.size(); ++ch) {
		for (size_t loc = 0; loc < gs1.m_numLoci[ch]; ++loc, ++locIdx) {
			double pos = gs1.m_lociPos[locIdx];
			index1.push_back(find(ret.m_lociPos.begin() + ret.m_chromIndex[ch],
					ret.m_lociPos.begin() + ret.m_chromIndex[ch + 1], pos) - ret.m_lociPos.begin());
		}
	}
	//
	// for index 2, mapped chromosome has to be found by name
	index2.clear();
	locIdx = 0;
	for (size_t ch2 = 0; ch2 < gs2.m_numLoci.size(); ++ch2) {
		size_t ch = find(ret.m_chromNames.begin(), ret.m_chromNames.end(), gs2.m_chromNames[ch2]) -
		            ret.m_chromNames.begin();
		for (size_t loc = 0; loc < gs2.m_numLoci[ch2]; ++loc, ++locIdx) {
			double pos = gs2.m_lociPos[locIdx];
			index2.push_back(find(ret.m_lociPos.begin() + ret.m_chromIndex[ch],
					ret.m_lociPos.begin() + ret.m_chromIndex[ch + 1], pos) - ret.m_lociPos.begin());
		}
	}
	return ret;
}


const GenoStructure GenoStruTrait::gsRemoveLoci(const vectoru & kept)
{
	GenoStructure & gs = s_genoStruRepository[m_genoStruIdx];
	// loci are now remainining loci
	vectoru numLoci(numChrom(), 0);
	vectorf lociPos;
	vectorstr lociNames;
	matrixstr alleleNames;
	vectoru::const_iterator loc = kept.begin();

	for (; loc != kept.end(); ++loc) {
		size_t ch = chromLocusPair(*loc).first;
		numLoci[ch]++;
		lociPos.push_back(locusPos(*loc));
		if (!gs.m_lociNames.empty()) {
			lociNames.push_back(locusName(*loc));
			// for locus-specific allele names
			if (gs.m_alleleNames.size() > 1)
				alleleNames.push_back(gs.m_alleleNames[*loc]);
		}
	}
	// for common allele names
	if (alleleNames.empty() && gs.m_alleleNames.size() == 1)
		alleleNames = gs.m_alleleNames;
	return GenoStructure(gs.m_ploidy, numLoci, gs.m_chromTypes, isHaplodiploid(),
		lociPos, gs.m_chromNames, alleleNames, lociNames, gs.m_infoFields);
}


const GenoStructure GenoStruTrait::gsAddChrom(const vectorf & lociPos, const vectorstr & lociNames,
                                              const string & chromName, const matrixstr & alleleNames, size_t chromType) const
{
	DBG_ASSERT(lociNames.empty() || lociPos.size() == lociNames.size(), ValueError,
		"Please specify locus name for all inserted loci.");

	DBG_ASSERT(alleleNames.size() <= 1 || alleleNames.size() == lociPos.size(), ValueError,
		"Please specify allele names for none or all inserted loci.");

	for (size_t i = 1; i < lociPos.size(); ++i) {
		DBG_ASSERT(lociPos[i - 1] < lociPos[i], ValueError,
			"Loci position hsould be distinct, and in increasing order.");
	}

	GenoStructure & gs = s_genoStruRepository[m_genoStruIdx];

	// original structure
	vectoru newLoci = gs.m_numLoci;
	newLoci.push_back(lociPos.size());
	//
	vectorf newLociPos = gs.m_lociPos;
	newLociPos.insert(newLociPos.end(), lociPos.begin(), lociPos.end());
	//
	vectorstr newLociNames = gs.m_lociNames;
	if (newLociNames.empty()) {
		if (!lociNames.empty()) {
			newLociNames.resize(gs.m_totNumLoci, string());
			newLociNames.insert(newLociNames.end(), lociNames.begin(), lociNames.end());
		}
	} else {
		if (lociNames.empty())
			newLociNames.resize(newLociNames.size() + lociPos.size(), string());
		else
			newLociNames.insert(newLociNames.end(), lociNames.begin(), lociNames.end());
	}
	//
	matrixstr newAlleleNames = gs.m_alleleNames;
	if (!newAlleleNames.empty() || !alleleNames.empty()) {
		if (newAlleleNames.empty())
			newAlleleNames.resize(gs.m_totNumLoci, vectorstr());
		else if (newAlleleNames.size() == 1)
			newAlleleNames.resize(gs.m_totNumLoci, newAlleleNames[0]);
		if (alleleNames.empty())
			newAlleleNames.resize(newAlleleNames.size() + lociPos.size(), vectorstr());
		else if (alleleNames.size() == 1)
			newAlleleNames.resize(newAlleleNames.size() + lociPos.size(), alleleNames[0]);
		else
			newAlleleNames.insert(newAlleleNames.end(), alleleNames.begin(), alleleNames.end());
	}
	//
	vectorstr newChromNames = gs.m_chromNames;
	newChromNames.push_back(chromName.empty() ? (boost::format("chrom%1%") % (gs.m_numLoci.size() + 1)).str() : chromName);
	//
	vectoru newChromTypes = gs.m_chromTypes;
	newChromTypes.push_back(chromType);

	return GenoStructure(gs.m_ploidy, newLoci, newChromTypes, gs.m_haplodiploid,
		newLociPos, newChromNames, newAlleleNames, newLociNames, gs.m_infoFields);
}


const GenoStructure GenoStruTrait::gsSetAlleleNames(const lociList & loci_, const matrixstr & alleleNames)
{
	GenoStructure gs = s_genoStruRepository[m_genoStruIdx];

	// if alleleNames is used totally redefined ...
	if (loci_.allAvail())
		return GenoStructure(gs.m_ploidy, gs.m_numLoci, gs.m_chromTypes, gs.m_haplodiploid,
			gs.m_lociPos, gs.m_chromNames, alleleNames, gs.m_lociNames, gs.m_infoFields);

	const vectoru & loci = loci_.elems(this);
	matrixstr names = gs.m_alleleNames;

	// expand
	if (names.empty()) {
		names.resize(gs.m_totNumLoci, vectorstr());
		for (size_t i = 0; i < loci.size(); ++i) {
			DBG_FAILIF(loci[i] >= names.size(), IndexError, "Allele name index out of range.");
			names[loci[i]] = alleleNames[ alleleNames.size() == 1 ? 0 : i];
		}
	} else if (names.size() == 1) {
		if (gs.m_totNumLoci == 1) {
			DBG_ASSERT(alleleNames.size() == 1 && loci.size() == 1 && loci[0] == 0, SystemError,
				"Something wrong with allele names.");
			names[0] = alleleNames[0];
		} else if (alleleNames.size() != 1 || alleleNames[0] != names[0]) {
			// expand
			for (size_t i = 1; i < gs.m_totNumLoci; ++i)
				names.push_back(names[0]);
			for (size_t i = 0; i < loci.size(); ++i) {
				DBG_FAILIF(loci[i] >= names.size(), IndexError, "Allele name index out of range.");
				names[loci[i]] = alleleNames[ alleleNames.size() == 1 ? 0 : i];
			}
		}
	} else {
		DBG_ASSERT(names.size() == gs.m_totNumLoci, SystemError,
			"Inconsistent length of allele names.");
		for (size_t i = 0; i < loci.size(); ++i) {
			DBG_FAILIF(loci[i] >= names.size(), IndexError, "Allele name index out of range.");
			names[loci[i]] = alleleNames[ alleleNames.size() == 1 ? 0 : i];
		}
	}
	// replace common alleles
	return GenoStructure(gs.m_ploidy, gs.m_numLoci, gs.m_chromTypes, gs.m_haplodiploid,
		gs.m_lociPos, gs.m_chromNames, names, gs.m_lociNames, gs.m_infoFields);
}


const GenoStructure GenoStruTrait::gsAddLoci(const vectoru & chrom, const vectorf & lociPos,
                                             const vectorstr & lociNames, const matrixstr & alleleNames,
                                             vectoru & newIndex) const
{
	DBG_ASSERT(chrom.size() == lociPos.size(), ValueError,
		"Please specify chromosome and position for all inserted loci.");

	DBG_ASSERT(lociNames.empty() || lociPos.size() == lociNames.size(), ValueError,
		"Please specify locus name for none or all inserted loci.");

	DBG_ASSERT(alleleNames.empty() || alleleNames.size() == 1 || alleleNames.size() == lociPos.size(), ValueError,
		"Please speicfy allele name for none or all inserted loci.");

	GenoStructure & gs = s_genoStruRepository[m_genoStruIdx];

	// original names
	vectorstr newLociNames = gs.m_lociNames;
	if (newLociNames.empty() && !lociNames.empty())
		newLociNames.resize(gs.m_totNumLoci, string());

	matrixstr newAlleleNames = gs.m_alleleNames;
	if (newAlleleNames.empty())
		newAlleleNames.resize(gs.m_totNumLoci, vectorstr());
	else if (newAlleleNames.size() == 1)
		newAlleleNames.resize(gs.m_totNumLoci, newAlleleNames[0]);

	// the original structure...
	vectoru newLoci = gs.m_numLoci;
	vectorf newLociPos = gs.m_lociPos;
	for (size_t i = 0; i < lociPos.size(); ++i) {
		size_t ch = chrom[i];
		double pos = lociPos[i];
		string name = lociNames.empty() ? string() : lociNames[i];
		vectorstr alleleName = alleleNames.size() > 1 ? alleleNames[i] : (alleleNames.size() == 1 ? alleleNames[0] : vectorstr());
		DBG_ASSERT(ch < newLoci.size(), ValueError, "Chromosome index out of range\n"
			                                        "Please use addChrom function if a new chromosome is added");
		//
		// append to the last
		if (newLociPos.empty() || ((newLoci[ch] == 0 || pos > newLociPos.back()) && ch == numChrom() - 1)) {
			newLociPos.push_back(pos);
			if (!lociNames.empty())
				newLociNames.push_back(name);
			newAlleleNames.push_back(alleleName);
			newLoci[ch]++;
			continue;
		}

		// find beginning and end of chromosome
		size_t chBegin = 0;
		size_t chEnd = newLoci[0];
		for (size_t j = 1; j <= ch; ++j) {
			chBegin += newLoci[j - 1];
			chEnd += newLoci[j];
		}
		size_t insertPos = 0;
		// the second condition assumes that there is at least one locus on the chromosome so
		// newLociPos[chEnd-1] points to the last locus on that chromosome.
		if (newLoci[ch] == 0 || pos > newLociPos[chEnd - 1])
			insertPos = chEnd;
		else {
			for (size_t j = chBegin; j < chEnd; ++j) {
				if (pos < newLociPos[j]) {
					insertPos = j;
					break;
				}
			}
		}
		newLoci[ch]++;
		// insert here
		newLociPos.insert(newLociPos.begin() + insertPos, pos);
		if (!lociNames.empty())
			newLociNames.insert(newLociNames.begin() + insertPos, name);
		newAlleleNames.insert(newAlleleNames.begin() + insertPos, alleleName);
	}

	// set newIndex
	GenoStructure ret = GenoStructure(gs.m_ploidy, newLoci, gs.m_chromTypes, gs.m_haplodiploid,
		newLociPos, gs.m_chromNames, newAlleleNames, newLociNames, gs.m_infoFields);
	newIndex.clear();
	for (size_t i = 0; i < lociPos.size(); ++i) {
		size_t ch = chrom[i];
		double pos = lociPos[i];
		newIndex.push_back(find(ret.m_lociPos.begin() + ret.m_chromIndex[ch],
				ret.m_lociPos.begin() + ret.m_chromIndex[ch + 1], pos) - ret.m_lociPos.begin());
	}
	return ret;
}


string GenoStruTrait::ploidyName() const
{
	DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
		"PloidyName: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

	if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 1)
		return "haploid";
	else if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 2) {
		if (s_genoStruRepository[m_genoStruIdx].m_haplodiploid)
			return "haplodiploid";
		else
			return "diploid";
	} else if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 3)
		return "triploid";
	else if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 4)
		return "tetraploid";
	else
		return (boost::format("%1%-ploid") % s_genoStruRepository[m_genoStruIdx].m_ploidy).str();
}


pairu GenoStruTrait::chromLocusPair(size_t locus) const
{
	CHECKRANGEABSLOCUS(locus);

	pairu loc;

	for (size_t i = 1, iEnd = numChrom(); i <= iEnd; ++i) {
		if (s_genoStruRepository[m_genoStruIdx].m_chromIndex[i] > locus) {
			loc.first = i - 1;
			loc.second = locus - s_genoStruRepository[m_genoStruIdx].m_chromIndex[i - 1];
			break;
		}
	}
	return loc;
}


string GenoStruTrait::alleleName(const ULONG allele, const size_t locus) const
{
	DBG_FAILIF(allele > ModuleMaxAllele, IndexError,
		(boost::format("Allele %1% out of range of 0 ~ %2%") % allele % ModuleMaxAllele).str());
	if (g_alleleName.empty())
		// try to avoid excessive format call...
		for (size_t i = 0; i < 16; i++)
			g_alleleName.push_back((boost::format("%1%") % i).str());
	const matrixstr & allNames = s_genoStruRepository[m_genoStruIdx].m_alleleNames;
	if (allNames.empty())
		return allele < 16 ? g_alleleName[allele] : (boost::format("%1%") % allele).str();
	const vectorstr & names = allNames.size() > locus ? allNames[locus] : allNames[0];
	if (allele < names.size()) {
		return names[allele];
	} else
		return allele < 16 ? g_alleleName[allele] : (boost::format("%1%") % allele).str();
}


vectoru GenoStruTrait::lociByNames(const vectorstr & names) const
{
	vectoru indexes(names.size());

	const map<string, size_t> & nameMap = s_genoStruRepository[m_genoStruIdx].m_lociNameMap;

	vectorstr::const_iterator name = names.begin();
	vectorstr::const_iterator nameEnd = names.end();

	for (size_t i = 0; name != nameEnd; ++name, ++i) {
		map<string, size_t>::const_iterator it = nameMap.find(*name);

		if (it == nameMap.end())
			throw ValueError("Failed to find locus with name " + *name);

		indexes[i] = it->second;
	}

	return indexes;
}


vectoru GenoStruTrait::indexesOfLoci(const lociList & loci) const
{
	return loci.elems(this);
}


vectoru GenoStruTrait::lociByPos(const vectorpos & positions) const
{
	vectoru indexes(positions.size());

	s_genoStruRepository[m_genoStruIdx].buildLociPosMap();
	const map<genomic_pos, size_t> & lociPosMap = s_genoStruRepository[m_genoStruIdx].m_lociPosMap;

	vectorpos::const_iterator pos = positions.begin();
	vectorpos::const_iterator posEnd = positions.end();

	for (size_t i = 0; pos != posEnd; ++pos, ++i) {
		map<genomic_pos, size_t>::const_iterator it = lociPosMap.find(
			genomic_pos((*pos).first, PRECISION((*pos).second)));

		if (it == lociPosMap.end())
			throw ValueError((boost::format("Failed to find locus with chromsome %1% and position %2%")
				 % (*pos).first % size_t((*pos).second)).str());

		indexes[i] = it->second;
	}

	return indexes;
}


vectorstr GenoStruTrait::alleleNames(const size_t locus) const
{
	const matrixstr & names = s_genoStruRepository[m_genoStruIdx].m_alleleNames;

	if (names.empty())
		return vectorstr();
	return names.size() > locus ? names[locus] : names[0];
}


size_t GenoStruTrait::infoIdx(const string & name) const
{
	vectorstr & names = s_genoStruRepository[m_genoStruIdx].m_infoFields;

	for (size_t i = 0; i < names.size(); ++i)
		if (names[i] == name)
			return i;
	throw IndexError("Info field '" + name + "' is not found. Plese use infoFields=['"
		+ name + "'] option of population() during construction\n" +
		"or use addInfoField('" + name + "') to add to an existing population.");
	// this should never be reached.
	return 0;
}


const GenoStructure GenoStruTrait::gsAddInfoFields(const vectorstr & fields)
{
	GenoStructure gs = GenoStructure(s_genoStruRepository[m_genoStruIdx]);

	vectorstr::const_iterator it = fields.begin();
	vectorstr::const_iterator it_end = fields.end();

	for (; it != it_end; ++it) {
		if (std::find(gs.m_infoFields.begin(), gs.m_infoFields.end(), *it) == gs.m_infoFields.end())
			gs.m_infoFields.push_back(*it);
	}
	gs.m_refCount = 0;
	return gs;
}


const GenoStructure GenoStruTrait::gsSetInfoFields(const vectorstr & fields)
{
	GenoStructure gs = GenoStructure(s_genoStruRepository[m_genoStruIdx]);

	gs.m_infoFields.clear();
	vectorstr::const_iterator it = fields.begin();
	vectorstr::const_iterator it_end = fields.end();
	for (; it != it_end; ++it) {
		if (std::find(gs.m_infoFields.begin(), gs.m_infoFields.end(), *it) == gs.m_infoFields.end())
			gs.m_infoFields.push_back(*it);
	}
	gs.m_refCount = 0;
	return gs;
}


}
