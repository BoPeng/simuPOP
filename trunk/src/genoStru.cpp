/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu                                                        *
*                                                                         *
*   $LastChangedDate$
*   $Rev$
*
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

#include "genoStru.h"

namespace simuPOP {
GenoStructure::GenoStructure(UINT ploidy, const vectoru & loci, const vectoru & chromTypes, bool haplodiploid,
	const vectorf & lociPos, const vectorstr & chromNames, const vectorstr & alleleNames,
	const vectorstr & lociNames, const vectorstr & infoFields)
	: m_ploidy(ploidy), m_numLoci(loci), m_chromTypes(),
	m_haplodiploid(haplodiploid), m_lociPos(lociPos), m_chromIndex(loci.size() + 1),
	m_chromNames(chromNames), m_alleleNames(alleleNames), m_lociNames(lociNames),
	m_infoFields(infoFields)
{
	DBG_ASSERT(ploidy >= 1, ValueError,
		"Ploidy must be >= 1. Given " + toStr(ploidy) );

	// default: one chromosome, one locus
	// otherwise, Loci copies from loci
	if (loci.empty()) {
		m_numLoci.resize(1);
		m_numLoci[0] = 1;
		m_chromIndex.resize(2);
	}

	// chromosome type
	setChromTypes(chromTypes);

	// build chromosome index
	ULONG i, j;
	for (m_chromIndex[0] = 0, i = 1; i <= m_numLoci.size(); ++i)
		m_chromIndex[i] = m_chromIndex[i - 1] + m_numLoci[i - 1];

	m_totNumLoci = m_chromIndex[m_numLoci.size()];

	// if lociPos not specified, use 1,2,3.. 1,2,3. .. on each chromosome.
	if (m_lociPos.empty() ) {
		m_lociPos.resize(m_totNumLoci);
		for (i = 0; i < m_numLoci.size(); ++i)
			for (j = 0; j < m_numLoci[i]; j++)
				m_lociPos[m_chromIndex[i] + j] = j + 1;
	}
#ifndef OPTIMIZED
	else {                                                                            // check loci distance
		// loci distance, if specified, chould have length of chromosome.
		DBG_FAILIF(m_lociPos.size() != m_totNumLoci, ValueError,
			"You should specify loci distance for every locus (" + toStr(m_totNumLoci) + ")");

		for (i = 0; i < m_numLoci.size(); ++i)
			for (j = 0; j < m_numLoci[i]; ++j)
				DBG_FAILIF(j > 0 && fcmp_le(m_lociPos[m_chromIndex[i] + j], m_lociPos[m_chromIndex[i] + j - 1]),
					ValueError, "Loci position should be distinct, and in increasing order.");
	}
#endif

	DBG_ASSERT(m_chromNames.empty() || m_chromNames.size() == m_numLoci.size(), ValueError,
		"Chromosome names, if specified, should be given to every chromosomes");

	if (m_chromNames.empty()) {
		m_chromNames.resize(m_numLoci.size());
		for (i = 0; i < m_numLoci.size(); ++i)
			m_chromNames[i] = "chrom" + toStr(i + 1);
	}
#ifndef OPTIMIZED
	else {
		map<string, int> nameMap;
		// check uniqueness of the names
		for (i = 0; i < m_numLoci.size(); ++i) {
			if (nameMap.find(m_chromNames[i]) != nameMap.end())
				throw ValueError("Given chromosome names should be unique");
			else
				nameMap[m_chromNames[i]] = 0;
		}
	}
#endif

	DBG_ASSERT(m_lociNames.empty() || m_lociNames.size() == m_totNumLoci, ValueError,
		"Loci names, if specified, should be given to every loci");
	if (m_lociNames.empty()) {
		m_lociNames.resize(m_totNumLoci);
		for (i = 0; i < m_numLoci.size(); ++i)
			for (j = 0; j < m_numLoci[i]; j++)
				m_lociNames[m_chromIndex[i] + j] = "loc" + toStr(i + 1) + "-" + toStr(j + 1);
	}
#ifndef OPTIMIZED
	else {
		map<string, int> nameMap;
		// check uniqueness of the names
		for (i = 0; i < m_totNumLoci; ++i) {
			if (nameMap.find(m_lociNames[i]) != nameMap.end())
				throw ValueError("Given loci names should be unique");
			else
				nameMap[m_lociNames[i]] = 0;
		}
	}
#endif
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
		m_chromTypes.resize(m_numLoci.size(), Autosome);
	else
		m_chromTypes = chromTypes;
	// has only one chromX?
	m_chromX = -1;
	// check if the type is valid.
	for (size_t i = 0; i < m_chromTypes.size(); ++i) {
		UINT type = m_chromTypes[i];
		DBG_ASSERT(type == Autosome || type == ChromosomeX || type == ChromosomeY || type == Mitochondrial,
			ValueError, "Chromsome type can only be one of Autosome, ChromosomeX, ChromosomeY and Mitochondrial");
	}
	for (size_t i = 0; i < m_chromTypes.size(); ++i) {
		if (m_chromTypes[i] == ChromosomeX) {
			DBG_ASSERT(m_chromX == -1, ValueError,
				"Only one chromosome X can be specified");
			m_chromX = i;
		}
	}
	m_chromY = -1;
	for (size_t i = 0; i < m_chromTypes.size(); ++i) {
		if (m_chromTypes[i] == ChromosomeY) {
			DBG_ASSERT(m_chromY == -1, ValueError,
				"Only one chromosome Y can be specified");
			m_chromY = i;
		}
	}
	DBG_FAILIF(m_chromX * m_chromY < 0, ValueError,
		"It is invalid to set only chromosome X or Y.");
	//
	m_mitochondrial = -1;
	for (size_t i = 0; i < m_chromTypes.size(); ++i) {
		if (m_chromTypes[i] == Mitochondrial) {
			DBG_ASSERT(m_mitochondrial == -1, ValueError,
				"Only one mitochondria chromosome can be specified");
			m_mitochondrial = i;
		}
	}
}


// initialize static variable s)genoStruRepository.
vector<GenoStructure> GenoStruTrait::s_genoStruRepository = vector<GenoStructure>();


double GenoStruTrait::lociDist(UINT loc1, UINT loc2) const
{
	// right now, it is assumed that locus is not the first
	// on a chromosome
	DBG_FAILIF(chromLocusPair(loc1).first != chromLocusPair(loc2).first,
		ValueError, "locusDist assumes that two loci are on the same chromosome");
	return locusPos(loc2) - locusPos(loc1);
}


UINT GenoStruTrait::lociLeft(UINT loc) const
{
	CHECKRANGEABSLOCUS(loc);

	for (UINT i = 1, iEnd = numChrom(); i <= iEnd;  ++i)
		if (s_genoStruRepository[m_genoStruIdx].m_chromIndex[i] > loc)
			return s_genoStruRepository[m_genoStruIdx].m_chromIndex[i] - loc;
	DBG_ASSERT(false, SystemError, "This should not be reached.");
	return 0;
}


double GenoStruTrait::distLeft(UINT loc) const
{
	CHECKRANGEABSLOCUS(loc);

	for (UINT i = 1, iEnd = numChrom(); i <= iEnd;  ++i)
		if (s_genoStruRepository[m_genoStruIdx].m_chromIndex[i] > loc)
			return locusPos(s_genoStruRepository[m_genoStruIdx].m_chromIndex[i] - 1) - locusPos(loc);
	DBG_ASSERT(false, SystemError, "This should not be reached.");
	return 0;
}


UINT GenoStruTrait::lociCovered(UINT loc, double dist) const
{
	DBG_FAILIF(fcmp_lt(dist, 0.), ValueError,
		"Distance has to be positive in function lociCovered");

	const vectorf & lociPos = s_genoStruRepository[m_genoStruIdx].m_lociPos;

	int chrom = chromLocusPair(loc).first;
	UINT endLoc = chromEnd(chrom);
	double beginPos = lociPos[loc];

	for (size_t i = loc + 1; i < endLoc; ++i)
		if (lociPos[i] - beginPos > dist)
			return i - loc;
	return endLoc - loc;
}


void GenoStruTrait::setGenoStructure(UINT ploidy, const vectoru & loci, const vectoru & chromTypes, bool haplodiploid,
                                     const vectorf & lociPos, const vectorstr & chromNames, const vectorstr & alleleNames,
                                     const vectorstr & lociNames, const vectorstr & infoFields)
{
	// only allow for TraitMaxIndex-1 different genotype structures
	// As a matter of fact, most simuPOP scripts have only one
	// population type.
	if (s_genoStruRepository.size() == TraitMaxIndex - 1) {
		throw SystemError("This simuPOP library only allows " + toStr(TraitMaxIndex - 1)
			+ " different genotype structures. \n" +
			+ "If you do need more structures, modify individual.h/TraitMaxType and " +
			+ "recompile simuPOP.");
	}

	GenoStructure tmp = GenoStructure(ploidy, loci, chromTypes, haplodiploid,
		lociPos, chromNames, alleleNames, lociNames, infoFields);

	for (TraitIndexType it = 0; it < s_genoStruRepository.size();
	     ++it) {
		// object comparison
		if (s_genoStruRepository[it] == tmp) {
			m_genoStruIdx = it;
			return;
		}
	}
	// if not found
	s_genoStruRepository.push_back(tmp);
	// the last one.
	m_genoStruIdx = s_genoStruRepository.size() - 1;
}


GenoStructure & GenoStruTrait::gsAddChromFromStru(size_t idx) const
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
	vectorstr lociNames = gs1.m_lociNames;
	lociNames.insert(lociNames.end(), gs2.m_lociNames.begin(), gs2.m_lociNames.end());
	//
	vectoru chromTypes = gs1.m_chromTypes;
	chromTypes.insert(chromTypes.end(), gs2.m_chromTypes.begin(), gs2.m_chromTypes.end());
	//
	return *new GenoStructure(gs1.m_ploidy, numLoci, chromTypes, gs1.m_haplodiploid, lociPos,
		chromNames, gs1.m_alleleNames, lociNames, gs1.m_infoFields);
}


GenoStructure & GenoStruTrait::gsAddLociFromStru(size_t idx) const
{
#define addLocusName(name); \
	if (std::find(lociNames.begin(), lociNames.end(), name) == lociNames.end()) \
		lociNames.push_back(name);\
	else \
	{ \
		int n = 1; \
		while (true) \
		{ \
			string name_ = name + "_" + toStr(n++); \
			if (std::find(lociNames.begin(), lociNames.end(), name_) == lociNames.end()) \
			{ \
				lociNames.push_back(name_); \
				break; \
			} \
		} \
	}


	GenoStructure & gs1 = s_genoStruRepository[m_genoStruIdx];
	GenoStructure & gs2 = s_genoStruRepository[idx];

	// identify another
	DBG_FAILIF(gs1.m_ploidy != gs2.m_ploidy || gs1.m_haplodiploid != gs2.m_haplodiploid,
		ValueError, "Merged population should have the same ploidy");

	// which pop has more chromosomes?
	vectoru loci(std::max(gs1.m_numLoci.size(), gs2.m_numLoci.size()));
	vectorstr chromNames;
	vectoru chromTypes;
	for (size_t ch = 0; ch < loci.size(); ++ch) {
		DBG_FAILIF(ch < gs1.m_numLoci.size() && ch < gs2.m_numLoci.size()
			&& gs1.m_chromTypes[ch] != gs2.m_chromTypes[ch],
			ValueError, "Chromosomes of different types can not be merged.");
		if (ch < gs1.m_numLoci.size()) {
			loci[ch] += gs1.m_numLoci[ch];
			chromNames.push_back(gs1.m_chromNames[ch]);
			chromTypes.push_back(gs1.m_chromTypes[ch]);
		}
		if (ch < gs2.m_numLoci.size()) {
			loci[ch] += gs2.m_numLoci[ch];
			chromNames.push_back(gs2.m_chromNames[ch]);
			chromTypes.push_back(gs2.m_chromTypes[ch]);
		}
	}

	// loci pos and loci name
	vectorf lociPos;
	vectorstr lociNames;

	for (size_t ch = 0; ch < loci.size(); ++ch) {
		size_t idx1 = 0;
		size_t idx2 = 0;
		for (size_t loc = 0; loc < loci[ch]; ++loc) {
			if (ch < gs1.m_numLoci.size() && ch < gs2.m_numLoci.size()) {
				// index 1 done
				double pos1 = idx1 < gs1.m_numLoci[ch] ? gs1.m_lociPos[gs1.m_chromIndex[ch] + idx1] : 1.e9;
				double pos2 = idx2 < gs2.m_numLoci[ch] ? gs2.m_lociPos[gs2.m_chromIndex[ch] + idx2] : 1.e9;
				if (idx2 >= gs2.m_numLoci[ch] || pos1 < pos2) {
					// push this in
					lociPos.push_back(pos1);
					string name = gs1.m_lociNames[gs1.m_chromIndex[ch] + idx1];
					addLocusName(name);
					idx1++;
				} else if (idx1 >= gs1.m_numLoci[ch] || pos1 > pos2) {
					// push this in
					lociPos.push_back(pos2);
					string name = gs2.m_lociNames[gs2.m_chromIndex[ch] + idx2];
					addLocusName(name);
					idx2++;
				} else
					throw ValueError("Duplicate loci position. Can not merge");
			} else if (ch < gs1.m_numLoci.size() && ch >= gs2.m_numLoci.size()) {
				// add idx 1
				lociPos.push_back(gs1.m_lociPos[gs1.m_chromIndex[ch] + idx1]);
				string name = gs1.m_lociNames[gs1.m_chromIndex[ch] + idx1];
				addLocusName(name);
				idx1++;
			} else if (ch >= gs1.m_numLoci.size() && ch < gs2.m_numLoci.size()) {
				// add idx 2
				lociPos.push_back(gs2.m_lociPos[gs2.m_chromIndex[ch] + idx2]);
				string name = gs2.m_lociNames[gs2.m_chromIndex[ch] + idx2];
				addLocusName(name);
				idx2++;
			} else
				throw SystemError("This should not happen");
		}
	}
	//
	return *new GenoStructure(gs1.m_ploidy, loci, chromTypes, gs1.m_haplodiploid, lociPos,
		chromNames, gs1.m_alleleNames, lociNames, gs1.m_infoFields);
#undef addLocusName
}


GenoStructure & GenoStruTrait::gsRemoveLoci(const vectoru & loci,
                                            vectoru & kept)
{
	if (kept.empty()) {
		for (size_t loc = 0; loc < totNumLoci(); ++loc) {
			if (find(loci.begin(), loci.end(), loc) == loci.end())
				kept.push_back(loc);
		}
	}

	// loci are now remainining loci
	vectoru numLoci(numChrom(), 0);
	vectorf lociPos;
	vectorstr lociNames;
	vectoru::iterator loc = kept.begin();
	for (; loc != kept.end(); ++loc) {
		UINT ch = chromLocusPair(*loc).first;
		numLoci[ch]++;
		lociPos.push_back(locusPos(*loc));
		lociNames.push_back(locusName(*loc));
	}
	return *new GenoStructure(ploidy(), numLoci, chromTypes(), isHaplodiploid(),
		lociPos, chromNames(), alleleNames(), lociNames, infoFields());
}


GenoStructure & GenoStruTrait::gsAddChrom(const vectorf & lociPos, const vectorstr & lociNames,
                                          const string & chromName, UINT chromType) const
{
	DBG_ASSERT(lociNames.empty() || lociPos.size() == lociNames.size(), ValueError,
		"Please specify locus name for all inserted loci.");

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
	if (lociNames.empty()) {
		for (size_t i = 0; i < lociPos.size(); ++i)
			newLociNames.push_back("loc" + toStr(gs.m_numLoci.size() + 1) + "-" + toStr(i + 1));
	} else
		newLociNames.insert(newLociNames.end(), lociNames.begin(), lociNames.end());
	//
	vectorstr newChromNames = gs.m_chromNames;
	newChromNames.push_back(chromName.empty() ? "chrom" + toStr(gs.m_numLoci.size() + 1) : chromName);
	//
	vectoru newChromTypes = gs.m_chromTypes;
	newChromTypes.push_back(chromType);

	return *new GenoStructure(gs.m_ploidy, newLoci, newChromTypes, gs.m_haplodiploid,
		newLociPos, newChromNames, gs.m_alleleNames, newLociNames, gs.m_infoFields);
}


GenoStructure & GenoStruTrait::gsAddLoci(const vectoru & chrom, const vectorf & lociPos,
                                         const vectorstr & lociNames, vectoru & newIndex) const
{
	DBG_ASSERT(chrom.size() == lociPos.size(), ValueError,
		"Please specify chromosome and position for all inserted loci.");

	DBG_ASSERT(lociNames.empty() || lociPos.size() == lociNames.size(), ValueError,
		"Please specify locus name for all inserted loci.");

	GenoStructure & gs = s_genoStruRepository[m_genoStruIdx];

	// original names
	vectorstr newLociNames = gs.m_lociNames;
	// first, new names
	vectorstr newNames = lociNames;
	if (lociNames.empty()) {
		for (size_t i = 0, j = 0; i < lociPos.size(); ++i) {
			while (true) {
				++j;
				string name = "ins" + toStr(j);
				if (std::find(newLociNames.begin(), newLociNames.end(), name) == newLociNames.end()) {
					newNames.push_back(name);
					break;
				}
			}
		}
	}

	// the original structure...
	vectoru newLoci = gs.m_numLoci;
	vectorf newLociPos = gs.m_lociPos;
	for (size_t i = 0; i < lociPos.size(); ++i) {
		UINT ch = chrom[i];
		double pos = lociPos[i];
		string name = newNames[i];
		DBG_ASSERT(ch < newLoci.size(), ValueError, "Chromosome index out of range\n"
			                                        "Please use addChrom function if a new chromosome is added");
		//
		newLoci[ch]++;
		// find beginning and end of chromosome
		size_t chBegin = 0;
		size_t chEnd = newLoci[0];
		for (size_t i = 1; i <= ch; ++i) {
			chBegin += newLoci[i - 1];
			chEnd += newLoci[i];
		}
		// append to the last
		if (pos > lociPos[chEnd - 1] && ch == numChrom() - 1) {
			newLociPos.push_back(pos);
			newLociNames.push_back(name);
			continue;
		}
		size_t insertPos = 0;
		if (pos > lociPos[chEnd - 1])
			insertPos = chEnd;
		else {
			for (size_t i = chBegin; i < chEnd; ++i) {
				if (pos < newLociPos[i]) {
					insertPos = i;
					break;
				}
			}
		}
		// insert here
		newLociPos.insert(newLociPos.begin() + insertPos, pos);
		newLociNames.insert(newLociNames.begin() + insertPos, name);
	}

	// set newIndex
	newIndex.clear();
	for (vectorstr::const_iterator name = newNames.begin();
	     name != newNames.end(); ++name)
		newIndex.push_back(find(newLociNames.begin(), newLociNames.end(), *name)
			- newLociNames.begin());
	return *new GenoStructure(gs.m_ploidy, newLoci, gs.m_chromTypes, gs.m_haplodiploid,
		newLociPos, gs.m_chromNames, gs.m_alleleNames, newLociNames, gs.m_infoFields);
}


void GenoStruTrait::setGenoStructure(GenoStructure & rhs)
{
	for (TraitIndexType it = 0; it < s_genoStruRepository.size();
	     ++it) {
		// object comparison
		if (s_genoStruRepository[it] == rhs) {
			m_genoStruIdx = it;
			return;
		}
	}

	// if not found, make a copy and store it.
	s_genoStruRepository.push_back(rhs);
	m_genoStruIdx = s_genoStruRepository.size() - 1;
}


string GenoStruTrait::ploidyName() const
{
	DBG_FAILIF(m_genoStruIdx == TraitMaxIndex, SystemError,
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
		return toStr(s_genoStruRepository[m_genoStruIdx].m_ploidy) + "-polid";
}


std::pair<UINT, UINT> GenoStruTrait::chromLocusPair(UINT locus) const
{
	CHECKRANGEABSLOCUS(locus);

	pair<UINT, UINT> loc;

	for (UINT i = 1, iEnd = numChrom(); i <= iEnd;  ++i) {
		if (s_genoStruRepository[m_genoStruIdx].m_chromIndex[i] > locus) {
			loc.first = i - 1;
			loc.second = locus - s_genoStruRepository[m_genoStruIdx].m_chromIndex[i - 1];
			break;
		}
	}
	return loc;
}


string GenoStruTrait::alleleName(const UINT allele) const
{
	DBG_FAILIF(allele > ModuleMaxAllele,
		IndexError, "Allele " + toStr(allele) + " out of range of 0 ~ " +
		toStr(ModuleMaxAllele));
	if (allele < s_genoStruRepository[m_genoStruIdx].m_alleleNames.size() ) {
		DBG_FAILIF(allele >= s_genoStruRepository[m_genoStruIdx].m_alleleNames.size(),
			IndexError, "No name for allele " + toStr(allele));
		return s_genoStruRepository[m_genoStruIdx].m_alleleNames[allele];
	} else
		return toStr(allele);
}


UINT GenoStruTrait::infoIdx(const string & name) const
{
	vectorstr & names = s_genoStruRepository[m_genoStruIdx].m_infoFields;

	for (UINT i = 0; i < names.size(); ++i)
		if (names[i] == name)
			return i;
	throw IndexError("Info field '" + name + "' is not found. Plese use infoFields=['"
		+ name + "'] option of population() during construction\n" +
		"or use addInfoField('" + name + "') to add to an existing population.");
	// this should never be reached.
	return 0;
}


GenoStructure & GenoStruTrait::struAddInfoFields(const vectorstr & fields)
{
	GenoStructure * gs = new GenoStructure(s_genoStruRepository[m_genoStruIdx]);

	gs->m_infoFields.insert(gs->m_infoFields.end(), fields.begin(), fields.end());
	return *gs;
}


GenoStructure & GenoStruTrait::struSetInfoFields(const vectorstr & fields)
{
	GenoStructure * gs = new GenoStructure(s_genoStruRepository[m_genoStruIdx]);

	gs->m_infoFields = fields;
	return *gs;
}


}
