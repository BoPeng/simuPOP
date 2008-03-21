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
GenoStructure::GenoStructure(UINT ploidy, const vectoru & loci, bool sexChrom,
                             const vectorf & lociPos, const vectorstr & chromNames, const vectorstr & alleleNames,
                             const vectorstr & lociNames, UINT maxAllele, const vectorstr & infoFields)
	: m_ploidy(ploidy),  m_numChrom(loci.size()), m_numLoci(loci), m_sexChrom(sexChrom),
	m_lociPos(lociPos), m_chromIndex(loci.size() + 1),
	m_chromNames(chromNames), m_alleleNames(alleleNames), m_lociNames(lociNames),
	m_maxAllele(maxAllele),
	m_infoFields(infoFields)
{
#ifdef BINARYALLELE
	DBG_ASSERT(maxAllele == 1, ValueError,
		"max allele must be 1 for binary modules");
#endif
	DBG_ASSERT(ploidy >= 1, ValueError,
		"Ploidy must be >= 1. Given " + toStr(ploidy) );

	// default: one chromosome, one locus
	// otherwise, Loci copies from loci
	if (loci.empty()) {
		m_numChrom = 1;
		m_numLoci.resize(1);
		m_numLoci[0] = 1;
		m_chromIndex.resize(2);
	}

	// build chromosome index
	ULONG i, j;
	for (m_chromIndex[0] = 0, i = 1; i <= m_numChrom; ++i)
		m_chromIndex[i] = m_chromIndex[i - 1] + m_numLoci[i - 1];

	m_totNumLoci = m_chromIndex[m_numChrom];
	m_genoSize = m_totNumLoci * m_ploidy;

	// if lociPos not specified, use 1,2,3.. 1,2,3. .. on each chromosome.
	if (m_lociPos.empty() ) {
		m_lociPos.resize(m_totNumLoci);
		for (i = 0; i < m_numChrom; ++i)
			for (j = 0; j < m_numLoci[i]; j++)
				m_lociPos[m_chromIndex[i] + j] = j + 1;
	}
#ifndef OPTIMIZED
	else {                                                                            // check loci distance
		// loci distance, if specified, chould have length of chromosome.
		DBG_FAILIF(m_lociPos.size() != m_totNumLoci, ValueError,
			"You should specify loci distance for every locus (" + toStr(m_totNumLoci) + ")");

		for (i = 0; i < m_numChrom; ++i)
			for (j = 0; j < m_numLoci[i]; ++j)
				DBG_FAILIF(j > 0 && fcmp_lt(m_lociPos[m_chromIndex[i] + j], m_lociPos[m_chromIndex[i] + j - 1]),
					ValueError, "Loci distance must be in order.");
	}
#endif

	DBG_ASSERT(m_chromNames.empty() || m_chromNames.size() == m_numChrom, ValueError,
		"Chromosome names, if specified, should be given to every chromosomes");

	if (m_chromNames.empty()) {
		m_chromNames.resize(m_numChrom);
		for (i = 0; i < m_numChrom; ++i)
			m_chromNames[i] = "chrom" + toStr(i + 1);
	}

	DBG_ASSERT(m_lociNames.empty() || m_lociNames.size() == m_totNumLoci, ValueError,
		"Loci names, if specified, should be given to every loci");
	if (m_lociNames.empty()) {
		m_lociNames.resize(m_totNumLoci);
		for (i = 0; i < m_numChrom; ++i)
			for (j = 0; j < m_numLoci[i]; j++)
				m_lociNames[m_chromIndex[i] + j] = "loc" + toStr(i + 1) + "-" + toStr(j + 1);
	}
#ifndef OPTIMIZED
	else {
		// check uniqueness of the names
		for (i = 0; i < m_totNumLoci; ++i)
			for (j = i + 1; j < m_totNumLoci; ++j)
				if (m_lociNames[i] == m_lociNames[j])
					throw ValueError("Given loci names should be unique");
	}
#endif

	DBG_WARNING( (!m_alleleNames.empty()) && m_alleleNames.size() != m_maxAllele + 1,
		"Not all allele names are given. ");

}


bool GenoStructure::operator==(const GenoStructure & rhs)
{
	// compare pointer directly will be fastest
	if (this == &rhs || (
	                     (m_ploidy == rhs.m_ploidy) &&
	                     (m_numLoci == rhs.m_numLoci) &&
	                     (m_sexChrom == rhs.m_sexChrom) &&
	                     (m_lociPos == rhs.m_lociPos) &&
	                     (m_chromNames == rhs.m_chromNames) &&
	                     (m_alleleNames == rhs.m_alleleNames) &&
	                     (m_lociNames == rhs.m_lociNames) &&
	                     (m_maxAllele == rhs.m_maxAllele) &&
	                     (m_infoFields == rhs.m_infoFields)
	                     ))
		return true;
	else
		return false;
}


// initialize static variable s)genoStruRepository.
vector<GenoStructure> GenoStruTrait::s_genoStruRepository = vector<GenoStructure>();


/// distance between loci \c loc1 and \c loc2. These two loci should be
/// on the same chromosome.
double GenoStruTrait::lociDist(UINT loc1, UINT loc2) const
{
	// right now, it is assumed that locus is not the first
	// on a chromosome
	DBG_FAILIF(chromLocusPair(loc1).first != chromLocusPair(loc2).first,
		ValueError, "locusDist assumes that two loci are on the same chromosome");
	return locusPos(loc2) - locusPos(loc1);
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


void GenoStruTrait::setGenoStructure(UINT ploidy, const vectoru & loci, bool sexChrom,
                                     const vectorf & lociPos, const vectorstr & chromNames, const vectorstr & alleleNames,
                                     const vectorstr & lociNames, UINT maxAllele, const vectorstr & infoFields)
{
	/// only allow for TraitMaxIndex-1 different genotype structures
	/// As a matter of fact, most simuPOP scripts have only one
	/// population type.
	if (s_genoStruRepository.size() == TraitMaxIndex - 1) {
		throw SystemError("This simuPOP library only allows " + toStr(TraitMaxIndex - 1)
			+ " different genotype structures. \n" +
			+ "If you do need more structures, modify individual.h/TraitMaxType and " +
			+ "recompile simuPOP.");
	}

	GenoStructure tmp = GenoStructure(ploidy, loci, sexChrom,
		lociPos, chromNames, alleleNames, lociNames, maxAllele, infoFields);

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


GenoStructure & GenoStruTrait::mergeGenoStru(size_t idx, bool byChromosome) const
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
	DBG_FAILIF(gs1.m_ploidy != gs2.m_ploidy, ValueError,
		"Merged population should have the same ploidy");
	// which pop has more chromosomes?
	if (byChromosome) {
		// loci
		vectoru loci(std::max(gs1.m_numLoci.size(), gs2.m_numLoci.size()));
		for (size_t ch = 0; ch < loci.size(); ++ch) {
			if (ch < gs1.m_numLoci.size())
				loci[ch] += gs1.m_numLoci[ch];
			if (ch < gs2.m_numLoci.size())
				loci[ch] += gs2.m_numLoci[ch];
		}
		vectorstr chromNames;
		for (size_t ch = 0; ch < loci.size(); ++ch) {
			if (ch < gs1.m_numLoci.size())
				chromNames.push_back(gs1.m_chromNames[ch]);
			else
				chromNames.push_back(gs2.m_chromNames[ch]);
		}
		DBG_DO(DBG_POPULATION, cout << "New number of loci " << loci << endl);
		// sex chrom
		DBG_FAILIF(gs1.m_sexChrom && gs2.m_sexChrom && gs1.m_numLoci.size() != gs2.m_numLoci.size(),
			ValueError, "If both population has sex chromosome, they should be the same chromosome");
		DBG_FAILIF(gs1.m_sexChrom && !gs2.m_sexChrom && gs1.m_numLoci.size() <= gs2.m_numLoci.size(),
			ValueError, "The same chromosome of another population should also be the sex chromosome");
		DBG_FAILIF(!gs1.m_sexChrom && gs2.m_sexChrom && gs2.m_numLoci.size() <= gs1.m_numLoci.size(),
			ValueError, "The same chromosome of another population should also be the sex chromosome");
		bool sexChrom = gs1.m_sexChrom || gs2.m_sexChrom;
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
		}                                                                                 // each chromosome
		DBG_DO(DBG_POPULATION, cout << "New loci positions: " << lociPos << endl);
		DBG_DO(DBG_POPULATION, cout << "New loci names: " << lociNames << endl);
		//
		UINT maxAllele = std::max(gs1.m_maxAllele, gs2.m_maxAllele);
		return *new GenoStructure(gs1.m_ploidy, loci, sexChrom, lociPos,
			chromNames, gs1.m_alleleNames, lociNames, maxAllele, gs1.m_infoFields);
	} else {
		vectoru loci = gs1.m_numLoci;
		loci.insert(loci.end(), gs2.m_numLoci.begin(), gs2.m_numLoci.end());
		DBG_FAILIF(gs1.m_sexChrom, ValueError,
			"Population with sex chromosome has to be at the end");
		vectorf lociPos = gs1.m_lociPos;
		lociPos.insert(lociPos.end(), gs2.m_lociPos.begin(), gs2.m_lociPos.end());
		DBG_FAILIF(gs1.m_alleleNames != gs2.m_alleleNames, ValueError,
			"Merged population should have the same allele names (sorry, no allele names at each locus for now)");
		vectorstr lociNames = gs1.m_lociNames;
		vectorstr chromNames = gs1.m_chromNames;
		chromNames.insert(chromNames.end(), gs2.m_chromNames.begin(), gs2.m_chromNames.end());
		// add locus name, if there is no duplicate, fine. Otherwise, add '_' to the names.
		for (vectorstr::const_iterator it = gs2.m_lociNames.begin(); it != gs2.m_lociNames.end(); ++it) {
			addLocusName(*it);
		}
		UINT maxAllele = std::max(gs1.m_maxAllele, gs2.m_maxAllele);
		return *new GenoStructure(gs1.m_ploidy, loci, gs2.m_sexChrom, lociPos,
			chromNames, gs1.m_alleleNames, lociNames, maxAllele, gs1.m_infoFields);
	}
#undef addLocusName
}

/// CPPONLY
GenoStructure & GenoStruTrait::removeLociFromGenoStru(const vectoru & remove, const vectoru & keep)
{
	vectoru loci;

	DBG_FAILIF(remove.empty() && keep.empty(), ValueError,
		"Please specify either remove or keep");
	if (remove.empty())
		loci = keep;
	else {
		for (size_t loc = 0; loc < totNumLoci(); ++loc)
			// if not removed
			if (find(remove.begin(), remove.end(), loc) == remove.end())
				loci.push_back(loc);
	}
	// loci are now remainining loci
	vectoru newNumLoci;
	vectorf newLociDist;
	vectorstr newLociNames;
	vectorstr newChromNames;
	UINT curCh = 999999;                                              // not 0, will be set to 0 soon.
	for (vectoru::iterator loc = loci.begin();
	     loc != loci.end(); ++loc) {
		UINT ch = chromLocusPair(*loc).first;
		if (newNumLoci.empty() || curCh != ch) {
			newNumLoci.push_back(1);
			newChromNames.push_back(chromName(ch));
			curCh = ch;
		} else
			newNumLoci.back()++;
		newLociDist.push_back(locusPos(*loc));
		newLociNames.push_back(locusName(*loc));
	}
	return *new GenoStructure(ploidy(), newNumLoci, sexChrom(), newLociDist,
		newChromNames, alleleNames(), newLociNames, maxAllele(), infoFields());
}


/// CPPONLY add some loci to genotype structure
GenoStructure & GenoStruTrait::insertBeforeLociToGenoStru(const vectoru & idx, const vectorf & pos, const vectorstr & names) const
{
	DBG_FAILIF(idx.size() != pos.size(), ValueError,
		"Index and position should have the same length");
	DBG_FAILIF(!names.empty() && idx.size() != names.size(), ValueError,
		"LociNames, if given, should be the same length as idx");

	GenoStructure & gs = s_genoStruRepository[m_genoStruIdx];

	vectoru loci = gs.m_numLoci;
	vectorf lociPos = gs.m_lociPos;
	for (size_t i = 0; i < idx.size(); ++i) {
		CHECKRANGEABSLOCUS(idx[i]);
		lociPos.insert(lociPos.begin() + idx[i] + i, pos[i]);
		loci[chromLocusPair(idx[i]).first]++;
	}
	vectorstr lociNames = gs.m_lociNames;
	// add locus name, if there is no duplicate, fine. Otherwise, add '_' to the names.
	if (names.empty()) {
		for (size_t i = 0; i < idx.size(); ++i) {
			int ch = chromLocusPair(idx[i]).first + 1;
			int loc = chromLocusPair(idx[i]).second + 1;
			string base_name = "ins" + toStr(ch) + "_" + toStr(loc) + "_";
			int n = 1;
			while (true) {
				string name = base_name + toStr(n++);
				if (std::find(lociNames.begin(), lociNames.end(), name) == lociNames.end()) {
					// +i because the array has shifted.
					lociNames.insert(lociNames.begin() + idx[i] + i, name);
					break;
				}
			}
		}
	} else {                                                                          // given names
		for (size_t i = 0; i < idx.size(); ++i) {
			if (std::find(lociNames.begin(), lociNames.end(), names[i]) == lociNames.end())
				lociNames.insert(lociNames.begin() + idx[i] + i, names[i]);
			else
				throw ValueError("Locus name " + names[i] + " already exists.");
		}
	}
	return *new GenoStructure(gs.m_ploidy, loci, gs.m_sexChrom, lociPos,
		gs.m_chromNames, gs.m_alleleNames, lociNames, gs.m_maxAllele, gs.m_infoFields);
}


/// CPPONLY append some loci to genotype structure
GenoStructure & GenoStruTrait::insertAfterLociToGenoStru(const vectoru & idx, const vectorf & pos, const vectorstr & names) const
{
	DBG_FAILIF(idx.size() != pos.size(), ValueError,
		"Index and position should have the same length");
	DBG_FAILIF(!names.empty() && idx.size() != names.size(), ValueError,
		"LociNames, if given, should be the same length as idx");

	GenoStructure & gs = s_genoStruRepository[m_genoStruIdx];

	vectoru loci = gs.m_numLoci;
	vectorf lociPos = gs.m_lociPos;
	for (size_t i = 0; i < idx.size(); ++i) {
		CHECKRANGEABSLOCUS(idx[i]);
		lociPos.insert(lociPos.begin() + idx[i] + 1 + i, pos[i]);
		loci[chromLocusPair(idx[i]).first]++;
	}
	vectorstr lociNames = gs.m_lociNames;
	// add locus name, if there is no duplicate, fine. Otherwise, add '_' to the names.
	if (names.empty()) {
		for (size_t i = 0; i < idx.size(); ++i) {
			int ch = chromLocusPair(idx[i]).first + 1;
			int loc = chromLocusPair(idx[i]).second + 1;
			string base_name = "app" + toStr(ch) + "_" + toStr(loc) + "_";
			int n = 1;
			while (true) {
				string name = base_name + toStr(n++);
				if (std::find(lociNames.begin(), lociNames.end(), name) == lociNames.end()) {
					lociNames.insert(lociNames.begin() + idx[i] + i + 1, name);
					break;
				}
			}
		}
	} else {                                                                          // given names
		for (size_t i = 0; i < idx.size(); ++i) {
			if (std::find(lociNames.begin(), lociNames.end(), names[i]) == lociNames.end())
				lociNames.insert(lociNames.begin() + idx[i] + i + 1, names[i]);
			else
				throw ValueError("Locus name " + names[i] + " already exists.");
		}
	}
	return *new GenoStructure(gs.m_ploidy, loci, gs.m_sexChrom, lociPos,
		gs.m_chromNames, gs.m_alleleNames, lociNames, gs.m_maxAllele, gs.m_infoFields);
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
	else if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 2)
		return "diploid";
	else if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 3)
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


string GenoStruTrait::alleleName(const Allele allele) const
{
#ifndef BINARYALLELE
	DBG_FAILIF(allele > s_genoStruRepository[m_genoStruIdx].m_maxAllele,
		IndexError, "Allele out of range of 0 ~ " +
		toStr(s_genoStruRepository[m_genoStruIdx].m_maxAllele));
	if (allele < s_genoStruRepository[m_genoStruIdx].m_alleleNames.size() ) {
		DBG_FAILIF(allele >= s_genoStruRepository[m_genoStruIdx].m_alleleNames.size(),
			IndexError, "No name for allele " + toStr(static_cast<UINT>(allele)));

		return s_genoStruRepository[m_genoStruIdx].m_alleleNames[allele];
	} else
		return toStr(static_cast<int>(allele));
#else
	if (static_cast<unsigned>(allele) < s_genoStruRepository[m_genoStruIdx].m_alleleNames.size() )
		return s_genoStruRepository[m_genoStruIdx].m_alleleNames[allele];
	else if (allele)
		return "1";
	else
		return "0";
#endif

}


UINT GenoStruTrait::infoIdx(const string & name) const
{
	vectorstr & names = s_genoStruRepository[m_genoStruIdx].m_infoFields;

	for (UINT i = 0; i < names.size(); ++i)
		if (names[i] == name)
			return i;
	throw IndexError("Info field '" + name + "' is not found. "
		                                     "Plese use infoFields=['" + name + "'] option of population() during construction\n"
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


/// CPPONLY should should only be called from population
GenoStructure & GenoStruTrait::struSetInfoFields(const vectorstr & fields)
{
	GenoStructure * gs = new GenoStructure(s_genoStruRepository[m_genoStruIdx]);

	gs->m_infoFields = fields;
	return *gs;
}


}
