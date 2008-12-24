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

#include "individual.h"
#include <sstream>
using std::ostringstream;
using std::setprecision;

namespace simuPOP {


individual & individual::operator=(const individual & rhs)
{
	m_flags = rhs.m_flags;
	setGenoPtr(rhs.genoPtr());
	setInfoPtr(rhs.infoPtr());
	// also copy genoStru pointer...
	this->setGenoStruIdx(rhs.genoStruIdx());
	return *this;
}


individual & individual::copyFrom(const individual & rhs)
{
	m_flags = rhs.m_flags;
	copy(rhs.genoBegin(), rhs.genoEnd(), genoBegin());
	copy(rhs.infoBegin(), rhs.infoEnd(), infoBegin());
	// also copy genoStru pointer...
	this->setGenoStruIdx(rhs.genoStruIdx());
	return *this;
}


bool individual::operator==(const individual & rhs) const
{
	if (genoStruIdx() != rhs.genoStruIdx() ) {
		DBG_DO(DBG_POPULATION, cout << "Geno stru different" << endl);
		return false;
	}

	if (ISSETFLAG(m_flags, m_flagFemale) != ISSETFLAG(rhs.m_flags, m_flagFemale)
	    || ISSETFLAG(m_flags, m_flagAffected) != ISSETFLAG(rhs.m_flags, m_flagAffected) ) {
		DBG_DO(DBG_POPULATION, cout << "Flags different: sex "
			                        << ISSETFLAG(m_flags, m_flagFemale) << " vs " << ISSETFLAG(rhs.m_flags, m_flagFemale) << ", aff "
			                        << ISSETFLAG(m_flags, m_flagAffected) << " vs " << ISSETFLAG(rhs.m_flags, m_flagAffected)
			                        << endl);
		return false;
	}

	for (UINT i = 0, iEnd = genoSize(); i < iEnd;  ++i)
		if (*(m_genoPtr + i) != *(rhs.m_genoPtr + i) )
			return false;

	for (UINT i = 0, iEnd = infoSize(); i < iEnd;  ++i)
		if (*(m_infoPtr + i) != *(rhs.m_infoPtr + i) ) {
			DBG_DO(DBG_POPULATION, cout << "Information field " << infoField(i) << " differ" << endl);
			return false;
		}
	return true;
}


int individual::__cmp__(const individual & rhs) const
{
	return (*this == rhs) ? 0 : 1;
}


bool individual::validIndex(UINT idx) const
{
	UINT cnt = totNumLoci();

	return validIndex(idx % cnt, idx / cnt);
}


bool individual::validIndex(UINT idx, UINT p) const
{
	std::pair<UINT, UINT> chIdx = chromLocusPair(idx);
	return validIndex(chIdx.second, p, chIdx.first);
}


bool individual::validIndex(UINT idx, UINT p, UINT ch) const
{
	// well, this might change later.
	if (ploidy() != 2)
		return true;

	if (p == 1 && isHaplodiploid() && sex() == Male)
		return false;

	Sex s = sex();
	UINT t = chromType(ch);
	if ((s == Female && t == ChromosomeY) ||    // female chromsome Y
	    (s == Male &&                           // second copy of chromosome X and first copy of chromosome Y
	     ((p == 1 && t == ChromosomeX) || (p == 0 && t == ChromosomeY))))
		return false;

	return true;
}


string individual::alleleChar(UINT idx) const
{
	CHECKRANGEGENOSIZE(idx);

	return validIndex(idx) ? alleleName(allele(idx)) : "_";
}


string individual::alleleChar(UINT idx, UINT p) const
{
	CHECKRANGEABSLOCUS(idx);
	CHECKRANGEPLOIDY(p);

	return validIndex(idx, p) ? alleleName(allele(idx, p)) : "_";
}


string individual::alleleChar(UINT idx, UINT p, UINT ch) const
{
	CHECKRANGELOCUS(ch, idx);
	CHECKRANGEPLOIDY(p);
	CHECKRANGECHROM(ch);

	return validIndex(idx, p, ch) ? alleleName(allele(idx, p, ch)) : "_";
}


void individual::swap(individual & ind, bool swapContent)
{
	if (genoStruIdx() != ind.genoStruIdx() )
		throw SystemError("Can only swap individuals with different geno structure.");

	std::swap(m_infoPtr, ind.m_infoPtr);

	if (swapContent) {
		Allele tmp;
		for (UINT i = 0, iEnd = genoSize(); i < iEnd;  i++) {
			tmp = m_genoPtr[i];
			m_genoPtr[i] = ind.m_genoPtr[i];
			ind.m_genoPtr[i] = tmp;
		}
	} else {
		std::swap(m_genoPtr, ind.m_genoPtr);
	}
}


string individual::__repr__()
{
	ostringstream os;
	int width = 1;

	if (MaxAllele() == 1)
		width = 1;
	else
		width = 3;
	display(os, width, vectori(), vectori());
	const string & str = os.str();
	if (str.size() < 100)
		return str;
	else
		return str.substr(0, 100) + "...";
}


void individual::display(ostream & out, int width, const vectori & chrom, const vectori & loci)
{
	out << sexChar() << affectedChar() << " ";
	UINT pEnd = ploidy();
	for (UINT p = 0; p < pEnd;  ++p) {
		if (chrom.empty() && loci.empty()) {
			for (UINT ch = 0, chEnd = numChrom(); ch < chEnd; ++ch) {
				for (UINT j = 0, jEnd = numLoci(ch); j < jEnd;  ++j)
					out << setw(width) << alleleChar(j, p, ch);
				out << " ";
			}
		} else if (!chrom.empty() && loci.empty()) {
			for (vectori::const_iterator ch = chrom.begin(); ch != chrom.end(); ++ch) {
				for (UINT j = 0, jEnd = numLoci(*ch); j < jEnd;  ++j)
					out << setw(width) << alleleChar(j, p, *ch);
				out << " ";
			}
		} else if (chrom.empty() && !loci.empty()) {
			for (vectori::const_iterator loc = loci.begin(); loc != loci.end(); ++loc)
				out << setw(width) << alleleChar(*loc, p);
			out << " ";
		} else  // both specified
			throw ValueError("Please specify only one of chrom and loci.");

		if (p != pEnd - 1)
			out << "| ";
	}
	if (infoSize() != 0) {
		out << "| ";
		for (vectorinfo::const_iterator info = infoBegin(); info != infoEnd(); ++info)
			out << " " << setprecision(2) << *info;
	}
}


individual & pyIndIterator::next()
{
	// this is the easy (and faster) case
	if (m_allInds) {
		if (m_index == m_end)
			throw StopIteration("");
		else
			return *m_index++;
	}
	// check the visibility of individuals
	do {
		if (m_index == m_end)
			throw StopIteration("");
		else if (m_index->visible()) {
			if (m_allVisibles || m_index->iteratable())
				return *m_index++;
			else
				++m_index;
		} else
			++m_index;
	} while (true);
}


}
