/**
 *  $File: individual.cpp $
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
		DBG_DO(DBG_POPULATION, cerr << "Geno stru different" << endl);
		return false;
	}

	if (ISSETFLAG(m_flags, m_flagFemale) != ISSETFLAG(rhs.m_flags, m_flagFemale)
	    || ISSETFLAG(m_flags, m_flagAffected) != ISSETFLAG(rhs.m_flags, m_flagAffected) ) {
		DBG_DO(DBG_POPULATION, cerr << "Flags different: sex "
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
			DBG_DO(DBG_POPULATION, cerr << "Information field " << infoField(i) << " differ" << endl);
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


UINT individual::allele(UINT idx, int p, int chrom) const
{
	DBG_FAILIF(p < 0 && chrom >= 0, ValueError,
		"A valid ploidy index has to be specified if chrom is non-positive");
	if (p < 0) {
		CHECKRANGEGENOSIZE(idx);
		return static_cast<UINT>(*(m_genoPtr + idx));
	} else if (chrom < 0) {
		CHECKRANGEABSLOCUS(idx);
		CHECKRANGEPLOIDY(static_cast<UINT>(p));
		return static_cast<UINT>(*(m_genoPtr + idx + p * totNumLoci() ));
	} else {
		CHECKRANGELOCUS(chrom, idx);
		CHECKRANGEPLOIDY(static_cast<UINT>(p));
		CHECKRANGECHROM(static_cast<UINT>(chrom));
		return static_cast<UINT>(*(m_genoPtr + idx + p * totNumLoci() + chromBegin(chrom)));
	}
}


string individual::alleleChar(UINT idx, int p, int chrom) const
{
	DBG_FAILIF(p < 0 && chrom >= 0, ValueError,
		"A valid ploidy index has to be specified if chrom is non-positive");
	if (p < 0) {
		CHECKRANGEGENOSIZE(idx);
		return validIndex(idx) ? alleleName(allele(idx), idx % totNumLoci()) : "_";
	} else if (chrom < 0) {
		CHECKRANGEABSLOCUS(idx);
		CHECKRANGEPLOIDY(static_cast<UINT>(p));
		return validIndex(idx, p) ? alleleName(allele(idx, p), idx) : "_";
	} else {
		CHECKRANGELOCUS(static_cast<UINT>(chrom), idx);
		CHECKRANGEPLOIDY(static_cast<UINT>(p));
		CHECKRANGECHROM(static_cast<UINT>(chrom));

		return validIndex(idx, p, chrom) ? alleleName(allele(idx, p, chrom), idx + chromBegin(chrom)) : "_";
	}
}


PyObject * individual::genotype(int p, int chrom)
{
	DBG_FAILIF(p < 0 && chrom >= 0, ValueError,
		"A valid ploidy index has to be specified if chrom is non-positive");
	if (p < 0) {
		return Allele_Vec_As_NumArray(m_genoPtr, m_genoPtr + genoSize());
	} else if (chrom < 0) {
		CHECKRANGEPLOIDY(static_cast<UINT>(p));
		return Allele_Vec_As_NumArray(m_genoPtr + p * totNumLoci(),
			m_genoPtr + (p + 1) * totNumLoci() );
	} else {
		CHECKRANGEPLOIDY(static_cast<UINT>(p));
		CHECKRANGECHROM(static_cast<UINT>(chrom));
		return Allele_Vec_As_NumArray(m_genoPtr + p * totNumLoci() + chromBegin(chrom),
			m_genoPtr + p * totNumLoci() + chromEnd(chrom));
	}
}


void individual::setAllele(Allele allele, UINT idx, int p, int chrom)
{
	DBG_FAILIF(p < 0 && chrom >= 0, ValueError,
		"A valid ploidy index has to be specified if chrom is non-positive");
	if (p < 0) {
		CHECKRANGEGENOSIZE(idx);
		*(m_genoPtr + idx) = allele;
	} else if (chrom < 0) {
		CHECKRANGEABSLOCUS(idx);
		CHECKRANGEPLOIDY(static_cast<UINT>(p));
		*(m_genoPtr + idx + p * totNumLoci()) = allele;
	} else {
		CHECKRANGELOCUS(static_cast<UINT>(chrom), idx);
		CHECKRANGEPLOIDY(static_cast<UINT>(p));
		CHECKRANGECHROM(static_cast<UINT>(chrom));
		*(m_genoPtr + idx + p * totNumLoci() + chromBegin(chrom) ) = allele;
	}
}


void individual::setGenotype(const vectora & geno, int p, int chrom)
{
	DBG_FAILIF(p < 0 && chrom >= 0, ValueError,
		"A valid ploidy index has to be specified if chrom is non-positive");
	UINT sz = geno.size();
	if (p < 0) {
		for (UINT i = 0; i < totNumLoci() * ploidy(); i++)
			*(m_genoPtr + i) = geno[i % sz];
	} else if (chrom < 0) {
		CHECKRANGEPLOIDY(static_cast<UINT>(p));
		GenoIterator ptr = m_genoPtr + p * totNumLoci();

		UINT sz = geno.size();
		for (UINT i = 0; i < totNumLoci(); i++)
			*(ptr + i) = geno[i % sz];

	} else {
		CHECKRANGEPLOIDY(static_cast<UINT>(p));
		CHECKRANGECHROM(static_cast<UINT>(chrom));
		GenoIterator ptr = m_genoPtr + p * totNumLoci() + chromBegin(chrom);

		UINT sz = geno.size();
		for (UINT i = 0; i < numLoci(chrom); i++)
			*(ptr + i) = geno[i % sz];
	}
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

	if (ModuleMaxAllele == 1)
		width = 1;
	else
		width = 3;
	display(os, width, vectoru());
	const string & str = os.str();
	if (str.size() < 100)
		return str;
	else
		return str.substr(0, 100) + "...";
}


void individual::display(ostream & out, int width, const vectoru & loci)
{
	out << sexChar() << affectedChar() << " ";
	UINT pEnd = ploidy();
	for (UINT p = 0; p < pEnd;  ++p) {
		if (loci.empty()) {
			for (UINT ch = 0, chEnd = numChrom(); ch < chEnd; ++ch) {
				for (UINT j = 0, jEnd = numLoci(ch); j < jEnd;  ++j)
					out << setw(width) << alleleChar(j, p, ch);
				out << " ";
			}
		} else {
			for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc)
				out << setw(width) << alleleChar(*loc, p);
			out << " ";
		}

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
