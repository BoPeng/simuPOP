/**
 *  $File: Individual.cpp $
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

#include "individual.h"
#include <sstream>
using std::ostringstream;
using std::setprecision;

#if PY_VERSION_HEX >= 0x03000000
#  define PyInt_FromLong(x) PyLong_FromLong(x)
#endif

namespace simuPOP {


Individual & Individual::operator=(const Individual & rhs)
{
	m_flags = rhs.m_flags;
	setGenoPtr(rhs.genoPtr());
	setInfoPtr(rhs.infoPtr());
	// also copy genoStru pointer...
	this->setGenoStruIdx(rhs.genoStruIdx());
	return *this;
}


Individual & Individual::copyFrom(const Individual & rhs)
{
	m_flags = rhs.m_flags;
	copy(rhs.genoBegin(), rhs.genoEnd(), genoBegin());
	copy(rhs.infoBegin(), rhs.infoEnd(), infoBegin());
	// also copy genoStru pointer...
	this->setGenoStruIdx(rhs.genoStruIdx());
	return *this;
}


bool Individual::operator==(const Individual & rhs) const
{
	if (genoStruIdx() != rhs.genoStruIdx()) {
		DBG_DO(DBG_POPULATION, cerr << "Geno stru different" << endl);
		return false;
	}

	if (ISSETFLAG(m_flags, m_flagFemale) != ISSETFLAG(rhs.m_flags, m_flagFemale)
	    || ISSETFLAG(m_flags, m_flagAffected) != ISSETFLAG(rhs.m_flags, m_flagAffected)) {
		DBG_DO(DBG_POPULATION, cerr << "Flags different: sex "
			                        << ISSETFLAG(m_flags, m_flagFemale) << " vs " << ISSETFLAG(rhs.m_flags, m_flagFemale) << ", aff "
			                        << ISSETFLAG(m_flags, m_flagAffected) << " vs " << ISSETFLAG(rhs.m_flags, m_flagAffected)
			                        << endl);
		return false;
	}

	for (UINT i = 0, iEnd = genoSize(); i < iEnd; ++i)
		if (*(m_genoPtr + i) != *(rhs.m_genoPtr + i))
			return false;

	for (UINT i = 0, iEnd = infoSize(); i < iEnd; ++i)
		if (*(m_infoPtr + i) != *(rhs.m_infoPtr + i)) {
			DBG_DO(DBG_POPULATION, cerr << "Information field " << infoField(i) << " differ" << endl);
			return false;
		}
	return true;
}


int Individual::__cmp__(const Individual & rhs) const
{
	return (*this == rhs) ? 0 : 1;
}


bool Individual::validIndex(UINT idx) const
{
	UINT cnt = totNumLoci();

	return validIndex(idx % cnt, idx / cnt);
}


bool Individual::validIndex(UINT idx, UINT p) const
{
	pairu chIdx = chromLocusPair(idx);

	return validIndex(chIdx.second, p, chIdx.first);
}


bool Individual::validIndex(UINT idx, UINT p, UINT ch) const
{
	// well, this might change later.
	if (ploidy() != 2)
		return true;

	if (p == 1 && isHaplodiploid() && sex() == MALE)
		return false;

	Sex s = sex();
	UINT t = chromType(ch);
	if ((s == FEMALE && t == CHROMOSOME_Y) ||       // female chromsome Y
	    (s == MALE &&                               // second copy of chromosome X and first copy of chromosome Y
	     ((p == 1 && t == CHROMOSOME_X) || (p == 0 && t == CHROMOSOME_Y))))
		return false;

	return true;
}


UINT Individual::allele(UINT idx, int p, int chrom) const
{
	DBG_FAILIF(p < 0 && chrom >= 0, ValueError,
		"A valid ploidy index has to be specified if chrom is non-positive");
	if (p < 0) {
		CHECKRANGEGENOSIZE(idx);
		return static_cast<UINT>(*(m_genoPtr + idx));
	} else if (chrom < 0) {
		CHECKRANGEABSLOCUS(idx);
		CHECKRANGEPLOIDY(static_cast<UINT>(p));
		return static_cast<UINT>(*(m_genoPtr + idx + p * totNumLoci()));
	} else {
		CHECKRANGELOCUS(chrom, idx);
		CHECKRANGEPLOIDY(static_cast<UINT>(p));
		CHECKRANGECHROM(static_cast<UINT>(chrom));
		return static_cast<UINT>(*(m_genoPtr + idx + p * totNumLoci() + chromBegin(chrom)));
	}
}


string Individual::alleleChar(UINT idx, int p, int chrom) const
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


PyObject * Individual::genotype(const uintList & ply, const uintList & ch)
{
	DBG_WARNIF(true, "The returned object of function Individual.genotype() is a special "
		             "carray object that reflects the underlying genotype of an individual. "
		             "It will become invalid once the population changes. Please use "
		             "list(ind.genotype()) if you would like to keep a copy of genotypes");

	size_t beginP = 0;
	size_t endP = 0;
	size_t beginCh = 0;
	size_t endCh = 0;

	if (ply.allAvail())
		endP = ploidy();
	else {
		const vectoru & ploidys = ply.elems();
		if (ploidys.empty()) {
			Py_INCREF(Py_None);
			return Py_None;
		}
		beginP = ploidys[0];
		endP = ploidys[0];
		CHECKRANGEPLOIDY(static_cast<UINT>(beginP));
		for (size_t i = 1; i < ploidys.size(); ++i) {
			CHECKRANGEPLOIDY(static_cast<UINT>(ploidys[i]));
			if (beginP > ploidys[i])
				beginP = ploidys[i];
			if (endP < ploidys[i])
				endP = ploidys[i];
		}
		++endP;
	}
	if (ch.allAvail())
		endCh = numChrom();
	else {
		const vectoru & chroms = ch.elems();
		if (chroms.empty()) {
			Py_INCREF(Py_None);
			return Py_None;
		}
		beginCh = chroms[0];
		endCh = chroms[0];
		CHECKRANGECHROM(static_cast<UINT>(beginCh));
		for (size_t i = 1; i < chroms.size(); ++i) {
			CHECKRANGECHROM(static_cast<UINT>(chroms[i]));
			if (beginCh > chroms[i])
				beginCh = chroms[i];
			if (endCh < chroms[i])
				endCh = chroms[i];
		}
		++endCh;
	}

	if (endP > beginP + 1) {
		// has to be all chromosomes
		DBG_FAILIF(beginCh != 0 || endCh != numChrom(), ValueError,
			"If multiple ploidy are chosen, all chromosomes has to be chosen.");
		return Allele_Vec_As_NumArray(m_genoPtr + beginP * totNumLoci(),
			m_genoPtr + endP * totNumLoci());
	} else
		return Allele_Vec_As_NumArray(m_genoPtr + beginP * totNumLoci() + chromBegin(beginCh),
			m_genoPtr + beginP * totNumLoci() + chromEnd(endCh - 1));
}


PyObject * Individual::genoAtLoci(const lociList & lociList)
{
	size_t ply = ploidy();

	if (isHaplodiploid() && sex() == MALE)
		ply = 1;

	vectori alleles;

	if (lociList.allAvail()) {

		alleles.reserve(ply * totNumLoci());

		for (size_t ch = 0; ch < numChrom(); ++ch) {
			UINT chType = chromType(ch);
			if (chType == CHROMOSOME_Y && sex() == FEMALE)
				continue;
			for (size_t idx = 0; idx < numLoci(ch); ++idx) {
				for (size_t p = 0; p < ply; ++p) {
					if (((chType == CHROMOSOME_X && p == 1) ||
					     (chType == CHROMOSOME_Y && p == 0)) && sex() == MALE)
						continue;
					alleles.push_back(allele(idx, p, ch));
				}
			}
		}
	} else {
		const vectoru & loci = lociList.elems(this);

		vectoru chromTypes;

		for (size_t j = 0; j < loci.size(); ++j)
			chromTypes.push_back(chromType(chromLocusPair(loci[j]).first));

		alleles.reserve(ply * loci.size());

		for (size_t idx = 0; idx < loci.size(); ++idx) {
			for (size_t p = 0; p < ply; ++p) {
				if (chromTypes[idx] == CHROMOSOME_Y && sex() == FEMALE)
					continue;
				if (((chromTypes[idx] == CHROMOSOME_X && p == 1) ||
				     (chromTypes[idx] == CHROMOSOME_Y && p == 0)) && sex() == MALE)
					continue;
				alleles.push_back(allele(loci[idx], p));
			}
		}
	}
	PyObject * genoObj = PyTuple_New(alleles.size());
	// set value
	for (size_t j = 0; j < alleles.size(); ++j)
		PyTuple_SET_ITEM(genoObj, j, PyInt_FromLong(alleles[j]));
	return genoObj;
}


void Individual::setAllele(Allele allele, UINT idx, int p, int chrom)
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
		*(m_genoPtr + idx + p * totNumLoci() + chromBegin(chrom)) = allele;
	}
}


void Individual::setGenotype(const uintList & genoList, const uintList & ply, const uintList & ch)
{
	const vectoru & geno = genoList.elems();

	UINT sz = geno.size();
	size_t idx = 0;

	vectoru ploidys = ply.elems();

	if (ply.allAvail()) {
		for (size_t i = 0; i < ploidy(); ++i)
			ploidys.push_back(i);
	} else {
#ifndef OPTIMIZED
		for (size_t i = 0; i < ploidys.size(); ++i) {
			CHECKRANGEPLOIDY(static_cast<UINT>(ploidys[i]));
		}
#endif
	}
	vectoru chroms = ch.elems();
	if (ch.allAvail()) {
		for (size_t i = 0; i < numChrom(); ++i)
			chroms.push_back(i);
	} else {
#ifndef OPTIMIZED
		for (size_t i = 0; i < chroms.size(); ++i) {
			CHECKRANGECHROM(static_cast<UINT>(chroms[i]));
		}
#endif
	}
	for (size_t i = 0; i < ploidys.size(); ++i) {
		size_t p = ploidys[i];
		for (size_t j = 0; j < chroms.size(); ++j) {
			size_t chrom = chroms[j];
			GenoIterator ptr = m_genoPtr + p * totNumLoci() + chromBegin(chrom);

			for (UINT i = 0; i < numLoci(chrom); i++, ++idx)
				*(ptr + i) = ToAllele(geno[idx % sz]);
		}
	}
}


void Individual::swap(Individual & ind, bool swapContent)
{
	if (genoStruIdx() != ind.genoStruIdx())
		throw SystemError("Can only swap individuals with different geno structure.");

	std::swap(m_infoPtr, ind.m_infoPtr);

	if (swapContent) {
		Allele tmp;
		for (UINT i = 0, iEnd = genoSize(); i < iEnd; i++) {
			tmp = m_genoPtr[i];
			m_genoPtr[i] = ind.m_genoPtr[i];
			ind.m_genoPtr[i] = tmp;
		}
	} else {
		std::swap(m_genoPtr, ind.m_genoPtr);
	}
}


void Individual::display(ostream & out, int width, const vectoru & loci)
{
	out << (sex() == MALE ? 'M' : 'F') << (affected() ? 'A' : 'U') << " ";
	UINT pEnd = ploidy();
	for (UINT p = 0; p < pEnd; ++p) {
		if (loci.empty()) {
			for (UINT ch = 0, chEnd = numChrom(); ch < chEnd; ++ch) {
				for (UINT j = 0, jEnd = numLoci(ch); j < jEnd; ++j)
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
		for (vectorf::const_iterator info = infoBegin(); info != infoEnd(); ++info)
			out << " " << *info;
	}
}


}
