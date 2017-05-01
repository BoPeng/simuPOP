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
#include <set>
#include <map>
#include <algorithm>
using std::ostringstream;
using std::setprecision;

#if PY_VERSION_HEX >= 0x03000000
#  define PyInt_FromLong(x) PyLong_FromLong(x)
#endif
extern "C" PyObject * PyDefDict_New();

namespace simuPOP {


Individual & Individual::operator=(const Individual & rhs)
{
	m_flags = rhs.m_flags;
	setGenoPtr(rhs.genoPtr());
	setInfoPtr(rhs.infoPtr());
	LINEAGE_EXPR(setLineagePtr(rhs.lineagePtr()));
	// also copy genoStru pointer...
	this->setGenoStruIdx(rhs.genoStruIdx());
	return *this;
}


Individual & Individual::copyFrom(const Individual & rhs)
{
	m_flags = rhs.m_flags;
#ifdef MUTANTALLELE
	copyGenotype(rhs.genoBegin(), rhs.genoEnd(), genoBegin());
#else
	copy(rhs.genoBegin(), rhs.genoEnd(), genoBegin());
#endif
	copy(rhs.infoBegin(), rhs.infoEnd(), infoBegin());
	LINEAGE_EXPR(copy(rhs.lineageBegin(), rhs.lineageEnd(), lineageBegin()));
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

#ifndef MUTANTALLELE
	for (size_t i = 0, iEnd = genoSize(); i < iEnd; ++i)
		if (*(m_genoPtr + i) != *(rhs.m_genoPtr + i))
			return false;
#else
	vectorm::const_val_iterator it = m_genoPtr.get_val_iterator();
	vectorm::const_val_iterator it_end = (m_genoPtr + genoSize()).get_val_iterator();
	size_t s_idx = m_genoPtr.index();
	vectorm::const_val_iterator rit = rhs.m_genoPtr.get_val_iterator();
	vectorm::const_val_iterator rit_end = (rhs.m_genoPtr + genoSize()).get_val_iterator();
	size_t rs_idx = rhs.m_genoPtr.index();
	for (; it != it_end; ++it, ++rit)
		if (rit == rit_end || ((it->first - s_idx) != (rit->first - rs_idx)) || it->second != rit->second)
			return false;
	if (rit != rit_end)
		return false;
#endif

#ifdef LINEAGE
	for (size_t i = 0, iEnd = genoSize(); i < iEnd; ++i)
		if (*(m_lineagePtr + i) != *(rhs.m_lineagePtr + i))
			return false;
#endif

	for (size_t i = 0, iEnd = infoSize(); i < iEnd; ++i)
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


bool Individual::validIndex(size_t idx) const
{
	size_t cnt = totNumLoci();

	return validIndex(idx % cnt, idx / cnt);
}


bool Individual::validIndex(size_t idx, size_t p) const
{
	pairu chIdx = chromLocusPair(idx);

	return validIndex(chIdx.second, p, chIdx.first);
}


bool Individual::validIndex(size_t /* idx */, size_t p, size_t ch) const
{
	// well, this might change later.
	if (ploidy() != 2)
		return true;

	if (p == 1 && isHaplodiploid() && sex() == MALE)
		return false;

	Sex s = sex();
	size_t t = chromType(ch);
	if ((s == FEMALE && t == CHROMOSOME_Y) ||       // female chromsome Y
	    (s == MALE &&                               // second copy of chromosome X and first copy of chromosome Y
	     ((p == 1 && t == CHROMOSOME_X) || (p == 0 && t == CHROMOSOME_Y))) ||
	    (p > 0 && t == MITOCHONDRIAL))              // mitochondiral DNA only uses the first homologous copy
		return false;

	return true;
}


ULONG Individual::allele(size_t idx, ssize_t p, ssize_t chrom) const
{
	DBG_FAILIF(p < 0 && chrom >= 0, ValueError,
		"A valid ploidy index has to be specified if chrom is non-positive");
	if (p < 0) {
		CHECKRANGEGENOSIZE(idx);
		return static_cast<ULONG>(DEREF_ALLELE(m_genoPtr + idx));
	} else if (chrom < 0) {
		CHECKRANGEABSLOCUS(idx);
		CHECKRANGEPLOIDY(static_cast<size_t>(p));
		return static_cast<ULONG>(DEREF_ALLELE(m_genoPtr + idx + p * totNumLoci()));
	} else {
		CHECKRANGELOCUS(chrom, idx);
		CHECKRANGEPLOIDY(static_cast<size_t>(p));
		CHECKRANGECHROM(static_cast<size_t>(chrom));
		return static_cast<ULONG>(DEREF_ALLELE(m_genoPtr + idx + p * totNumLoci() + chromBegin(chrom)));
	}
}


long Individual::alleleLineage(size_t idx, ssize_t p, ssize_t chrom) const
{
#ifdef LINEAGE
	DBG_FAILIF(p < 0 && chrom >= 0, ValueError,
		"A valid ploidy index has to be specified if chrom is non-positive");
	if (p < 0) {
		CHECKRANGEGENOSIZE(idx);
		return *(m_lineagePtr + idx);
	} else if (chrom < 0) {
		CHECKRANGEABSLOCUS(idx);
		CHECKRANGEPLOIDY(static_cast<size_t>(p));
		return *(m_lineagePtr + idx + p * totNumLoci());
	} else {
		CHECKRANGELOCUS(chrom, idx);
		CHECKRANGEPLOIDY(static_cast<size_t>(p));
		CHECKRANGECHROM(static_cast<size_t>(chrom));
		return *(m_lineagePtr + idx + p * totNumLoci() + chromBegin(chrom));
	}
#else
	(void)idx;
	(void)p;
	(void)chrom;
	return 0;
#endif
}


string Individual::alleleChar(size_t idx, ssize_t p, ssize_t chrom) const
{
	DBG_FAILIF(p < 0 && chrom >= 0, ValueError,
		"A valid ploidy index has to be specified if chrom is non-positive");
	if (p < 0) {
		CHECKRANGEGENOSIZE(idx);
		return validIndex(idx) ? alleleName(allele(idx), idx % totNumLoci()) : "_";
	} else if (chrom < 0) {
		CHECKRANGEABSLOCUS(idx);
		CHECKRANGEPLOIDY(static_cast<size_t>(p));
		return validIndex(idx, p) ? alleleName(allele(idx, p), idx) : "_";
	} else {
		CHECKRANGELOCUS(static_cast<size_t>(chrom), idx);
		CHECKRANGEPLOIDY(static_cast<size_t>(p));
		CHECKRANGECHROM(static_cast<size_t>(chrom));

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
		CHECKRANGEPLOIDY(static_cast<size_t>(beginP));
		for (size_t i = 1; i < ploidys.size(); ++i) {
			CHECKRANGEPLOIDY(static_cast<size_t>(ploidys[i]));
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
		CHECKRANGECHROM(static_cast<size_t>(beginCh));
		for (size_t i = 1; i < chroms.size(); ++i) {
			CHECKRANGECHROM(static_cast<size_t>(chroms[i]));
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


pyMutantIterator Individual::mutants(const uintList & ply, const uintList & ch)
{
	size_t beginP = 0;
	size_t endP = 0;
	size_t beginCh = 0;
	size_t endCh = 0;

	if (ply.allAvail())
		endP = ploidy();
	else {
		const vectoru & ploidys = ply.elems();
		if (ploidys.empty())
			return pyMutantIterator(m_genoPtr, 0, 0, 1);
		beginP = ploidys[0];
		endP = ploidys[0];
		CHECKRANGEPLOIDY(static_cast<size_t>(beginP));
		for (size_t i = 1; i < ploidys.size(); ++i) {
			CHECKRANGEPLOIDY(static_cast<size_t>(ploidys[i]));
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
		if (chroms.empty())
			return pyMutantIterator(m_genoPtr, 0, 0, 1);
		beginCh = chroms[0];
		endCh = chroms[0];
		CHECKRANGECHROM(static_cast<size_t>(beginCh));
		for (size_t i = 1; i < chroms.size(); ++i) {
			CHECKRANGECHROM(static_cast<size_t>(chroms[i]));
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
		return pyMutantIterator(m_genoPtr, beginP * totNumLoci(),
			endP * totNumLoci(), genoSize());
	} else
		return pyMutantIterator(m_genoPtr, beginP * totNumLoci() + chromBegin(beginCh),
			beginP * totNumLoci() + chromEnd(endCh - 1), genoSize());
}


PyObject * Individual::lineage(const uintList & ply, const uintList & ch)
{
#ifdef LINEAGE
	DBG_WARNIF(true, "The returned object of function Individual.lineage() is a special "
		             "carray_lineage object that reflects the underlying genotype lineage of an individual. "
		             "It will become invalid once the population changes. Please use "
		             "list(ind.lineage()) if you would like to keep a copy of lineages");

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
		CHECKRANGEPLOIDY(static_cast<size_t>(beginP));
		for (size_t i = 1; i < ploidys.size(); ++i) {
			CHECKRANGEPLOIDY(static_cast<size_t>(ploidys[i]));
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
		CHECKRANGECHROM(static_cast<size_t>(beginCh));
		for (size_t i = 1; i < chroms.size(); ++i) {
			CHECKRANGECHROM(static_cast<size_t>(chroms[i]));
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
		return Lineage_Vec_As_NumArray(m_lineagePtr + beginP * totNumLoci(),
			m_lineagePtr + endP * totNumLoci());
	} else
		return Lineage_Vec_As_NumArray(m_lineagePtr + beginP * totNumLoci() + chromBegin(beginCh),
			m_lineagePtr + beginP * totNumLoci() + chromEnd(endCh - 1));
#else
	(void)ply;
	(void)ch;
	Py_INCREF(Py_None);
	return Py_None;
#endif
}


PyObject * Individual::mutAtLoci(const lociList & lociList)
{
	ssize_t ply = ploidy();

	PyObject * mutDict = PyDefDict_New();
	PyObject * dkey = NULL;
	PyObject * dval = NULL;

	if (isHaplodiploid() && sex() == MALE)
		ply = 1;

	if (lociList.allAvail()) {
#ifdef MUTANTALLELE
		bool autosome_only = chromX() == -1 && chromY() == -1;
		if (autosome_only) {

			vectorm::const_val_iterator m_ptr = m_genoPtr.get_val_iterator();
			vectorm::const_val_iterator m_end = (m_genoPtr + genoSize()).get_val_iterator();
			for (; m_ptr != m_end; ++m_ptr) {
				PyDict_SetItem(mutDict,
					dkey = PyInt_FromLong(m_ptr->first % genoSize()),
					dval = PyInt_FromLong(m_ptr->second));
				Py_DECREF(dkey);
				Py_DECREF(dval);
			}
		} else {
			vector<ULONG> tmp_mutants(ply * totNumLoci(), 0);

			vectorm::const_val_iterator m_ptr = m_genoPtr.get_val_iterator();
			vectorm::const_val_iterator m_end = (m_genoPtr + genoSize()).get_val_iterator();
			for (; m_ptr != m_end; ++m_ptr)
				tmp_mutants[ m_ptr->first % genoSize() ] = m_ptr->second;

			for (ssize_t p = 0; p < ply; ++p) {
				for (ssize_t ch = 0; ch < static_cast<ssize_t>(numChrom()); ++ch) {
					size_t chType = chromType(ch);
					if (chType == CHROMOSOME_Y && sex() == FEMALE)
						continue;
					if (((chType == CHROMOSOME_X && p == 1) ||
					     (chType == CHROMOSOME_Y && p == 0)) && sex() == MALE)
						continue;
					if (chType == MITOCHONDRIAL && p > 0)
						continue;
					size_t k = p * totNumLoci() + chromBegin(ch);
					for (size_t idx = 0; idx < numLoci(ch); ++idx)
						if (tmp_mutants[k + idx] != 0) {
							PyDict_SetItem(mutDict,
								dkey = PyInt_FromLong(k + idx),
								dval = PyInt_FromLong(tmp_mutants[k + idx]));
							Py_DECREF(dkey);
							Py_DECREF(dval);
						}
				}
			}
		}
#else
		for (ssize_t p = 0; p < ply; ++p) {
			for (ssize_t ch = 0; ch < static_cast<ssize_t>(numChrom()); ++ch) {
				size_t chType = chromType(ch);
				if (chType == CHROMOSOME_Y && sex() == FEMALE)
					continue;
				if (((chType == CHROMOSOME_X && p == 1) ||
				     (chType == CHROMOSOME_Y && p == 0)) && sex() == MALE)
					continue;
				if (chType == MITOCHONDRIAL && p > 0)
					continue;
				size_t k = p * totNumLoci() + chromBegin(ch);
				for (size_t idx = 0; idx < numLoci(ch); ++idx) {
					if (*(m_genoPtr + k + idx) != 0) {
						PyDict_SetItem(mutDict,
							dkey = PyInt_FromLong(k + idx),
							dval = PyInt_FromLong(*(m_genoPtr + k + idx)));
						Py_DECREF(dkey);
						Py_DECREF(dval);
					}
				}
			}
		}
#endif
	} else {
#ifdef MUTANTALLELE
		bool autosome_only = chromX() == -1 && chromY() == -1;
		if (autosome_only) {
			const vectoru & loci = lociList.elems(this);

			size_t nLoci = loci.size();

			// if there is a small number of loci, or if a small fraction of loci
			// get the loci directly.
			if (nLoci * 100 < totNumLoci()) {
				for (int p = 0; p < ply; ++p) {
					size_t k = p * totNumLoci();
					for (size_t idx = 0; idx < loci.size(); ++idx) {
						ULONG a = allele(loci[idx], p);
						if (a != 0) {
							PyDict_SetItem(mutDict,
								dkey = PyInt_FromLong(k + loci[idx]),
								dval = PyInt_FromLong(a));
							Py_DECREF(dkey);
							Py_DECREF(dval);
						}
					}
				}
			} else {
				// check if the loci is ordered
				bool ordered = true;
				for (size_t i = 1; i < loci.size(); ++i)
					if (loci[i] < loci[i - 1]) {
						ordered = false;
						break;
					}
				const vectoru * sorted = &loci;
				if (!ordered) {
					vectoru * sloci = new vectoru(loci.begin(), loci.end());
					std::sort(sloci->begin(), sloci->end());
					sorted = sloci;
				}

				size_t nLoci = totNumLoci();
				size_t l = sorted->front();
				size_t h = sorted->back();

				vectoru::const_iterator eit = sorted->end();
				for (size_t p = 0; p < ploidy(); ++p) {
					vectoru::const_iterator it = sorted->begin();
					vectorm::const_val_iterator m_ptr = (m_genoPtr + (l + p * nLoci)).get_val_iterator();
					vectorm::const_val_iterator m_end = (m_genoPtr + (h + 1 + p * nLoci)).get_val_iterator();
					for (; m_ptr != m_end; ++m_ptr) {
						size_t loc = m_ptr->first % nLoci;
						// move it
						while (*it < loc && it != eit)
							++it;
						// do not look any further
						if (it == eit)
							break;
						else if (*it == loc) {
							PyDict_SetItem(mutDict,
								dkey = PyInt_FromLong(loc + p * nLoci),
								dval = PyInt_FromLong(m_ptr->second));
							Py_DECREF(dkey);
							Py_DECREF(dval);
						}
					}
				}
				if (!ordered)
					delete sorted;
			}
		} else {
			const vectoru & loci = lociList.elems(this);

			size_t nLoci = loci.size();

			vectoru chromTypes;

			for (size_t j = 0; j < loci.size(); ++j)
				chromTypes.push_back(chromType(chromLocusPair(loci[j]).first));

			// if there is a small number of loci, or if a small fraction of loci
			// get the loci directly.
			if (nLoci * 100 < totNumLoci()) {
				for (int p = 0; p < ply; ++p) {
					size_t k = p * totNumLoci();
					for (size_t idx = 0; idx < loci.size(); ++idx) {
						if (chromTypes[idx] == CHROMOSOME_Y && sex() == FEMALE)
							continue;
						if (((chromTypes[idx] == CHROMOSOME_X && p == 1) ||
						     (chromTypes[idx] == CHROMOSOME_Y && p == 0)) && sex() == MALE)
							continue;
						if (chromTypes[idx] == MITOCHONDRIAL && p > 0)
							continue;
						ULONG a = allele(loci[idx], p);
						if (a) {
							PyDict_SetItem(mutDict,
								dkey = PyInt_FromLong(k + loci[idx]),
								dval = PyInt_FromLong(a));
							Py_DECREF(dkey);
							Py_DECREF(dval);
						}
					}
				}
			} else {
				// otherwise put all mutants in a temporary array
				vector<ULONG> tmp_mutants(ply * nLoci, 0);
				std::map<size_t, size_t> loci_map;
				for (size_t i = 0; i < nLoci; ++i)
					loci_map[loci[i]] = i;

				vectorm::const_val_iterator m_ptr = m_genoPtr.get_val_iterator();
				vectorm::const_val_iterator m_end = (m_genoPtr + genoSize()).get_val_iterator();
				for (; m_ptr != m_end; ++m_ptr) {
					size_t loc = m_ptr->first % genoSize();
					size_t p = loc / totNumLoci();
					// find it
					std::map<size_t, size_t>::iterator it = loci_map.find(loc % totNumLoci());
					if (it != loci_map.end())
						tmp_mutants[it->second + p * nLoci] = m_ptr->second;
				}

				for (int p = 0; p < ply; ++p) {
					for (size_t idx = 0; idx < loci.size(); ++idx) {
						if (chromTypes[idx] == CHROMOSOME_Y && sex() == FEMALE)
							continue;
						if (((chromTypes[idx] == CHROMOSOME_X && p == 1) ||
						     (chromTypes[idx] == CHROMOSOME_Y && p == 0)) && sex() == MALE)
							continue;
						if (chromTypes[idx] == MITOCHONDRIAL && p > 0)
							continue;
						if (tmp_mutants[idx + p * nLoci] != 0) {
							PyDict_SetItem(mutDict,
								dkey = PyInt_FromLong(loci[idx] + p * totNumLoci()),
								dval = PyInt_FromLong(tmp_mutants[idx + p * nLoci]));
							Py_DECREF(dkey);
							Py_DECREF(dval);
						}
					}
				}
			}
		}
#else
		const vectoru & loci = lociList.elems(this);

		vectoru chromTypes;

		for (size_t j = 0; j < loci.size(); ++j)
			chromTypes.push_back(chromType(chromLocusPair(loci[j]).first));


		for (int p = 0; p < ply; ++p) {
			for (size_t idx = 0; idx < loci.size(); ++idx) {
				if (chromTypes[idx] == CHROMOSOME_Y && sex() == FEMALE)
					continue;
				if (((chromTypes[idx] == CHROMOSOME_X && p == 1) ||
				     (chromTypes[idx] == CHROMOSOME_Y && p == 0)) && sex() == MALE)
					continue;
				if (chromTypes[idx] == MITOCHONDRIAL && p > 0)
					continue;
				if (allele(loci[idx], p) != 0) {
					PyDict_SetItem(mutDict,
						dkey = PyInt_FromLong(loci[idx] + p * totNumLoci()),
						dval = PyInt_FromLong(allele(loci[idx], p)));
					Py_DECREF(dkey);
					Py_DECREF(dval);
				}
			}
		}
#endif
	}
	return mutDict;
}


PyObject * Individual::genoAtLoci(const lociList & lociList)
{
	ssize_t ply = ploidy();

	if (isHaplodiploid() && sex() == MALE)
		ply = 1;

	vector<ULONG> alleles;

	if (lociList.allAvail()) {
#ifdef MUTANTALLELE
		alleles.reserve(ply * totNumLoci());
		vector<ULONG> tmp_alleles(ply * totNumLoci(), 0);

		vectorm::const_val_iterator m_ptr = m_genoPtr.get_val_iterator();
		vectorm::const_val_iterator m_end = (m_genoPtr + genoSize()).get_val_iterator();
		for (; m_ptr != m_end; ++m_ptr)
			tmp_alleles[ m_ptr->first % genoSize() ] = m_ptr->second;

		for (ssize_t ch = 0; ch < static_cast<ssize_t>(numChrom()); ++ch) {
			size_t chType = chromType(ch);
			if (chType == CHROMOSOME_Y && sex() == FEMALE)
				continue;
			for (size_t idx = 0; idx < numLoci(ch); ++idx) {
				for (ssize_t p = 0; p < ply; ++p) {
					if (((chType == CHROMOSOME_X && p == 1) ||
					     (chType == CHROMOSOME_Y && p == 0)) && sex() == MALE)
						continue;
					if (chType == MITOCHONDRIAL && p > 0)
						continue;
					alleles.push_back(tmp_alleles[p * totNumLoci() + chromBegin(ch) + idx]);
				}
			}
		}
#else
		alleles.reserve(ply * totNumLoci());

		for (ssize_t ch = 0; ch < static_cast<ssize_t>(numChrom()); ++ch) {
			size_t chType = chromType(ch);
			if (chType == CHROMOSOME_Y && sex() == FEMALE)
				continue;
			for (size_t idx = 0; idx < numLoci(ch); ++idx) {
				for (ssize_t p = 0; p < ply; ++p) {
					if (((chType == CHROMOSOME_X && p == 1) ||
					     (chType == CHROMOSOME_Y && p == 0)) && sex() == MALE)
						continue;
					if (chType == MITOCHONDRIAL && p > 0)
						continue;
					alleles.push_back(allele(idx, p, ch));
				}
			}
		}

#endif
	} else {
#ifdef MUTANTALLELE
		const vectoru & loci = lociList.elems(this);

		size_t nLoci = loci.size();
		alleles.reserve(ply * nLoci);

		vectoru chromTypes;

		for (size_t j = 0; j < loci.size(); ++j)
			chromTypes.push_back(chromType(chromLocusPair(loci[j]).first));

		// if there is a small number of loci, or if a small fraction of loci
		// get the loci directly.
		if (nLoci * 100 < totNumLoci()) {
			for (size_t idx = 0; idx < loci.size(); ++idx) {
				for (int p = 0; p < ply; ++p) {
					if (chromTypes[idx] == CHROMOSOME_Y && sex() == FEMALE)
						continue;
					if (((chromTypes[idx] == CHROMOSOME_X && p == 1) ||
					     (chromTypes[idx] == CHROMOSOME_Y && p == 0)) && sex() == MALE)
						continue;
					if (chromTypes[idx] == MITOCHONDRIAL && p > 0)
						continue;
					alleles.push_back(allele(loci[idx], p));
				}
			}
		} else {
			// otherwise put all alleles in a temporary array
			vector<ULONG> tmp_alleles(ply * nLoci, 0);
			std::map<size_t, size_t> loci_map;
			for (size_t i = 0; i < nLoci; ++i)
				loci_map[loci[i]] = i;

			vectorm::const_val_iterator m_ptr = m_genoPtr.get_val_iterator();
			vectorm::const_val_iterator m_end = (m_genoPtr + genoSize()).get_val_iterator();
			for (; m_ptr != m_end; ++m_ptr) {
				size_t loc = m_ptr->first % genoSize();
				size_t p = loc / totNumLoci();
				// find it
				std::map<size_t, size_t>::iterator it = loci_map.find(loc % totNumLoci());
				if (it != loci_map.end())
					tmp_alleles[it->second + p * nLoci] = m_ptr->second;
			}

			for (size_t idx = 0; idx < loci.size(); ++idx) {
				for (int p = 0; p < ply; ++p) {
					if (chromTypes[idx] == CHROMOSOME_Y && sex() == FEMALE)
						continue;
					if (((chromTypes[idx] == CHROMOSOME_X && p == 1) ||
					     (chromTypes[idx] == CHROMOSOME_Y && p == 0)) && sex() == MALE)
						continue;
					if (chromTypes[idx] == MITOCHONDRIAL && p > 0)
						continue;
					alleles.push_back(tmp_alleles[idx + p * nLoci]);
				}
			}
		}
#else
		const vectoru & loci = lociList.elems(this);

		vectoru chromTypes;

		for (size_t j = 0; j < loci.size(); ++j)
			chromTypes.push_back(chromType(chromLocusPair(loci[j]).first));

		alleles.reserve(ply * loci.size());

		for (size_t idx = 0; idx < loci.size(); ++idx) {
			for (int p = 0; p < ply; ++p) {
				if (chromTypes[idx] == CHROMOSOME_Y && sex() == FEMALE)
					continue;
				if (((chromTypes[idx] == CHROMOSOME_X && p == 1) ||
				     (chromTypes[idx] == CHROMOSOME_Y && p == 0)) && sex() == MALE)
					continue;
				if (chromTypes[idx] == MITOCHONDRIAL && p > 0)
					continue;
				alleles.push_back(allele(loci[idx], p));
			}
		}
#endif
	}
	PyObject * genoObj = PyTuple_New(alleles.size());
	// set value
	for (size_t j = 0; j < alleles.size(); ++j)
		PyTuple_SET_ITEM(genoObj, j, PyInt_FromLong(alleles[j]));
	return genoObj;
}


void Individual::setAllele(ULONG allele, size_t idx, int p, int chrom)
{
	DBG_FAILIF(p < 0 && chrom >= 0, ValueError,
		"A valid ploidy index has to be specified if chrom is non-positive");
	if (p < 0) {
		CHECKRANGEGENOSIZE(idx);
#ifdef MUTANTALLELE
		(m_genoPtr + idx).assignIfDiffer(allele);
#else
		*(m_genoPtr + idx) = TO_ALLELE(allele);
#endif
	} else if (chrom < 0) {
		CHECKRANGEABSLOCUS(idx);
		CHECKRANGEPLOIDY(static_cast<size_t>(p));
#ifdef MUTANTALLELE
		(m_genoPtr + idx + p * totNumLoci()).assignIfDiffer(allele);
#else
		*(m_genoPtr + idx + p * totNumLoci()) = TO_ALLELE(allele);
#endif
	} else {
		CHECKRANGELOCUS(static_cast<size_t>(chrom), idx);
		CHECKRANGEPLOIDY(static_cast<size_t>(p));
		CHECKRANGECHROM(static_cast<size_t>(chrom));
#ifdef MUTANTALLELE
		(m_genoPtr + idx + p * totNumLoci() + chromBegin(chrom)).assignIfDiffer(allele);
#else
		*(m_genoPtr + idx + p * totNumLoci() + chromBegin(chrom)) = TO_ALLELE(allele);
#endif
	}
}


void Individual::setAlleleLineage(long lineage, size_t idx, int p, int chrom)
{
#ifdef LINEAGE
	DBG_FAILIF(p < 0 && chrom >= 0, ValueError,
		"A valid ploidy index has to be specified if chrom is non-positive");
	if (p < 0) {
		CHECKRANGEGENOSIZE(idx);
		*(m_lineagePtr + idx) = lineage;
	} else if (chrom < 0) {
		CHECKRANGEABSLOCUS(idx);
		CHECKRANGEPLOIDY(static_cast<size_t>(p));
		*(m_lineagePtr + idx + p * totNumLoci()) = lineage;
	} else {
		CHECKRANGELOCUS(static_cast<size_t>(chrom), idx);
		CHECKRANGEPLOIDY(static_cast<size_t>(p));
		CHECKRANGECHROM(static_cast<size_t>(chrom));
		*(m_lineagePtr + idx + p * totNumLoci() + chromBegin(chrom)) = lineage;
	}
#else
	(void)lineage;
	(void)idx;
	(void)p;
	(void)chrom;
#endif
}


void Individual::setGenotype(const uintList & genoList, const uintList & ply, const uintList & ch)
{
	const vectoru & geno = genoList.elems();

	size_t sz = geno.size();
	size_t idx = 0;

	vectoru ploidys = ply.elems();

	if (ply.allAvail()) {
		for (size_t i = 0; i < ploidy(); ++i)
			ploidys.push_back(i);
	} else {
#ifndef OPTIMIZED
		for (size_t i = 0; i < ploidys.size(); ++i) {
			CHECKRANGEPLOIDY(static_cast<size_t>(ploidys[i]));
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
			CHECKRANGECHROM(static_cast<size_t>(chroms[i]));
		}
#endif
	}
	for (size_t i = 0; i < ploidys.size(); ++i) {
		size_t p = ploidys[i];
		for (size_t j = 0; j < chroms.size(); ++j) {
			size_t chrom = chroms[j];
			GenoIterator ptr = m_genoPtr + p * totNumLoci() + chromBegin(chrom);

#ifdef MUTANTALLELE
			vectorm ctmp;
			ctmp.resize(numLoci(chrom));
			for (size_t i = 0; i < numLoci(chrom); i++, ++idx) {
				Allele atmp = TO_ALLELE(geno[idx % sz]);
				if (atmp != 0)
					ctmp.push_back(i, atmp);
			}
			copyGenotype(ctmp.begin(), ctmp.end(), ptr);
#else
			for (size_t i = 0; i < numLoci(chrom); i++, ++idx)
				*(ptr + i) = TO_ALLELE(geno[idx % sz]);
#endif

		}
	}
}


void Individual::setLineage(const uintList & lineageList, const uintList & ply, const uintList & ch)
{
#ifdef LINEAGE
	const vectoru & lineage = lineageList.elems();

	size_t sz = lineage.size();
	size_t idx = 0;

	vectoru ploidys = ply.elems();

	if (ply.allAvail()) {
		for (size_t i = 0; i < ploidy(); ++i)
			ploidys.push_back(i);
	} else {
#  ifndef OPTIMIZED
		for (size_t i = 0; i < ploidys.size(); ++i) {
			CHECKRANGEPLOIDY(static_cast<size_t>(ploidys[i]));
		}
#  endif
	}
	vectoru chroms = ch.elems();
	if (ch.allAvail()) {
		for (size_t i = 0; i < numChrom(); ++i)
			chroms.push_back(i);
	} else {
#  ifndef OPTIMIZED
		for (size_t i = 0; i < chroms.size(); ++i) {
			CHECKRANGECHROM(static_cast<size_t>(chroms[i]));
		}
#  endif
	}
	for (size_t i = 0; i < ploidys.size(); ++i) {
		size_t p = ploidys[i];
		for (size_t j = 0; j < chroms.size(); ++j) {
			size_t chrom = chroms[j];
			LineageIterator ptr = m_lineagePtr + p * totNumLoci() + chromBegin(chrom);

			for (size_t i = 0; i < numLoci(chrom); i++, ++idx)
				*(ptr + i) = static_cast<long>(lineage[idx % sz]);
		}
	}
#else
	(void)lineageList;
	(void)ply;
	(void)ch;
#endif
}


void Individual::swap(Individual & ind, bool swapContent)
{
	if (genoStruIdx() != ind.genoStruIdx())
		throw SystemError("Can only swap individuals with different geno structure.");

	std::swap(m_infoPtr, ind.m_infoPtr);

	if (swapContent) {
		Allele tmp;
		LINEAGE_EXPR(long tmpLineage);
		for (size_t i = 0, iEnd = genoSize(); i < iEnd; i++) {
			tmp = DEREF_ALLELE(m_genoPtr + i);
			REF_ASSIGN_ALLELE(m_genoPtr + i, DEREF_ALLELE(ind.m_genoPtr + i));
			REF_ASSIGN_ALLELE(ind.m_genoPtr + i, tmp);
		}
#ifdef LINEAGE
		for (size_t i = 0, iEnd = genoSize(); i < iEnd; i++) {
			tmpLineage = m_lineagePtr[i];
			m_lineagePtr[i] = ind.m_lineagePtr[i];
			ind.m_lineagePtr[i] = tmpLineage;
		}
#endif
	} else {
		std::swap(m_genoPtr, ind.m_genoPtr);
		LINEAGE_EXPR(std::swap(m_lineagePtr, ind.m_lineagePtr));
	}
}


void Individual::display(ostream & out, int width, const vectoru & loci, const vectoru & infoIdx)
{
	out << (sex() == MALE ? 'M' : 'F') << (affected() ? 'A' : 'U') << " ";
	int pEnd = ploidy();
	for (int p = 0; p < pEnd; ++p) {
		if (loci.empty()) {
			for (ssize_t ch = 0, chEnd = numChrom(); ch < chEnd; ++ch) {
				for (size_t j = 0, jEnd = numLoci(ch); j < jEnd; ++j)
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
		for (vectoru::const_iterator info = infoIdx.begin();
			info != infoIdx.end(); ++info)
			out << " " << *(infoBegin() + *info);
	}
}


}
