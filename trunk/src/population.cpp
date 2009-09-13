/**
 *  $File: population.cpp $
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

#include "population.h"
#include "pedigree.h"
#include "virtualSubPop.h"

// for file compression
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

namespace simuPOP {

population::population(const uintList & size,
	float ploidy,
	const uintList & loci,
	const uintList & chromTypes,
	const floatList & lociPos,
	int ancGen,
	const stringList & chromNames,
	const stringMatrix & alleleNames,
	const stringList & lociNames,
	const stringList & subPopNames,
	const stringList & infoFields)
	:
	GenoStruTrait(),
	m_popSize(0),
	m_subPopSize(size.elems()),
	m_subPopNames(),
	m_subPopIndex(size.elems().size() + 1),
	m_vspSplitter(NULL),
	m_genotype(0),
	m_info(0),
	m_inds(0),
	m_ancestralGens(ancGen),
	m_vars(NULL, true),
	m_ancestralPops(0),
	m_rep(-1),
	m_gen(0),
	m_curAncestralGen(0),
	m_indOrdered(true),
	m_selectionFlags()
{
	DBG_DO(DBG_POPULATION, cerr << "Constructor of population is called\n");

	DBG_FAILIF(m_subPopSize.size() > MaxSubPopID, ValueError,
		"Number of subpopulations exceed maximum allowed subpopulation numbers");

	// get a GenoStructure with parameters. GenoStructure may be shared by some populations
	// a whole set of functions ploidy() etc in GenoStruTriat can be used after this step.
	DBG_FAILIF(static_cast<UINT>(ploidy) * 1.0 != ploidy && fcmp_ne(ploidy, Haplodiploid),
		ValueError, "Only integer ploidy number or Haplodiploid can be specified");

	setGenoStructure(fcmp_eq(ploidy, Haplodiploid) ? 2 : static_cast<UINT>(ploidy),
		loci.elems(), chromTypes.elems(), fcmp_eq(ploidy, Haplodiploid), lociPos.elems(),
		chromNames.elems(), alleleNames.elems(), lociNames.elems(), infoFields.elems());

	DBG_DO(DBG_DEVEL, cerr << "individual size is " << sizeof(individual) << '+'
		                   << sizeof(Allele) << '*' << genoSize() << endl
		                   << ", infoPtr: " << sizeof(double *)
		                   << ", GenoPtr: " << sizeof(Allele *) << ", Flag: " << sizeof(unsigned char)
		                   << ", plus genoStru"
		                   << "\ngenoSize " << genoSize()
		                   << endl);

	// m_popSize will be defined in fitSubPopStru
	if (m_subPopSize.empty())
		m_subPopSize.resize(1, 0);

	fitSubPopStru(m_subPopSize, subPopNames.elems());
}


population::~population()
{
	DBG_DO(DBG_POPULATION,
		cerr << "Destructor of population is called" << endl);

	if (m_vspSplitter)
		delete m_vspSplitter;

	decGenoStruRef();
}


population::population(const population & rhs) :
	GenoStruTrait(rhs),
	m_popSize(rhs.m_popSize),
	m_subPopSize(rhs.m_subPopSize),
	m_subPopNames(rhs.m_subPopNames),
	m_subPopIndex(rhs.m_subPopIndex),
	m_vspSplitter(NULL),
	m_genotype(0),
	m_info(0),
	m_inds(0),
	m_ancestralGens(rhs.m_ancestralGens),
	m_vars(rhs.m_vars),                                                                     // variables will be copied
	m_rep(-1),                                                                              // rep is set to -1 for new pop (until simulator really set them
	m_gen(0),
	m_curAncestralGen(rhs.m_curAncestralGen),
	m_indOrdered(true),
	m_selectionFlags()
{
	DBG_DO(DBG_POPULATION,
		cerr << "Copy constructor of population is called" << endl);

	try {
		m_inds.resize(rhs.m_popSize);
		m_genotype.resize(m_popSize * genoSize());
		// have 0 length for mpi/non-head node
		m_info.resize(rhs.m_popSize * infoSize());
	} catch (...) {
		throw RuntimeError("Failed to copy population, likely a memory allocation failure.");
	}

	// individuals will always have the correct genostructure
	// by using their copied pointer
	// population, however, need to set this pointer correctly
	//
	setGenoStruIdx(rhs.genoStruIdx());
	incGenoStruRef();   // inc ref to the new

	// copy genotype one by one so individual genoPtr will not
	// point outside of subpopulation region.
	GenoIterator ptr = m_genotype.begin();
	InfoIterator infoPtr = m_info.begin();
	UINT step = genoSize();
	UINT infoStep = infoSize();
	for (ULONG i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		m_inds[i].setGenoPtr(ptr);
		m_inds[i].setInfoPtr(infoPtr);
		m_inds[i].copyFrom(rhs.m_inds[i]);
	}

	// copy ancestral populations
	try {
		// copy all. individual will be shallow copied
		m_ancestralPops = rhs.m_ancestralPops;
		// need to setGenoPtr
		for (size_t ap = 0; ap < m_ancestralPops.size(); ++ap) {
			popData & lp = m_ancestralPops[ap];
			const popData & rp = rhs.m_ancestralPops[ap];

			vector<individual> & linds = lp.m_inds;
			const vector<individual> & rinds = rp.m_inds;

			GenoIterator lg = lp.m_genotype.begin();
			ConstGenoIterator rg = rp.m_genotype.begin();

			InfoIterator li = lp.m_info.begin();
			ConstInfoIterator ri = rp.m_info.begin();

			ULONG ps = rinds.size();

			for (ULONG i = 0; i < ps; ++i) {
				linds[i].setGenoPtr(rinds[i].genoPtr() - rg + lg);
				linds[i].setInfoPtr(rinds[i].infoPtr() - ri + li);
			}
		}
	} catch (...) {
		cerr << "Unable to copy ancestral populations. "
		     << "The popolation size may be too big." << endl
		     << "The population will still be usable but without any ancestral population stored." << endl;
		m_ancestralGens = 0;
		m_ancestralPops.clear();
	}

	// copy virtual subpop splitters
	setVirtualSplitter(rhs.virtualSplitter());
}


void population::popData::swap(population & pop)
{
	pop.m_subPopSize.swap(m_subPopSize);
	pop.m_subPopNames.swap(m_subPopNames);
	pop.m_genotype.swap(m_genotype);
	pop.m_info.swap(m_info);
	pop.m_inds.swap(m_inds);
	std::swap(pop.m_indOrdered, m_indOrdered);
}


population * population::clone() const
{
	return new population(*this);
}


UINT population::subPopByName(const string & name) const
{
	vectorstr::const_iterator it = find(m_subPopNames.begin(), m_subPopNames.end(), name);

	DBG_FAILIF(it == m_subPopNames.end(), ValueError,
		"Subpopulation " + name + " not found.");
	return it - m_subPopNames.begin();
}


string population::subPopName(vspID subPop) const
{
	DBG_ASSERT(m_subPopNames.empty() || m_subPopNames.size() == numSubPop(), SystemError,
		"subpopulation names can either be empty, or be specified for all subpopulations.");
	CHECKRANGESUBPOP(subPop.subPop());
	string name = m_subPopNames.empty() ? UnnamedSubPop : m_subPopNames[subPop.subPop()];;
	if (subPop.isVirtual()) {
		if (name.empty())
			return m_vspSplitter->name(subPop.virtualSubPop());
		else
			return name + " - " + m_vspSplitter->name(subPop.virtualSubPop());
	} else
		return name;
}


vectorstr population::subPopNames() const
{
	DBG_ASSERT(m_subPopNames.empty() || m_subPopNames.size() == numSubPop(), SystemError,
		"subpopulation names can either be empty, or be specified for all subpopulations.");
	return m_subPopNames.empty() ? vectorstr(numSubPop(), UnnamedSubPop) : m_subPopNames;
}


void population::setSubPopName(const string & name, SubPopID subPop)
{
	CHECKRANGESUBPOP(subPop);
	if (m_subPopNames.empty())
		m_subPopNames = vectorstr(numSubPop(), UnnamedSubPop);
	m_subPopNames[subPop] = name;
}


bool population::hasActivatedVirtualSubPop() const
{
	return m_vspSplitter != NULL && m_vspSplitter->activatedSubPop() != InvalidSubPopID;
}


bool population::hasActivatedVirtualSubPop(SubPopID subPop) const
{
	return m_vspSplitter != NULL && m_vspSplitter->activatedSubPop() == subPop;
}


bool population::hasVirtualSubPop() const
{
	return m_vspSplitter != NULL;
}


void population::setVirtualSplitter(vspSplitter * vsp)
{
	if (m_vspSplitter)
		delete m_vspSplitter;

	m_vspSplitter = vsp ? vsp->clone() : NULL;
}


UINT population::numVirtualSubPop() const
{
	return hasVirtualSubPop()
	       ? m_vspSplitter->numVirtualSubPop()
		   : 0;
}


void population::activateVirtualSubPop(vspID subPop,
                                       IterationType type)
{
	CHECKRANGESUBPOP(subPop.subPop());
	DBG_ASSERT(subPop.isVirtual(), ValueError, "Given virtual subpopulation ID is wrong");
	DBG_ASSERT(hasVirtualSubPop(), ValueError,
		"Population has no virtual subpopulations");
	m_vspSplitter->activate(*this, subPop.subPop(), subPop.virtualSubPop(), type);
	DBG_ASSERT(type != VisibleInds ||
		m_vspSplitter->activatedSubPop() == subPop.subPop(), SystemError,
		"Failed to activate virtual subpopulation");
}


void population::deactivateVirtualSubPop(SubPopID subPop)
{
	CHECKRANGESUBPOP(subPop);
	if (!hasActivatedVirtualSubPop(subPop))
		return;
	m_vspSplitter->deactivate(*this, subPop);
}


int population::__cmp__(const population & rhs) const
{
	if (genoStruIdx() != rhs.genoStruIdx() ) {
		DBG_DO(DBG_POPULATION, cerr << "Genotype structures are different" << endl);
		return 1;
	}

	if (popSize() != rhs.popSize() ) {
		DBG_DO(DBG_POPULATION, cerr << "Population sizes are different" << endl);
		return 1;
	}

	if (ancestralGens() != rhs.ancestralGens()) {
		DBG_DO(DBG_POPULATION, cerr << "Number of ancestral generations differ" << endl);
		return 1;
	}

	int curGen = m_curAncestralGen;
	int rhsCurGen = rhs.m_curAncestralGen;
	for (int depth = ancestralGens(); depth >= 0; --depth) {
		const_cast<population *>(this)->useAncestralGen(depth);
		const_cast<population &>(rhs).useAncestralGen(depth);
		for (ULONG i = 0, iEnd = popSize(); i < iEnd; ++i)
			if (m_inds[i] != rhs.m_inds[i]) {
				DBG_DO(DBG_POPULATION, cerr << "Individuals are different" << endl);
				const_cast<population *>(this)->useAncestralGen(curGen);
				const_cast<population &>(rhs).useAncestralGen(rhsCurGen);
				return 1;
			}
	}

	return 0;
}


individual & population::indByID(ULONG id, int ancGen, const string & idField)
{
	UINT idx = infoIdx(idField);
	DBG_FAILIF(ancGen != -1 && static_cast<UINT>(ancGen) > m_ancestralPops.size(), IndexError,
		"Ancestray generation " + toStr(ancGen) + " does not exist");
	
	for (UINT gen = 0; gen <= ancestralGens(); ++gen) {
		// only search specific generation
		if (ancGen != -1 && static_cast<UINT>(ancGen) != gen)
			continue;
		vector<individual> * inds = NULL;
		// search in current, not necessarily the present generation
		if (gen == m_curAncestralGen)
			inds = &m_inds;
		else
			inds = &m_ancestralPops[gen == 0 ? m_curAncestralGen -1 : gen - 1].m_inds;
		// first try our luck
		ULONG startID = (*inds)[0].intInfo(idx);
		if (idx >= startID && startID + (*inds).size() > id) {
			individual & ind = (*inds)[id - startID];
			if (static_cast<ULONG>(ind.intInfo(idx)) == id)
				return ind;
		}
		// now we have to search all individuals
		for (size_t i = 0; i < (*inds).size(); ++i) {
			if (static_cast<ULONG>((*inds)[i].intInfo(idx)) == id)
				return (*inds)[i];
		}
	}
	throw IndexError("No individual with ID " + toStr(id) + " could be found.");
	// this is just to suppress a warning.
	return m_inds[0];
}


individual & population::ancestor(ULONG idx, UINT gen, vspID vsp)
{
	DBG_FAILIF(vsp.isVirtual(), ValueError,
		"Function genotype currently does not support virtual subpopulation");

	DBG_FAILIF(gen > m_ancestralPops.size(), IndexError,
		"Ancestray generation " + toStr(gen) + " does not exist");
	if (!vsp.valid()) {
		if (gen == m_curAncestralGen)
			return this->ind(idx);
		UINT genIdx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
		DBG_FAILIF(idx > m_ancestralPops[genIdx].m_inds.size(),
			IndexError, "Individual index out of range");
		return m_ancestralPops[genIdx].m_inds[idx];
	} else {
		SubPopID subPop = vsp.subPop();
		if (gen == m_curAncestralGen)
			return this->ind(idx, subPop);
		UINT genIdx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
		DBG_FAILIF(static_cast<UINT>(subPop) > m_ancestralPops[genIdx].m_subPopSize.size(),
			IndexError, "subpopulation index out of range");
		DBG_FAILIF(idx > m_ancestralPops[genIdx].m_subPopSize[subPop],
			IndexError, "Individual index out of range");
		ULONG shift = 0;
		if (subPop > 0) {
			for (int i = 0; i < subPop; ++i)
				shift += m_ancestralPops[genIdx].m_subPopSize[i];
		}
		return m_ancestralPops[genIdx].m_inds[shift + idx];
	}
}


const individual & population::ancestor(ULONG idx, UINT gen, vspID vsp) const
{
	DBG_FAILIF(vsp.isVirtual(), ValueError,
		"Function genotype currently does not support virtual subpopulation");

	DBG_FAILIF(gen > m_ancestralPops.size(), IndexError,
		"Ancestray generation " + toStr(gen) + " does not exist");
	if (!vsp.valid()) {
		if (gen == m_curAncestralGen)
			return this->ind(idx);
		UINT genIdx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
		DBG_FAILIF(idx > m_ancestralPops[genIdx].m_inds.size(),
			IndexError, "Individual index out of range");
		return m_ancestralPops[genIdx].m_inds[idx];
	} else {
		SubPopID subPop = vsp.subPop();
		if (gen == m_curAncestralGen)
			return this->ind(idx, subPop);
		UINT genIdx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
		DBG_FAILIF(static_cast<UINT>(subPop) > m_ancestralPops[genIdx].m_subPopSize.size(),
			IndexError, "subpopulation index out of range");
		DBG_FAILIF(idx > m_ancestralPops[genIdx].m_subPopSize[subPop],
			IndexError, "Individual index out of range");
		ULONG shift = 0;
		if (subPop > 0) {
			for (int i = 0; i < subPop; ++i)
				shift += m_ancestralPops[genIdx].m_subPopSize[i];
		}
		return m_ancestralPops[genIdx].m_inds[shift + idx];
	}
}


IndAlleleIterator population::alleleIterator(UINT locus)
{
	CHECKRANGEABSLOCUS(locus);

	// if there is virtual subpop, use individual based iterator
	// or
	// if requires order, but the alleles are not ordered
	// use individual based
	int ct = chromType(chromLocusPair(locus).first);
	if (hasActivatedVirtualSubPop() || !indOrdered()
	    || (ct != Autosome && ct != Customized) || isHaplodiploid())
		// this is a complex case
		return IndAlleleIterator(locus, indIterator());
	else
		// a simplere case with straightforward iterator
		return IndAlleleIterator(m_genotype.begin() + locus, m_genotype.end() + locus,
			totNumLoci());
}


/// CPPONLY allele begin, for given subPop
IndAlleleIterator population::alleleIterator(UINT locus, UINT subPop)
{
	CHECKRANGEABSLOCUS(locus);
	CHECKRANGESUBPOP(subPop);

	int ct = chromType(chromLocusPair(locus).first);
	// this is a complex case
	if (hasActivatedVirtualSubPop() || !indOrdered()
	    || (ct != Autosome && ct != Customized) || isHaplodiploid())
		// this is a complex case
		return IndAlleleIterator(locus, indIterator(subPop));
	else
		// this is a complex case
		return IndAlleleIterator(
			m_genotype.begin() + m_subPopIndex[subPop] * genoSize() + locus,
			m_genotype.begin() + m_subPopIndex[subPop + 1] * genoSize() + locus,
			totNumLoci());
}


PyObject * population::genotype(vspID vsp)
{
	DBG_FAILIF(vsp.isVirtual(), ValueError,
		"Function genotype currently does not support virtual subpopulation");
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	sortIndividuals();
	if (!vsp.valid()) {
		// directly expose values. Do not copy data over.
		return Allele_Vec_As_NumArray(m_genotype.begin(), m_genotype.end());
	} else {
		SubPopID subPop = vsp.subPop();
		CHECKRANGESUBPOP(subPop);
		// directly expose values. Do not copy data over.
		return Allele_Vec_As_NumArray(genoBegin(subPop, true), genoEnd(subPop, true));
	}
	return NULL;
}


void population::setGenotype(const vectora & geno, vspID subPop)
{
	sortIndividuals();
	if (!subPop.valid()) {
		GenoIterator ptr = m_genotype.begin();
		ULONG sz = geno.size();
		for (ULONG i = 0; i < popSize() * genoSize(); ++i)
			*(ptr++) = geno[i % sz];
		return;
	}

	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	SubPopID sp = subPop.subPop();
	CHECKRANGESUBPOP(sp);

	ULONG sz = geno.size();
	if (!subPop.isVirtual()) {
		GenoIterator ptr = genoBegin(sp, true);
		for (ULONG i = 0; i < subPopSize(sp) * genoSize(); ++i)
			*(ptr++) = geno[i % sz];
	} else {
		activateVirtualSubPop(subPop, IteratableInds);
		IndIterator it = indIterator(sp, IteratableInds);
		ULONG i = 0;
		for (; it.valid(); ++it)
			for (GenoIterator git = it->genoBegin(); git != it->genoEnd(); ++git, ++i)
				*git = geno[i % sz];
	}
}


void population::validate(const string & msg) const
{
#ifndef OPTIMIZED
	DBG_ASSERT(m_info.size() == m_popSize * infoSize(), SystemError,
		msg + "Wrong information size");
	DBG_ASSERT(m_genotype.size() == m_popSize * genoSize(), SystemError,
		msg + "Wrong genotype size for this population");
	ConstInfoIterator ib = m_info.begin();
	ConstInfoIterator ie = m_info.end();
	ConstGenoIterator gb = m_genotype.begin();
	ConstGenoIterator ge = m_genotype.end();

	if (genoSize() > 0) {
		for (ConstIndIterator it = indIterator(); it.valid(); ++it) {
			DBG_ASSERT(it->genoPtr() >= gb && it->genoPtr() < ge, SystemError,
				msg + "Wrong genotype pointer");
		}
	}
	if (infoSize() > 0) {
		for (ConstIndIterator it = indIterator(); it.valid(); ++it) {
			DBG_ASSERT(it->infoPtr() >= ib && it->infoPtr() < ie, SystemError,
				msg + "Wrong information field pointer. (number of information fields: "
				+ toStr(infoSize()) + ")");
		}
	}
#endif
}


void population::fitSubPopStru(const vectoru & newSubPopSizes,
                               const vectorstr & newSubPopNames)
{
	ULONG newSize = accumulate(newSubPopSizes.begin(), newSubPopSizes.end(), 0UL);

	bool needsResize = m_popSize != newSize;

	if (needsResize) {
		UINT is = infoSize();
		UINT step = genoSize();
		m_popSize = newSize;
		try {
			m_genotype.resize(m_popSize * step);
			m_info.resize(m_popSize * is);
			m_inds.resize(m_popSize);
		} catch (...) {
			throw RuntimeError("Failed to create population (popSize=" + toStr(m_popSize) + ")");
		}
		// reset individual pointers
		GenoIterator ptr = m_genotype.begin();
		InfoIterator infoPtr = m_info.begin();
		for (ULONG i = 0; i < m_popSize; ++i, ptr += step, infoPtr += is) {
			m_inds[i].setGenoPtr(ptr);
			m_inds[i].setInfoPtr(infoPtr);
			m_inds[i].setGenoStruIdx(genoStruIdx());
		}
		setIndOrdered(true);
	}
	// help clear confusing
	std::fill(m_info.begin(), m_info.end(), 0.);

	if (newSubPopNames.empty() || newSubPopNames.size() == newSubPopSizes.size())
		setSubPopStru(newSubPopSizes, newSubPopNames);
	else {
		vectorstr spNames = newSubPopNames;
		spNames.resize(newSubPopSizes.size(), UnnamedSubPop);
		setSubPopStru(newSubPopSizes, spNames);
	}
}


void population::fitGenoStru(size_t stru)
{
	// set genotypic structure to a population.
	// This function will try not to change population size.
	UINT oldSize = genoSize();
	UINT oldInfoSize = infoSize();

	decGenoStruRef();   // dec ref to the old
	setGenoStruIdx(stru);
	incGenoStruRef();   // inc ref to the new
	UINT newSize = genoSize();
	UINT newInfoSize = infoSize();

	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		if (oldSize != newSize)
			m_genotype.resize(newSize * popSize());
		if (oldInfoSize != newInfoSize)
			m_info.resize(newInfoSize * popSize());
		// reset structure
		GenoIterator ptr = m_genotype.begin();
		InfoIterator infoPtr = m_info.begin();
		for (ULONG i = 0; i < m_popSize; ++i, ptr += newSize, infoPtr += newInfoSize) {
			m_inds[i].setGenoStruIdx(stru);
			m_inds[i].setGenoPtr(ptr);
			m_inds[i].setInfoPtr(infoPtr);
		}
	}
}


void population::setSubPopStru(const vectoru & newSubPopSizes,
                               const vectorstr & newSubPopNames)
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	DBG_ASSERT(accumulate(newSubPopSizes.begin(), newSubPopSizes.end(), 0UL) == m_popSize, ValueError,
		"Overall population size should not be changed in setSubPopStru.");

	DBG_ASSERT(newSubPopNames.empty() || newSubPopNames.size() == newSubPopSizes.size(), SystemError,
		"subpopulation names can either be empty, or be specified for all subpopulations.");

	if (newSubPopSizes.empty())
		m_subPopSize = vectoru(1, 0);
	else
		m_subPopSize = newSubPopSizes;
	m_subPopIndex.resize(numSubPop() + 1);
	m_subPopNames = newSubPopNames;

	// build subPop index
	UINT i = 1;
	for (m_subPopIndex[0] = 0; i <= numSubPop(); ++i)
		m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];
}


void population::setSubPopByIndInfo(const string & field)
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	UINT info = infoIdx(field);

	DBG_DO(DBG_POPULATION, cerr << "Sorting individuals." << endl);
	// sort individuals first
	std::sort(indIterator(), IndIterator(m_inds.end(), m_inds.end(), AllInds), indCompare(info));
	setIndOrdered(false);

	// sort individuals first
	// remove individuals with negative index.
	if (indIterator()->info(info) < 0) {
		// popsize etc will be changed.
		ULONG newPopSize = m_popSize;
		IndIterator it = indIterator();
		for (; it.valid();  ++it) {
			if (it->info(info) < 0)
				newPopSize-- ;
			else
				break;
		}
		// 'it' now point to the one with positive info(info)

		DBG_DO(DBG_POPULATION, cerr << "New pop size" << newPopSize << endl);

		// allocate new genotype and inds
		vectora newGenotype(genoSize() * newPopSize);
		vectorinfo newInfo(newPopSize * infoSize());
		vector<individual> newInds(newPopSize);

		// assign genotype location and set structure information for individuals
		GenoIterator ptr = newGenotype.begin();
		InfoIterator infoPtr = newInfo.begin();
		UINT step = genoSize();
		UINT infoStep = infoSize();
		for (ULONG i = 0; i < newPopSize; ++i, ptr += step, ++it, infoPtr += infoStep) {
			newInds[i].setGenoStruIdx(genoStruIdx());
			newInds[i].setGenoPtr(ptr);
			newInds[i].setInfoPtr(infoPtr);
			newInds[i].copyFrom(*it);                         // copy everything, with info value
		}

		// now, switch!
		m_genotype.swap(newGenotype);
		m_info.swap(newInfo);
		m_inds.swap(newInds);

		m_popSize = newPopSize;
		setIndOrdered(true);
	}

	if (m_inds.empty()) {
		m_subPopSize.resize(1, 0);
		m_subPopIndex.resize(2);
	} else {
		// reset indexes etc.
		UINT numSubPop = static_cast<UINT>(m_inds.back().info(info)) + 1;
		m_subPopSize.resize(numSubPop);
		m_subPopIndex.resize(numSubPop + 1);

		// check subpop size
		fill(m_subPopSize.begin(), m_subPopSize.end(), 0);
		for (IndIterator it = indIterator(); it.valid();  ++it)
			m_subPopSize[ static_cast<UINT>(it->info(info)) ]++;
	}
	// rebuild index
	size_t i = 1;
	for (m_subPopIndex[0] = 0; i <= numSubPop(); ++i)
		m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];
	// subpopulation names
	if (!m_subPopNames.empty())
		m_subPopNames.resize(numSubPop(), UnnamedSubPop);
}


vectoru population::splitSubPop(UINT subPop, const vectoru & sizes, const vectorstr & names)
{
	if (sizes.size() <= 1)
		return vectoru(1, subPop);

	DBG_FAILIF(accumulate(sizes.begin(), sizes.end(), 0LU) != subPopSize(subPop), ValueError,
		"Sum of parameter sizes should be the size of subpopulation " + toStr(subPop));

	DBG_ASSERT(names.empty() || sizes.size() == names.size(), ValueError,
		"Names should be given to none or all of the split subpopulations");

	if (!names.empty() && m_subPopNames.empty())
		m_subPopNames.resize(numSubPop(), UnnamedSubPop);

	vectoru subPopSizes;
	vectorstr subPopNames;
	vectoru ret(sizes.size());
	for (size_t sp = 0; sp < numSubPop(); ++sp) {
		if (sp != subPop) {
			subPopSizes.push_back(subPopSize(sp));
			if (!m_subPopNames.empty())
				subPopNames.push_back(m_subPopNames[sp]);
		} else {
			subPopSizes.insert(subPopSizes.end(), sizes.begin(), sizes.end());
			if (!m_subPopNames.empty()) {
				if (names.empty()) {
					for (size_t i = 0; i < sizes.size(); ++i)
						subPopNames.push_back(m_subPopNames[subPop]);
				} else
					subPopNames.insert(subPopNames.end(), names.begin(), names.end());
			}
			for (size_t i = 0; i < sizes.size(); ++i)
				ret[i] = sp + i;
		}
	}
	setSubPopStru(subPopSizes, subPopNames);
	return ret;
}


void population::removeSubPops(const uintList & subPopList)
{
	const vectoru & subPops = subPopList.elems();

#ifndef OPTIMIZED
	// check if subPops are valid
	for (vectoru::const_iterator sp = subPops.begin(); sp < subPops.end(); ++sp) {
		DBG_FAILIF(*sp >= numSubPop(), IndexError, "Subpopulation " + toStr(*sp) + " does not exist.");
	}
#endif
	sortIndividuals();
	vectoru new_size;
	vectorstr new_spNames;

	UINT step = genoSize();
	UINT infoStep = infoSize();
	vector<individual>::iterator oldInd = m_inds.begin();
	vector<individual>::iterator newInd = m_inds.begin();
	GenoIterator oldPtr = m_genotype.begin();
	GenoIterator newPtr = m_genotype.begin();
	InfoIterator oldInfoPtr = m_info.begin();
	InfoIterator newInfoPtr = m_info.begin();

	for (size_t sp = 0; sp < numSubPop(); ++sp) {
		ULONG spSize = subPopSize(sp);
		// do not remove
		if (find(subPops.begin(), subPops.end(), sp) == subPops.end()) {
			new_size.push_back(spSize);
			if (!m_subPopNames.empty())
				new_spNames.push_back(m_subPopNames[sp]);
			// do not remove.
			if (oldPtr != newPtr) {
				copy(oldInd, oldInd + spSize, newInd);
				copy(oldPtr, oldPtr + step * spSize, newPtr);
				copy(oldInfoPtr, oldInfoPtr + infoStep * spSize, newInfoPtr);
			}
			newInd += spSize;
			newPtr += step * spSize;
			newInfoPtr += infoStep * spSize;
		}
		oldInd += spSize;
		oldPtr += step * spSize;
		oldInfoPtr += infoStep * spSize;
	}
	//
	m_inds.erase(newInd, m_inds.end());
	m_genotype.erase(newPtr, m_genotype.end());
	m_info.erase(newInfoPtr, m_info.end());
	m_popSize = std::accumulate(new_size.begin(), new_size.end(), 0UL);
	setSubPopStru(new_size, new_spNames);
	//
	GenoIterator ptr = m_genotype.begin();
	InfoIterator infoPtr = m_info.begin();
	for (ULONG i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		m_inds[i].setGenoPtr(ptr);
		m_inds[i].setInfoPtr(infoPtr);
	}
}


void population::removeIndividuals(const uintList & indList)
{
	const vectoru & inds = indList.elems();

	sortIndividuals();
	vectoru new_size = m_subPopSize;

	UINT step = genoSize();
	UINT infoStep = infoSize();
	vector<individual>::iterator oldInd = m_inds.begin();
	vector<individual>::iterator newInd = m_inds.begin();
	GenoIterator oldPtr = m_genotype.begin();
	InfoIterator oldInfoPtr = m_info.begin();
	GenoIterator newPtr = m_genotype.begin();
	InfoIterator newInfoPtr = m_info.begin();

	// which ones are removed?
	vector<bool> removed(popSize(), false);
	vectoru::const_iterator it = inds.begin();
	vectoru::const_iterator it_end = inds.end();
	for (; it != it_end; ++it) {
		DBG_FAILIF(*it >= removed.size(), IndexError,
			"Individual index " + toStr(*it) + " out of range.");
		removed[*it] = true;
	}
	//
	for (size_t sp = 0; sp < numSubPop(); ++sp) {
		for (size_t ind = subPopBegin(sp); ind != subPopEnd(sp); ++ind) {
			if (removed[ind])
				--new_size[sp];
			else {
				// do not remove.
				if (oldPtr != newPtr) {
					copy(oldInd, oldInd + 1, newInd);
					copy(oldPtr, oldPtr + step, newPtr);
					copy(oldInfoPtr, oldInfoPtr + infoStep, newInfoPtr);
				}
				++newInd;
				newPtr += step;
				newInfoPtr += infoStep;
			}
			++oldInd;
			oldPtr += step;
			oldInfoPtr += infoStep;
		}
	}
	//
	m_inds.erase(newInd, m_inds.end());
	m_genotype.erase(newPtr, m_genotype.end());
	m_info.erase(newInfoPtr, m_info.end());
	m_popSize = std::accumulate(new_size.begin(), new_size.end(), 0UL);
	setSubPopStru(new_size, m_subPopNames);
	//
	GenoIterator ptr = m_genotype.begin();
	InfoIterator infoPtr = m_info.begin();
	for (ULONG i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		m_inds[i].setGenoPtr(ptr);
		m_inds[i].setInfoPtr(infoPtr);
	}
}


UINT population::mergeSubPops(const vectoru & subPops, const string & name)
{
	if (!name.empty() && m_subPopNames.empty())
		m_subPopNames.resize(numSubPop(), UnnamedSubPop);

	// merge all subpopulations
	if (subPops.empty()) {
		// [ popSize() ]
		vectoru sz(1, popSize());
		if (m_subPopNames.empty())
			setSubPopStru(sz, m_subPopNames);
		else
			setSubPopStru(sz, vectorstr(1, name.empty() ? m_subPopNames[0] : name));
		return 0;
	}
	if (subPops.size() == 1) {
		if (!name.empty())
			m_subPopNames[0] = name;
		return subPops[0];
	}

	// are they in order?
	bool consecutive = true;
	vectoru sps = subPops;
	sort(sps.begin(), sps.end());
	for (size_t i = 1; i < sps.size(); ++i)
		if (subPops[i] != subPops[i - 1] + 1) {
			consecutive = false;
			break;
		}
	// new subpopulation sizes and names
	vectoru new_size;
	vectorstr new_names;
	for (UINT sp = 0; sp < numSubPop(); ++sp) {
		if (find(subPops.begin(), subPops.end(), sp) != subPops.end()) {
			if (new_size.size() <= sps[0]) {
				new_size.push_back(subPopSize(sp));
				if (!m_subPopNames.empty())
					new_names.push_back(name.empty() ? m_subPopNames[sp] : name);
			} else
				new_size[sps[0]] += subPopSize(sp);
		} else {
			new_size.push_back(subPopSize(sp));
			if (!m_subPopNames.empty())
				new_names.push_back(m_subPopNames[sp]);
		}
	}
	// if consecutive, no need to move anyone
	if (consecutive) {
		setSubPopStru(new_size, new_names);
		return sps[0];
	}
	// difficult case.
	sortIndividuals();
	// find the new subpop order
	vectoru sp_order;
	for (size_t sp = 0; sp < sps[0]; ++sp)
		sp_order.push_back(sp);
	sp_order.insert(sp_order.end(), sps.begin(), sps.end());
	for (size_t sp = sps[0]; sp < numSubPop(); ++sp)
		if (find(subPops.begin(), subPops.end(), sp) == subPops.end())
			sp_order.push_back(sp);
	//
	DBG_ASSERT(sp_order.size() == numSubPop(), ValueError,
		"Incorrect resulting subpopulation number, maybe caused by duplicate entries in parameter subPops.");

	UINT step = genoSize();
	UINT infoStep = infoSize();
	vector<individual> new_inds;
	vectora new_genotype;
	vectorinfo new_info;
	new_inds.reserve(popSize());
	new_genotype.reserve(step * popSize());
	new_info.reserve(infoStep * popSize());

	for (size_t sp = 0; sp < numSubPop(); ++sp) {
		size_t src = sp_order[sp];
		if (subPopSize(src) == 0)
			continue;
		// do not remove.
		new_inds.insert(new_inds.end(), rawIndBegin(src), rawIndEnd(src));
		new_genotype.insert(new_genotype.end(), genoBegin(src, true), genoEnd(src, true));
		if (infoStep > 0)
			new_info.insert(new_info.end(),
				m_info.begin() + subPopBegin(src) * infoStep,
				m_info.begin() + (subPopEnd(src) - 1) * infoStep + 1);
	}
	//
	DBG_ASSERT(new_inds.size() == popSize() && new_genotype.size() == popSize() * step
		&& new_info.size() == popSize() * infoStep, SystemError,
		"Incorrect individual manipulation");
	//
	m_inds.swap(new_inds);
	m_genotype.swap(new_genotype);
	m_info.swap(new_info);
	setSubPopStru(new_size, new_names);
	//
	GenoIterator ptr = m_genotype.begin();
	InfoIterator infoPtr = m_info.begin();
	for (ULONG i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		m_inds[i].setGenoPtr(ptr);
		m_inds[i].setInfoPtr(infoPtr);
	}
	return sps[0];
}


void population::addChromFrom(const population & pop)
{
	UINT numLoci1 = totNumLoci();
	UINT numLoci2 = pop.totNumLoci();

	// obtain new genotype structure and set it
	setGenoStructure(gsAddChromFromStru(pop.genoStruIdx()));
	//
	DBG_FAILIF(ancestralGens() != pop.ancestralGens(), ValueError,
		"Can not add chromosomes from a population with different number of ancestral generations");
	//
	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		const_cast<population &>(pop).useAncestralGen(depth);
		//
		DBG_FAILIF(m_subPopSize != pop.m_subPopSize, ValueError,
			"Can not add chromosomes from a population with different subpopulation sizes");

		vectora newGenotype(genoSize() * m_popSize);

		// append pop2 chromosomes to the first one
		GenoIterator ptr = newGenotype.begin();
		UINT pEnd = ploidy();
		for (ULONG i = 0; i < m_popSize; ++i) {
			// set new geno structure
			m_inds[i].setGenoStruIdx(genoStruIdx());
			GenoIterator ptr1 = m_inds[i].genoPtr();
			GenoIterator ptr2 = pop.m_inds[i].genoPtr();
			m_inds[i].setGenoPtr(ptr);
			for (UINT p = 0; p < pEnd; ++p) {
				for (size_t j = 0; j < numLoci1; ++j)
					*(ptr++) = *(ptr1++);
				for (size_t j = 0; j < numLoci2; ++j)
					*(ptr++) = *(ptr2++);
			}
		}
		m_genotype.swap(newGenotype);
	}
	if (!indOrdered())
		// sort information only
		sortIndividuals(true);
}


void population::addIndFrom(const population & pop)
{
	DBG_FAILIF(genoStruIdx() != pop.genoStruIdx(), ValueError,
		"Cannot add individual from a population with different genotypic structure.");
	DBG_FAILIF(ancestralGens() != pop.ancestralGens(), ValueError,
		"Two populations should have the same number of ancestral generations.");
	// genotype pointers may be reset so this is needed.
	sortIndividuals();
	const_cast<population &>(pop).sortIndividuals();
	// go to the oldest generation
	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		const_cast<population &>(pop).useAncestralGen(depth);
		// calculate new population size
		m_subPopSize.insert(m_subPopSize.end(), pop.m_subPopSize.begin(), pop.m_subPopSize.end());
		// new population size
		m_popSize += pop.m_popSize;
		//
		m_inds.insert(m_inds.end(), pop.m_inds.begin(), pop.m_inds.end());
		m_genotype.insert(m_genotype.end(), pop.m_genotype.begin(), pop.m_genotype.end());
		m_info.insert(m_info.end(), pop.m_info.begin(), pop.m_info.end());
		// iterators ready
		GenoIterator ptr = m_genotype.begin();
		InfoIterator infoPtr = m_info.begin();
		UINT step = genoSize();
		UINT infoStep = infoSize();
		// set pointers
		for (ULONG i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
			m_inds[i].setGenoStruIdx(genoStruIdx());
			m_inds[i].setGenoPtr(ptr);
			m_inds[i].setInfoPtr(infoPtr);
		}
		// rebuild index
		m_subPopIndex.resize(numSubPop() + 1);
		size_t j = 1;
		for (m_subPopIndex[0] = 0; j <= numSubPop(); ++j)
			m_subPopIndex[j] = m_subPopIndex[j - 1] + m_subPopSize[j - 1];
	}
	if (!m_subPopNames.empty() && pop.m_subPopNames.empty()) {
		for (size_t i = 0; i < pop.numSubPop(); ++i)
			m_subPopNames.push_back(UnnamedSubPop);
	} else if (m_subPopNames.empty() && !pop.m_subPopNames.empty()) {
		m_subPopNames.resize(numSubPop(), UnnamedSubPop);
		m_subPopNames.insert(m_subPopNames.end(),
			pop.m_subPopNames.begin(), pop.m_subPopNames.end());
	} else {
		m_subPopNames.insert(m_subPopNames.end(),
			pop.m_subPopNames.begin(), pop.m_subPopNames.end());
	}
	DBG_ASSERT(m_subPopNames.empty() || m_subPopNames.size() == numSubPop(), SystemError,
		"subpopulation names can either be empty, or be specified for all subpopulations.");
}


void population::addLociFrom(const population & pop)
{
	DBG_FAILIF(ancestralGens() != pop.ancestralGens(), ValueError,
		"Can not add chromosomes from a population with different number of ancestral generations");

	UINT size1 = totNumLoci();
	UINT size2 = pop.totNumLoci();
	// obtain new genotype structure and set it
	vectoru indexes1;
	vectoru indexes2;
	setGenoStructure(gsAddLociFromStru(pop.genoStruIdx(),
			indexes1, indexes2));

	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		const_cast<population &>(pop).useAncestralGen(depth);
		//
		DBG_FAILIF(m_subPopSize != pop.m_subPopSize, ValueError,
			"Can not add chromosomes from a population with different subpopulation sizes");
		//
		vectora newGenotype(genoSize() * m_popSize);

		// merge chromosome by chromosome
		GenoIterator ptr = newGenotype.begin();
		UINT pEnd = ploidy();
		UINT newSize = totNumLoci();
		for (ULONG i = 0; i < m_popSize; ++i) {
			// set new geno structure
			m_inds[i].setGenoStruIdx(genoStruIdx());
			GenoIterator ptr1 = m_inds[i].genoPtr();
			GenoIterator ptr2 = pop.m_inds[i].genoPtr();
			// new genotype
			m_inds[i].setGenoPtr(ptr);
			// copy each allele
			for (UINT p = 0; p < pEnd; ++p) {
				for (size_t i = 0; i < size1; ++i)
					ptr[indexes1[i]] = *(ptr1++);
				for (size_t i = 0; i < size2; ++i)
					ptr[indexes2[i]] = *(ptr2++);
				ptr += newSize;
			}
		}
		m_genotype.swap(newGenotype);
	}

	// sort information only
	sortIndividuals(true);
}


void population::addChrom(const vectorf & lociPos, const vectorstr & lociNames,
                          const string & chromName, const stringMatrix & alleleNames,
                          UINT chromType)
{
	DBG_ASSERT(lociNames.empty() || lociPos.size() == lociNames.size(), ValueError,
		"Please specifiy locus name for all inserted loci.");

	size_t oldNumLoci = totNumLoci();
	// obtain new genotype structure and set it
	setGenoStructure(gsAddChrom(lociPos, lociNames, chromName, alleleNames.elems(), chromType));

	DBG_ASSERT(totNumLoci() - oldNumLoci == lociPos.size(), SystemError,
		"Failed to add chromosome.");

	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		//
		ULONG newPopGenoSize = genoSize() * m_popSize;
		vectora newGenotype(newPopGenoSize, 0);

		// copy data over
		GenoIterator newPtr = newGenotype.begin();
		UINT pEnd = ploidy();
		UINT gap = totNumLoci() - oldNumLoci;
		for (ULONG i = 0; i < m_popSize; ++i) {
			// set new geno structure
			m_inds[i].setGenoStruIdx(genoStruIdx());
			GenoIterator oldPtr = m_inds[i].genoPtr();
			// new genotype
			m_inds[i].setGenoPtr(newPtr);
			// copy each chromosome
			for (UINT p = 0; p < pEnd; ++p) {
				for (size_t i = 0; i < oldNumLoci; ++i)
					*(newPtr++) = *(oldPtr++);
				newPtr += gap;
			}
		}
		m_genotype.swap(newGenotype);
	}
	// if indOrdered is false:
	//   individual genotype is now sorted. If we do not do
	//   anything, genotype may be resorted. Sort info to
	//   so that the order is set to True.
	sortIndividuals(true);
}


vectoru population::addLoci(const uintList & chromList, const floatList & posList,
                            const vectorstr & lociNames, const stringMatrix & alleleNamesMatrix)
{
	const vectoru & chrom = chromList.elems();
	const vectorf & pos = posList.elems();
	const matrixstr & alleleNames = alleleNamesMatrix.elems();

	DBG_ASSERT(chrom.size() == pos.size(), ValueError,
		"Chromosome and position lists should have the same length");
	DBG_ASSERT(lociNames.empty() || pos.size() == lociNames.size(), ValueError,
		"Please specifiy locus name for all inserted loci.");

	vectoru newIndex;
	vectoru loci(totNumLoci());
	// obtain new genotype structure and set it
	setGenoStructure(gsAddLoci(chrom, pos, lociNames, alleleNames, newIndex));
	DBG_DO(DBG_POPULATION, cerr << "Indexes of inserted loci " << newIndex << endl);
	// loci at newIndex should have zero alleles...
	for (size_t i = 0, j = 0; j < totNumLoci(); ++j) {
		if (find(newIndex.begin(), newIndex.end(), j) == newIndex.end())
			loci[i++] = j;
	}
	DBG_DO(DBG_POPULATION, cerr << "Indexes of inserted loci " << newIndex
		                        << "\nIndexes of old loci in new structure " << loci << endl);

	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		//
		ULONG newPopGenoSize = genoSize() * m_popSize;
		vectora newGenotype(newPopGenoSize, 0);

		// copy data over
		GenoIterator newPtr = newGenotype.begin();
		UINT pEnd = ploidy();
		for (ULONG i = 0; i < m_popSize; ++i) {
			// set new geno structure
			m_inds[i].setGenoStruIdx(genoStruIdx());
			GenoIterator oldPtr = m_inds[i].genoPtr();
			// new genotype
			m_inds[i].setGenoPtr(newPtr);
			// copy each chromosome
			for (UINT p = 0; p < pEnd; ++p) {
				vectoru::iterator loc = loci.begin();
				for (; loc != loci.end(); ++loc)
					newPtr[*loc] = *(oldPtr++);
				newPtr += totNumLoci();
			}
		}
		m_genotype.swap(newGenotype);
	}
	// if indOrdered is false:
	//   individual genotype is now sorted. If we do not do
	//   anything, genotype may be resorted. Sort info to
	//   so that the order is set to True.
	sortIndividuals(true);
	return newIndex;
}


void population::resize(const uintList & sizeList, bool propagate)
{
	const vectoru & newSubPopSizes = sizeList.elems();

	DBG_FAILIF(newSubPopSizes.size() != numSubPop(), ValueError,
		"Resize should give subpopulation size for each subpopulation");

	ULONG newPopSize = accumulate(newSubPopSizes.begin(), newSubPopSizes.end(), 0UL);

	// prepare new population
	vector<individual> newInds(newPopSize);
	vectora newGenotype(genoSize() * newPopSize);
	vectorinfo newInfo(newPopSize * infoSize());
	// iterators ready
	GenoIterator ptr = newGenotype.begin();
	InfoIterator infoPtr = newInfo.begin();
	UINT step = genoSize();
	UINT infoStep = infoSize();
	// set pointers
	for (ULONG i = 0; i < newPopSize; ++i, ptr += step, infoPtr += infoStep) {
		newInds[i].setGenoStruIdx(genoStruIdx());
		newInds[i].setGenoPtr(ptr);
		newInds[i].setInfoPtr(infoPtr);
	}
	// copy stuff over
	ULONG startSP = 0;
	for (UINT sp = 0; sp < numSubPop(); ++sp) {
		ULONG spSize = subPopSize(sp);
		for (ULONG i = 0, j = 0; i < newSubPopSizes[sp]; ++j, ++i) {
			// repeating?
			if (spSize == 0 || ((j / spSize) > 0 && !propagate))
				break;
			newInds[startSP + i].copyFrom(ind(j % spSize, sp));
		}
		// point to the start of next subpopulation
		startSP += newSubPopSizes[sp];
	}
	// now, switch!
	m_genotype.swap(newGenotype);
	m_info.swap(newInfo);
	m_inds.swap(newInds);
	m_popSize = newPopSize;
	setIndOrdered(true);
	m_subPopSize = newSubPopSizes;
	// rebuild index
	size_t idx = 1;
	for (m_subPopIndex[0] = 0; idx <= numSubPop(); ++idx)
		m_subPopIndex[idx] = m_subPopIndex[idx - 1] + m_subPopSize[idx - 1];
}


population & population::extract(const string & field,
                                 const uintList & extractedLoci,
                                 const stringList & infoFieldList,
                                 int ancGen, pedigree * ped, const vectorstr & pedFields) const
{
	population & pop = *new population();
	bool removeInd = !field.empty();
	bool removeLoci = !extractedLoci.allAvail();
	const vectoru & loci = extractedLoci.elems();
	bool removeInfo = !infoFieldList.allAvail();
	const vectorstr & infoFields = infoFieldList.elems();

	DBG_DO(DBG_POPULATION, cerr << "Remove ind: " << removeInd
		                        << "\nRemove loci: " << removeLoci
		                        << "\nRemove info: " << removeInfo << endl);

	// will keep a sorted version of loci
	vectoru new_loci = loci;

	vectorstr keptInfoFields = removeLoci ? infoFields : this->infoFields();
	vectorstr allInfoFields(keptInfoFields.begin(), keptInfoFields.end());
	allInfoFields.insert(allInfoFields.end(), pedFields.begin(), pedFields.end());

	// population strcture.
	if (!removeLoci && !removeInfo && pedFields.empty()) {
		pop.setGenoStruIdx(genoStruIdx());
		incGenoStruRef();
	} else if (!removeLoci) {
		// only change information fields
		pop.setGenoStructure(ploidy(), numLoci(), chromTypes(), isHaplodiploid(),
			lociPos(), chromNames(), allAlleleNames(), lociNames(), allInfoFields);
	} else {
		// figure out number of loci.
		vectoru new_numLoci;
		vectorf new_lociPos;
		vectorstr new_lociNames;
		matrixstr new_alleleNames;
		vectoru new_chromTypes;
		vectorstr new_chromNames;
		if (removeLoci && !loci.empty()) {
			sort(new_loci.begin(), new_loci.end());
			vectoru::const_iterator it = new_loci.begin();
			vectoru::const_iterator it_end = new_loci.end();
			UINT last_ch = chromLocusPair(*it).first;
			// create the first chromosome
			new_numLoci.push_back(0);
			new_chromTypes.push_back(chromType(last_ch));
			new_chromNames.push_back(chromName(last_ch));
			for (; it != it_end; ++it) {
				DBG_FAILIF(*it >= totNumLoci(), IndexError,
					"Locus index " + toStr(*it) + " out of range.");
				// if new chromosome
				UINT ch = chromLocusPair(*it).first;
				if (ch != last_ch) {
					new_numLoci.push_back(0);
					new_chromTypes.push_back(chromType(ch));
					new_chromNames.push_back(chromName(ch));
					last_ch = ch;
				}
				// add a locus
				++new_numLoci.back();
				new_lociPos.push_back(locusPos(*it));
				new_lociNames.push_back(locusName(*it));
				new_alleleNames.push_back(alleleNames(*it));
			}
		}
		DBG_DO(DBG_POPULATION, cerr << "Extract population with \nnumLoci:" << new_numLoci
			                        << "\nchromType: " << new_chromTypes
			                        << "\nlociPos: " << new_lociPos
			                        << "\nchromNames: " << new_chromNames
			                        << "\nlociNames: " << new_lociNames
			                        << "\ninfoFields: " << allInfoFields
			                        << endl);
		pop.setGenoStructure(ploidy(), new_numLoci, new_chromTypes, isHaplodiploid(),
			new_lociPos, new_chromNames, new_alleleNames, new_lociNames, allInfoFields);
	}
	UINT step = pop.genoSize();
	UINT infoStep = pop.infoSize();
	// will be used to copy from *this
	UINT pStep = totNumLoci();
	UINT pEnd = ploidy() * totNumLoci();
	vectoru::iterator lociPtr = new_loci.begin();
	vectoru::iterator lociEnd = new_loci.end();
	//
	UINT info = removeInd ? (ped ? ped->infoIdx(field) : infoIdx(field)) : 0;
	vectoru infoList;
	vectorstr::const_iterator iit = keptInfoFields.begin();
	vectorstr::const_iterator iit_end = keptInfoFields.end();
	for (; iit != iit_end; ++iit)
		infoList.push_back(infoIdx(*iit));
	vectoru pedInfoList;
	if (ped != NULL && !pedFields.empty()) {
		iit = pedFields.begin();
		iit_end = pedFields.end();
		for (; iit != iit_end; ++iit)
			pedInfoList.push_back(ped->infoIdx(*iit));
	}
	vectoru::iterator infoPtr = infoList.begin();
	vectoru::iterator infoEnd = infoList.end();
	vectoru::iterator pedInfoPtr = pedInfoList.begin();
	vectoru::iterator pedInfoEnd = pedInfoList.end();
	//
	// copy individuals, from ancestor to current.
	int depth = ancestralGens();
	if (ancGen >= 0 && ancGen < depth)
		depth = ancGen;
	// ancestral depth can be -1
	pop.setAncestralDepth(m_ancestralGens);
	for (; depth >= 0; --depth) {
		const_cast<population *>(this)->useAncestralGen(depth);
		sortIndividuals();
		if (ped) {
			ped->useAncestralGen(depth);
			ped->sortIndividuals();
		}
		// determine the number of individuals
		vectoru spSizes;
		vector<vectoru> indIdx;
		ULONG size;
		if (!removeInd) {
			spSizes = subPopSizes();
			size = popSize();
		} else {
			for (ULONG i = 0; i < popSize(); ++i) {
				int sp = ped ? ped->ind(i).intInfo(info) : m_inds[i].intInfo(info);
				if (sp < 0)
					continue;
				if (spSizes.size() <= static_cast<UINT>(sp)) {
					spSizes.resize(sp + 1);
					indIdx.resize(sp + 1);
				}
				++spSizes[sp];
				indIdx[sp].push_back(i);
			}
			size = accumulate(spSizes.begin(), spSizes.end(), 0UL);
			DBG_DO(DBG_POPULATION, cerr << "New subpopulation size " << spSizes << endl);
		}

		vector<individual> new_inds;
		vectora new_genotype;
		vectorinfo new_info;

		new_inds.reserve(size);
		new_genotype.reserve(size * step);
		new_info.reserve(size * infoStep);
		// copy genotype and info...
		if (!removeInd) {
			new_inds.insert(new_inds.end(), m_inds.begin(), m_inds.end());
			// handle genotype
			if (!removeLoci)
				new_genotype.insert(new_genotype.end(), m_genotype.begin(), m_genotype.end());
			else {
				ConstRawIndIterator it = rawIndBegin();
				ConstRawIndIterator it_end = rawIndEnd();
				for (; it != it_end; ++it) {
					GenoIterator ptr = it->genoBegin();
					for (UINT p = 0; p < pEnd; p += pStep) {
						for (lociPtr = new_loci.begin(); lociPtr != lociEnd; ++lociPtr)
							new_genotype.push_back(*(ptr + *lociPtr + p));
					}
				}
			}
			// handle information fields
			if (!removeInfo && pedFields.empty())
				new_info.insert(new_info.end(), m_info.begin(), m_info.end());
			else if (ped != NULL && !pedFields.empty()) {
				// copy information fields from pop as well as pedigree
				ConstRawIndIterator it = rawIndBegin();
				ConstRawIndIterator it_end = rawIndEnd();
				ConstRawIndIterator pedIt = ped->rawIndBegin();
				for (; it != it_end; ++it, ++pedIt) {
					InfoIterator iPtr = it->infoBegin();
					for (infoPtr = infoList.begin(); infoPtr != infoEnd; ++infoPtr)
						new_info.push_back(*(iPtr + *infoPtr));
					iPtr = pedIt->infoBegin();
					for (pedInfoPtr = pedInfoList.begin(); pedInfoPtr != pedInfoEnd; ++pedInfoPtr)
						new_info.push_back(*(iPtr + *pedInfoPtr));
				}
			} else {
				ConstRawIndIterator it = rawIndBegin();
				ConstRawIndIterator it_end = rawIndEnd();
				for (; it != it_end; ++it) {
					InfoIterator iPtr = it->infoBegin();
					for (infoPtr = infoList.begin(); infoPtr != infoEnd; ++infoPtr)
						new_info.push_back(*(iPtr + *infoPtr));
				}
			}
		} else {
			// remove individual
			for (size_t sp = 0; sp < indIdx.size(); ++sp) {
				vectoru & idx = indIdx[sp];
				if (idx.empty())
					continue;
				vectoru::iterator it = idx.begin();
				vectoru::iterator it_end = idx.end();
				for (; it != it_end; ++it) {
					new_inds.push_back(m_inds[*it]);
					if (!removeLoci)
						new_genotype.insert(new_genotype.end(),
							indGenoBegin(*it), indGenoEnd(*it));
					else {
						GenoIterator ptr = indGenoBegin(*it);
						for (UINT p = 0; p < pEnd; p += pStep) {
							for (lociPtr = new_loci.begin(); lociPtr != lociEnd; ++lociPtr)
								new_genotype.push_back(*(ptr + *lociPtr + p));
						}
					}
					if (!removeInfo && pedFields.empty())
						new_info.insert(new_info.end(), m_inds[*it].infoBegin(),
							m_inds[*it].infoEnd());
					else if (ped != NULL && !pedFields.empty()) {
						InfoIterator iPtr = m_inds[*it].infoBegin();
						for (infoPtr = infoList.begin(); infoPtr != infoEnd; ++infoPtr)
							new_info.push_back(*(iPtr + *infoPtr));
						iPtr = ped->m_inds[*it].infoBegin();
						for (pedInfoPtr = pedInfoList.begin(); pedInfoPtr != pedInfoEnd; ++pedInfoPtr)
							new_info.push_back(*(iPtr + *pedInfoPtr));
					} else {
						InfoIterator iPtr = m_inds[*it].infoBegin();
						for (infoPtr = infoList.begin(); infoPtr != infoEnd; ++infoPtr)
							new_info.push_back(*(iPtr + *infoPtr));
					}
				}
			}
		}
		pop.m_popSize = size;
		if (m_subPopNames.empty() || m_subPopNames.size() != spSizes.size())
			pop.setSubPopStru(spSizes, vectorstr());
		else
			pop.setSubPopStru(spSizes, m_subPopNames);
		// set pointer
		vectora::iterator ptr = new_genotype.begin();
		vectorinfo::iterator infoPtr = new_info.begin();
		for (size_t i = 0; i < size; ++i, ptr += step, infoPtr += infoStep) {
			new_inds[i].setGenoStruIdx(pop.genoStruIdx());
			new_inds[i].setGenoPtr(ptr);
			new_inds[i].setInfoPtr(infoPtr);
		}
		// the arrays are ready, are they?
		DBG_ASSERT(new_inds.size() == size && (new_genotype.size() == size * step)
			&& (new_info.size() == size * infoStep), SystemError,
			"Failed to copy genotype:"
			"\ninds: " + toStr(new_inds.size()) + ", " + toStr(size) +
			"\ngenotype: " + toStr(new_genotype.size()) + ", " + toStr(size * step) +
			"\ninfo: " + toStr(new_info.size()) + ", " + toStr(size * infoStep));
		// now put them to use
		if (depth == 0) { // current generation
			pop.m_inds.swap(new_inds);
			pop.m_genotype.swap(new_genotype);
			pop.m_info.swap(new_info);
		} else {
			pop.m_ancestralPops.push_front(popData());
			popData & pd = pop.m_ancestralPops.front();
			pd.m_subPopSize.swap(spSizes);
			pd.m_genotype.swap(new_genotype);
			pd.m_info.swap(new_info);
			pd.m_inds.swap(new_inds);
		}
	}
	pop.m_curAncestralGen = 0;
	//
	return pop;
}


void population::removeLoci(const uintList & lociList, const uintList & keepList)
{
	const vectoru & loci = lociList.elems();
	const vectoru & keep = keepList.elems();

	DBG_FAILIF(!loci.empty() && !keep.empty(), ValueError,
		"Please specify only one of parameters loci and keep");

	if (loci.empty() && keep.empty())
		return;

	vectoru kept = keep;
	// kept must be in order so that genotypes could be copied correctly
	std::sort(kept.begin(), kept.end());
	UINT oldTotNumLoci = totNumLoci();
	// new geno structure is in effective now!
	setGenoStructure(gsRemoveLoci(loci, kept));

	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		//
		vectora newGenotype(genoSize() * m_popSize);
		// copy data over
		GenoIterator newPtr = newGenotype.begin();
		UINT pEnd = ploidy();
		for (ULONG i = 0; i < m_popSize; ++i) {
			// set new geno structure
			m_inds[i].setGenoStruIdx(genoStruIdx());
			GenoIterator oldPtr = m_inds[i].genoPtr();
			// new genotype
			m_inds[i].setGenoPtr(newPtr);
			for (UINT p = 0; p < pEnd; ++p) {
				vectoru::iterator loc = kept.begin();
				for (; loc != kept.end(); ++loc)
					// this line needs ordered kept array
					*(newPtr++) = oldPtr[*loc];
				oldPtr += oldTotNumLoci;
			}
		}
		m_genotype.swap(newGenotype);
	}
	setIndOrdered(true);
}


void population::recodeAlleles(const uintListFunc & newAlleles, const uintList & loci_,
                               const stringMatrix & alleleNamesMatrix)
{
	DBG_FAILIF(newAlleles.empty() && !newAlleles.func().isValid(), ValueError,
		"Please specify new alleles or a conversion function");

	const matrixstr & alleleNames = alleleNamesMatrix.elems();

	const vectoru & loci = loci_.elems();

	DBG_FAILIF(alleleNames.size() > 1 &&
		((loci_.allAvail() && alleleNames.size() != totNumLoci()) ||
		 (!loci_.allAvail() && alleleNames.size() != loci.size())),
		ValueError,
		"If locus-specific allele names are specified, they should be specified for all loci.");

	if (!alleleNames.empty()) {
		setGenoStructure(gsSetAlleleNames(loci_, alleleNames));
		for (int depth = ancestralGens(); depth >= 0; --depth) {
			useAncestralGen(depth);
			RawIndIterator it = rawIndBegin();
			RawIndIterator itEnd = rawIndEnd();
			for (; it != itEnd; ++it)
				it->setGenoStruIdx(genoStruIdx());
		}
	}

	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);

		GenoIterator ptr = m_genotype.begin();
		GenoIterator ptrEnd = m_genotype.end();
		if (!newAlleles.empty()) {
			const vectoru & map = newAlleles.elems();
			if (loci_.allAvail()) {
				for (; ptr != ptrEnd; ++ptr) {
					DBG_FAILIF(static_cast<UINT>(*ptr) >= map.size(),
						ValueError, "Allele " + toStr(int(*ptr)) + " can not be recoded");
					*ptr = ToAllele(map[*ptr]);
				}
			} else {
				UINT numLoci = totNumLoci();
				UINT iEnd = loci.size();
				for (; ptr != ptrEnd; ptr += numLoci) {
					for (size_t i = 0; i < iEnd; ++i) {
						DBG_FAILIF(loci[i] >= numLoci, IndexError, "Loci index out of range");
						DBG_FAILIF(static_cast<UINT>(*ptr) >= map.size(),
							ValueError, "Allele " + toStr(int(*ptr)) + " can not be recoded");
						*(ptr + loci[i]) = ToAllele(map[*(ptr + loci[i])]);
					}
				}
			}
			continue;
		} else {
			pyFunc func = newAlleles.func();
			UINT numLoci = totNumLoci();
			UINT iEnd = loci.size();
			for (; ptr != ptrEnd; ptr += numLoci) {
				if (loci_.allAvail()) {
					for (size_t i = 0; i < numLoci; ++i) {
						*(ptr + i) = ToAllele(func(PyObj_As_Int, "(ii)",
								static_cast<int>(*(ptr + i)), i));
					}
				} else {
					for (size_t i = 0; i < iEnd; ++i) {
						DBG_FAILIF(loci[i] >= numLoci, IndexError, "Loci index out of range");
						*(ptr + loci[i]) = ToAllele(func(PyObj_As_Int, "(ii)",
								static_cast<int>(*(ptr + loci[i])), loci[i]));
					}
				}
			}
		}
	}
}


void population::push(population & rhs)
{
	DBG_ASSERT(rhs.genoStruIdx() == genoStruIdx(), ValueError,
		"Evolution can not continue because the new generation has different \n"
		"genotypic structure.\n");

	DBG_FAILIF(!m_genotype.empty() && m_genotype.begin() == rhs.m_genotype.begin(), ValueError,
		"Passed population is a reference of current population, swapPop failed.");

	// front -1 pop, -2 pop, .... end
	//
	if (m_ancestralGens > 0
	    && ancestralGens() == static_cast<size_t>(m_ancestralGens))
		m_ancestralPops.pop_back();

	// save current population
	if (m_ancestralGens != 0) {
		// add a empty popData
		m_ancestralPops.push_front(popData());
		// get its reference
		popData & pd = m_ancestralPops.front();
		// swap with real data
		// current population may *not* be in order
		pd.swap(*this);
	}

	// then swap out data
	// can not use population::swap because it swaps too much data
	m_popSize = rhs.m_popSize;
	m_subPopSize.swap(rhs.m_subPopSize);
	m_subPopNames.swap(rhs.m_subPopNames);
	m_subPopIndex.swap(rhs.m_subPopIndex);
	std::swap(m_vspSplitter, rhs.m_vspSplitter);
	m_genotype.swap(rhs.m_genotype);
	m_info.swap(rhs.m_info);
	m_inds.swap(rhs.m_inds);
	std::swap(m_indOrdered, rhs.m_indOrdered);
	// current population should be working well
	// (with all datamember copied form rhs
	// rhs may not be working well since m_genotype etc
	// may be from ancestral pops
	if (rhs.m_popSize != rhs.m_inds.size()) {
		// keep size if pop size is OK.
		// remove all supopulation structure of rhs
		rhs.m_popSize = rhs.m_inds.size();
		rhs.setSubPopStru(rhs.m_subPopSize, rhs.m_subPopNames);
	}
	validate("Current population after push and discard:");
	rhs.validate("Outside population after push and discard:");
}


void population::addInfoFields(const stringList & fieldList, double init)
{
	const vectorstr & fields = fieldList.elems();

	DBG_ASSERT(m_info.size() == infoSize() * popSize(), SystemError,
		"Info size is wrong");

	vectorstr newfields;

	// oldsize, this is valid for rank 0
	UINT os = infoSize();
	for (vectorstr::const_iterator it = fields.begin(); it != fields.end(); ++it) {
		try {
			// has field
			UINT idx = infoIdx(*it);
			// only needs to initialize
			int oldAncPop = m_curAncestralGen;
			for (UINT anc = 0; anc <= m_ancestralPops.size(); anc++) {
				useAncestralGen(anc);

				for (IndIterator ind = indIterator(); ind.valid(); ++ind)
					ind->setInfo(init, idx);
			}
			useAncestralGen(oldAncPop);
		} catch (IndexError &) {
			newfields.push_back(*it);
		}
	}

	// add these fields
	if (!newfields.empty()) {
		setGenoStructure(gsAddInfoFields(newfields));

		// adjust information size.
		UINT is = infoSize();
		int oldAncPop = m_curAncestralGen;
		for (UINT anc = 0; anc <= m_ancestralPops.size(); anc++) {
			useAncestralGen(anc);
			vectorinfo newInfo(is * popSize(), 0.);
			// copy the old stuff in
			InfoIterator ptr = newInfo.begin();
			for (IndIterator ind = indIterator(); ind.valid(); ++ind) {
				copy(ind->infoBegin(), ind->infoBegin() + os, ptr);
				ind->setInfoPtr(ptr);
				ind->setGenoStruIdx(genoStruIdx());
				fill(ind->infoBegin() + os, ind->infoEnd(), init);
				ptr += is;
			}
			m_info.swap(newInfo);
		}
		useAncestralGen(oldAncPop);
	}
}


void population::setInfoFields(const stringList & fieldList, double init)
{
	const vectorstr & fields = fieldList.elems();

	setGenoStructure(gsSetInfoFields(fields));
	// reset info vector
	int oldAncPop = m_curAncestralGen;
	UINT is = infoSize();
	for (UINT anc = 0; anc <= m_ancestralPops.size(); anc++) {
		useAncestralGen(anc);
		vectorinfo newInfo(is * popSize(), init);
		InfoIterator ptr = newInfo.begin();
		for (IndIterator ind = indIterator(); ind.valid(); ++ind, ptr += is) {
			ind->setInfoPtr(ptr);
			ind->setGenoStruIdx(genoStruIdx());
		}
		m_info.swap(newInfo);
	}
	useAncestralGen(oldAncPop);
}


void population::updateInfoFieldsFrom(const stringList & fieldList, const population & pop,
                                      const stringList & fromFieldList, int ancGen)
{
	const vectorstr & fields = fieldList.elems();
	const vectorstr & fromFields = fromFieldList.elems();

	int depth = ancestralGens();

	if (ancGen > 0 && ancGen < depth)
		depth = ancGen;
	for (; depth >= 0; --depth) {
		useAncestralGen(depth);
		const_cast<population &>(pop).useAncestralGen(depth);
		DBG_FAILIF(subPopSizes() != pop.subPopSizes(), ValueError,
			"Two populations should have the same population structure.");
		for (UINT i = 0; i < fields.size(); ++i) {
			UINT fromIdx = fromFields.empty() ? pop.infoIdx(fields[i]) : pop.infoIdx(fromFields[i]);
			UINT toIdx = pop.infoIdx(fields[i]);
			// indInfo is supposed to be const, but it is troublesome to change that.
			setIndInfo(const_cast<population &>(pop).indInfo(fromIdx), toIdx);
		}
	}
}


void population::setIndInfo(const floatList & valueList, const uintString & field, vspID subPop)
{
	DBG_FAILIF(subPop.valid() && hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	UINT idx = field.empty() ? field.value() : infoIdx(field.name());

	CHECKRANGEINFO(idx);
	const vectorf & values = valueList.elems();
	size_t valueSize = values.size();
	if (subPop.valid()) {
		IndInfoIterator ptr = infoBegin(idx, subPop);
		for (size_t i = 0; ptr != infoEnd(idx, subPop); ++ptr, ++i)
			*ptr = static_cast<InfoType>(values[i % valueSize]);
	} else {
		IndInfoIterator ptr = infoBegin(idx);
		for (size_t i = 0; ptr != infoEnd(idx); ++ptr, ++i)
			*ptr = static_cast<InfoType>(values[i % valueSize]);
	}
}


// set ancestral depth, can be -1
void population::setAncestralDepth(int depth)
{
	// just to make sure.
	useAncestralGen(0);
	//
	if (depth >= 0 && m_ancestralPops.size() > static_cast<size_t>(depth)) {
		int numRemove = m_ancestralPops.size() - depth;
		while (numRemove-- > 0)
			m_ancestralPops.pop_back();
	}
	DBG_ASSERT(depth < 0 || m_ancestralPops.size() <= static_cast<size_t>(depth), SystemError,
		"Failed to change ancestral Depth");

	m_ancestralGens = depth;
}


void population::useAncestralGen(UINT idx)
{
	if (m_curAncestralGen >= 0 && idx == static_cast<UINT>(m_curAncestralGen))
		return;

	DBG_DO(DBG_POPULATION, cerr << "Use ancestral generation: " << idx <<
		" Current ancestral index: " << m_curAncestralGen << endl);

	if (idx == 0 || m_curAncestralGen != 0) {         // recover pop.
		popData & pd = m_ancestralPops[ m_curAncestralGen - 1];
		pd.swap(*this);
		m_curAncestralGen = 0;
		if (idx == 0) {                                               // restore key parameters from data
			m_popSize = m_inds.size();
			setSubPopStru(m_subPopSize, m_subPopNames);
			return;
		}
	}

	// now m_curAncestralGen is zero.
	DBG_ASSERT(idx <= m_ancestralPops.size(),
		ValueError, "Ancestry generation " + toStr(idx) + " does not exist.");

	// now idx should be at least 1
	m_curAncestralGen = idx;
	// swap  1 ==> 0, 2 ==> 1

	popData & pd = m_ancestralPops[m_curAncestralGen - 1];
	pd.swap(*this);
	m_popSize = m_inds.size();
	setSubPopStru(m_subPopSize, m_subPopNames);
}


void population::save(const string & filename) const
{
	boost::iostreams::filtering_ostream ofs;

	ofs.push(boost::iostreams::gzip_compressor());
	ofs.push(boost::iostreams::file_sink(filename, std::ios::binary));

	if (!ofs)
		throw ValueError("Can not open file " + filename);

	boost::archive::text_oarchive oa(ofs);
	oa << *this;
}


void population::load(const string & filename)
{
	boost::iostreams::filtering_istream ifs;

	ifs.push(boost::iostreams::gzip_decompressor());
	ifs.push(boost::iostreams::file_source(filename, std::ios::binary));
	// do not have to test again.
	if (!ifs)
		throw ValueError("Can not open file " + filename);

	// try to load the file
	try {
		boost::archive::text_iarchive ia(ifs);
		ia >> *this;
	} catch (...) {
		throw ValueError("Failed to load population " + filename + ".\n");
	}                                                                                               // try bin
}


PyObject * population::vars(vspID vsp)
{
	if (!vsp.valid()) {
		Py_INCREF(m_vars.dict());
		return m_vars.dict();
	}

	DBG_ASSERT(static_cast<UINT>(vsp.subPop()) < numSubPop(),
		IndexError, "Subpop index out of range of 0 ~ " + toStr(numSubPop() - 1) );

	DBG_ASSERT(m_vars.hasVar("subPop"), ValueError,
		"subPop statistics does not exist yet.");

	PyObject * spObj = m_vars.getVar("subPop");
	// vsp? A tube with (sp, vsp)
	PyObject * key = NULL;
	if (vsp.isVirtual())
		key = Py_BuildValue("(ii)", vsp.subPop(), vsp.virtualSubPop());
	else
		key = PyInt_FromLong(vsp.subPop());

	spObj = PyDict_GetItem(spObj, key);

	DBG_FAILIF(spObj == NULL, ValueError,
		"Statistics for specified (virtual) subpopulation does not exist.");

	Py_INCREF(spObj);
	Py_INCREF(key);
	return spObj;
}


// The same as vars(), but without increasing
// reference count.
PyObject * population::dict(int subPop)
{
	if (subPop < 0)
		return m_vars.dict();
	else {
		DBG_ASSERT(static_cast<UINT>(subPop) < numSubPop(),
			IndexError, "Subpop index out of range of 0 ~ " + toStr(numSubPop() - 1) );

		DBG_ASSERT(m_vars.hasVar("subPop"), ValueError,
			"subPop statistics does not exist yet.");

		PyObject * spObj = m_vars.getVar("subPop");
		spObj = PyList_GetItem(spObj, subPop);

		DBG_ASSERT(spObj != NULL, SystemError,
			"Something is wrong about the length of subPop list. ");

		return spObj;
	}
}


void population::sortIndividuals(bool infoOnly) const
{
	if (indOrdered())
		return;

	if (infoOnly) {
		DBG_DO(DBG_POPULATION, cerr << "Adjust info position " << endl);
		UINT is = infoSize();
		if (is == 0) {
			setIndOrdered(true);
			return;
		}
		vectorinfo tmpInfo(m_popSize * is);
		vectorinfo::iterator infoPtr = tmpInfo.begin();

		IndIterator ind = const_cast<population *>(this)->indIterator();
		for (; ind.valid(); ++ind) {
			copy(ind->infoBegin(), ind->infoEnd(), infoPtr);
			ind->setInfoPtr(infoPtr);
			infoPtr += is;
		}
		const_cast<population *>(this)->m_info.swap(tmpInfo);
	} else {
		DBG_DO(DBG_POPULATION, cerr << "Adjust geno and info position " << endl);

		size_t sz = genoSize();
		UINT is = infoSize();
		vectora tmpGenotype(m_popSize * genoSize());
		vectorinfo tmpInfo(m_popSize * infoSize());
		vectora::iterator it = tmpGenotype.begin();
		vectorinfo::iterator infoPtr = tmpInfo.begin();

		IndIterator ind = const_cast<population *>(this)->indIterator();
		for (; ind.valid(); ++ind) {
#ifdef BINARYALLELE
			copyGenotype(ind->genoBegin(), it, sz);
#else
			copy(ind->genoBegin(), ind->genoEnd(), it);
#endif
			ind->setGenoPtr(it);
			it += sz;

			copy(ind->infoBegin(), ind->infoEnd(), infoPtr);
			ind->setInfoPtr(infoPtr);
			infoPtr += is;
		}
		// discard original genotype
		const_cast<population *>(this)->m_genotype.swap(tmpGenotype);
		const_cast<population *>(this)->m_info.swap(tmpInfo);
	}
	setIndOrdered(true);
}


population & LoadPopulation(const string & file)
{
	population * p = new population();

	p->load(file);
	return *p;
}


}


