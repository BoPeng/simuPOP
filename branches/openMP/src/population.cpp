/**
 *  $File: Population.cpp $
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

#include "population.h"
#include "pedigree.h"
#include "virtualSubPop.h"

// for file compression
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#if PY_VERSION_HEX >= 0x03000000
#  define PyInt_FromLong(x) PyLong_FromLong(x)
#endif

namespace simuPOP {

Population::Population(const uintList & size,
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
	m_curAncestralGen(0),
	m_indOrdered(true),
	m_gen(0),
	m_rep(0)
{
	DBG_DO(DBG_POPULATION, cerr << "Constructor of population is called\n");

	// get a GenoStructure with parameters. GenoStructure may be shared by some populations
	// a whole set of functions ploidy() etc in GenoStruTriat can be used after this step.
	PARAM_FAILIF(static_cast<size_t>(ploidy) * 1.0 != double(ploidy) && fcmp_ne(ploidy, HAPLODIPLOID),
		ValueError, "Only integer ploidy number or HAPLODIPLOID can be specified");

	setGenoStructure(fcmp_eq(ploidy, HAPLODIPLOID) ? 2 : static_cast<size_t>(ploidy),
		loci.elems(), chromTypes.elems(), fcmp_eq(ploidy, HAPLODIPLOID), lociPos.elems(),
		chromNames.elems(), alleleNames.elems(), lociNames.elems(), infoFields.elems());

	DBG_DO(DBG_DEVEL, cerr	<< "individual size is " << sizeof(Individual) << '+'
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


Population::~Population()
{
	DBG_DO(DBG_POPULATION,
		cerr << "Destructor of population is called" << endl);

	if (m_vspSplitter)
		delete m_vspSplitter;

	decGenoStruRef();
}


Population::Population(const Population & rhs) :
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
	m_curAncestralGen(rhs.m_curAncestralGen),
	m_indOrdered(true),
	m_gen(rhs.m_gen),
	m_rep(rhs.m_rep)
{
	DBG_DO(DBG_POPULATION,
		cerr << "Copy constructor of population is called" << endl);

	try {
		m_inds.resize(rhs.m_popSize);
		m_genotype.resize(m_popSize * genoSize());
		// have 0 length for mpi/non-head node
		m_info.resize(rhs.m_popSize * infoSize());
	} catch (...) {
		throw RuntimeError("Failed to copy Population, likely a memory allocation failure.");
	}

	// individuals will always have the correct genostructure
	// by using their copied pointer
	// Population, however, need to set this pointer correctly
	//
	setGenoStruIdx(rhs.genoStruIdx());
	incGenoStruRef();   // inc ref to the new

	// copy genotype one by one so Individual genoPtr will not
	// point outside of subpopulation region.
	GenoIterator ptr = m_genotype.begin();
	InfoIterator infoPtr = m_info.begin();
	size_t step = genoSize();
	size_t infoStep = infoSize();
	for (size_t i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		m_inds[i].setGenoPtr(ptr);
		m_inds[i].setInfoPtr(infoPtr);
		m_inds[i].copyFrom(rhs.m_inds[i]);
	}

	// copy ancestral populations
	try {
		// copy all. Individual will be shallow copied
		m_ancestralPops = rhs.m_ancestralPops;
		// need to setGenoPtr
		for (size_t ap = 0; ap < m_ancestralPops.size(); ++ap) {
			popData & lp = m_ancestralPops[ap];
			const popData & rp = rhs.m_ancestralPops[ap];

			vector<Individual> & linds = lp.m_inds;
			const vector<Individual> & rinds = rp.m_inds;

			GenoIterator lg = lp.m_genotype.begin();
			ConstGenoIterator rg = rp.m_genotype.begin();

			InfoIterator li = lp.m_info.begin();
			ConstInfoIterator ri = rp.m_info.begin();

			size_t ps = rinds.size();

			for (size_t i = 0; i < ps; ++i) {
				linds[i].setGenoPtr(rinds[i].genoPtr() - rg + lg);
				linds[i].setInfoPtr(rinds[i].infoPtr() - ri + li);
			}
		}
	} catch (...) {
		cerr	<< "Unable to copy ancestral populations. "
		        << "The popolation size may be too big." << endl
		        << "The population will still be usable but without any ancestral population stored." << endl;
		m_ancestralGens = 0;
		m_ancestralPops.clear();
	}

	// copy virtual subpop splitters
	setVirtualSplitter(rhs.virtualSplitter());
}


void Population::popData::swap(Population & pop)
{
	pop.m_subPopSize.swap(m_subPopSize);
	pop.m_subPopNames.swap(m_subPopNames);
	pop.m_genotype.swap(m_genotype);
	pop.m_info.swap(m_info);
	pop.m_inds.swap(m_inds);
	std::swap(pop.m_indOrdered, m_indOrdered);
}


Population * Population::clone() const
{
	return new Population(*this);
}


size_t Population::subPopByName(const string & name) const
{
	vectorstr::const_iterator it = find(m_subPopNames.begin(), m_subPopNames.end(), name);

	if (it == m_subPopNames.end())
		throw ValueError("SubPopulation " + name + " not found.");
	return it - m_subPopNames.begin();
}


string Population::subPopName(vspID subPop) const
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


vectorstr Population::subPopNames() const
{
	DBG_ASSERT(m_subPopNames.empty() || m_subPopNames.size() == numSubPop(), SystemError,
		"subpopulation names can either be empty, or be specified for all subpopulations.");
	return m_subPopNames.empty() ? vectorstr(numSubPop(), UnnamedSubPop) : m_subPopNames;
}


void Population::setSubPopName(const string & name, size_t subPop)
{
	CHECKRANGESUBPOP(subPop);
	if (m_subPopNames.empty())
		m_subPopNames = vectorstr(numSubPop(), UnnamedSubPop);
	m_subPopNames[subPop] = name;
}


bool Population::hasActivatedVirtualSubPop() const
{
	return m_vspSplitter != NULL && m_vspSplitter->activatedSubPop() != InvalidValue;
}


bool Population::hasActivatedVirtualSubPop(size_t subPop) const
{
	return m_vspSplitter != NULL && m_vspSplitter->activatedSubPop() == subPop;
}


bool Population::hasVirtualSubPop() const
{
	return m_vspSplitter != NULL;
}


void Population::setVirtualSplitter(BaseVspSplitter * vsp)
{
	if (m_vspSplitter)
		delete m_vspSplitter;

	m_vspSplitter = vsp ? vsp->clone() : NULL;
}


size_t Population::numVirtualSubPop() const
{
	return hasVirtualSubPop()
	       ? m_vspSplitter->numVirtualSubPop()
		   : 0;
}


void Population::activateVirtualSubPop(vspID subPop) const
{
	CHECKRANGESUBPOP(subPop.subPop());
	if (!subPop.isVirtual())
		return;
	DBG_ASSERT(hasVirtualSubPop(), ValueError, "population has no virtual subpopulations");
	m_vspSplitter->activate(*this, subPop.subPop(), subPop.virtualSubPop());
	DBG_ASSERT(m_vspSplitter->activatedSubPop() == subPop.subPop(), SystemError,
		"Failed to activate virtual subpopulation");
}


void Population::deactivateVirtualSubPop(size_t subPop) const
{
	CHECKRANGESUBPOP(subPop);
	if (!hasActivatedVirtualSubPop(subPop))
		return;
	m_vspSplitter->deactivate(subPop);
}


int Population::__cmp__(const Population & rhs) const
{
	if (genoStruIdx() != rhs.genoStruIdx()) {
		DBG_DO(DBG_POPULATION, cerr << "Genotype structures are different" << endl);
		return 1;
	}

	if (popSize() != rhs.popSize()) {
		DBG_DO(DBG_POPULATION, cerr << "population sizes are different" << endl);
		return 1;
	}

	if (ancestralGens() != rhs.ancestralGens()) {
		DBG_DO(DBG_POPULATION, cerr << "Number of ancestral generations differ" << endl);
		return 1;
	}

	int curGen = m_curAncestralGen;
	int rhsCurGen = rhs.m_curAncestralGen;
	for (int depth = ancestralGens(); depth >= 0; --depth) {
		const_cast<Population *>(this)->useAncestralGen(depth);
		const_cast<Population &>(rhs).useAncestralGen(depth);
		for (size_t i = 0, iEnd = popSize(); i < iEnd; ++i)
			if (m_inds[i] != rhs.m_inds[i]) {
				DBG_DO(DBG_POPULATION, cerr << "individuals are different" << endl);
				const_cast<Population *>(this)->useAncestralGen(curGen);
				const_cast<Population &>(rhs).useAncestralGen(rhsCurGen);
				return 1;
			}
	}
	const_cast<Population *>(this)->useAncestralGen(curGen);
	const_cast<Population &>(rhs).useAncestralGen(rhsCurGen);

	return 0;
}


Individual & Population::indByID(double fid, const uintList & ancGens, const string & idField)
{
	size_t id = toID(fid);

	DBG_FAILIF(fabs(fid - id) > 1e-8, ValueError,
		"individual ID has to be integer (or a double round to full iteger).");

	size_t idx = infoIdx(idField);

	vectoru gens = ancGens.elems();
	if (ancGens.allAvail())
		for (int gen = 0; gen <= ancestralGens(); ++gen)
			gens.push_back(gen);
	else if (ancGens.unspecified())
		gens.push_back(m_curAncestralGen);

	for (size_t genIdx = 0; genIdx < gens.size(); ++genIdx) {
		ssize_t gen = gens[genIdx];
		vector<Individual> * inds = NULL;
		// search in current, not necessarily the present generation
		if (gen == m_curAncestralGen)
			inds = &m_inds;
		else
			inds = &m_ancestralPops[gen == 0 ? m_curAncestralGen - 1 : gen - 1].m_inds;
		// first try our luck
		size_t startID = (*inds)[0].intInfo(idx);
		if (idx >= startID && startID + (*inds).size() > id) {
			Individual & ind = (*inds)[id - startID];
			if (toID(ind.intInfo(idx)) == id)
				return ind;
		}
		// now we have to search all individuals
		for (size_t i = 0; i < (*inds).size(); ++i) {
			if (toID((*inds)[i].intInfo(idx)) == id)
				return (*inds)[i];
		}
	}
	// if still cannot be found, raise an IndexError.
	throw IndexError("No individual with ID " + toStr(id) + " could be found.");
	// this is just to suppress a warning.
	return m_inds[0];
}


pyIndIterator Population::individuals(vspID subPopID)
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), RuntimeError,
		"Can not call individuals when there is activated virtual subpopulation");

	vspID subPop = subPopID.resolve(*this);

	DBG_FAILIF(subPop.allAvailSP() || subPop.allAvailVSP(), ValueError, "Invalid (virtual) subpopulation ID.")
	if (!subPop.valid())
		return pyIndIterator(m_inds.begin(), m_inds.end(), true, vspFunctor());

	size_t spID = subPop.subPop();

#ifndef OPTIMIZED
	CHECKRANGESUBPOP(spID);
	size_t vspID = subPop.virtualSubPop();
	CHECKRANGEVIRTUALSUBPOP(vspID);
#endif
	if (subPop.isVirtual())
		return pyIndIterator(m_inds.begin() + subPopBegin(spID),
			m_inds.begin() + subPopEnd(spID), false,
			vspFunctor(*this, m_vspSplitter, subPop));
	else
		return pyIndIterator(m_inds.begin() + subPopBegin(spID),
			m_inds.begin() + subPopEnd(spID), true, vspFunctor());
}


Individual & Population::ancestor(double fidx, ssize_t gen, vspID vsp)
{
	size_t idx = toID(fidx);

	DBG_FAILIF(fabs(fidx - idx) > 1e-8, ValueError,
		"individual index has to be integer (or a double round to full iteger).");

	DBG_FAILIF(vsp.isVirtual(), ValueError,
		"Function genotype currently does not support virtual subpopulation");

	DBG_FAILIF(static_cast<size_t>(gen) > m_ancestralPops.size(), IndexError,
		"Ancestray generation " + toStr(gen) + " does not exist");
	if (!vsp.valid()) {
		if (gen == m_curAncestralGen)
			return m_inds[idx];
		ssize_t genIdx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
		DBG_FAILIF(idx > m_ancestralPops[genIdx].m_inds.size(),
			IndexError, "individual index out of range");
		return m_ancestralPops[genIdx].m_inds[idx];
	} else {
		size_t subPop = vsp.subPop();
		if (gen == m_curAncestralGen)
			return m_inds[idx + subPopBegin(subPop)];
		ssize_t genIdx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
		DBG_FAILIF(static_cast<size_t>(subPop) > m_ancestralPops[genIdx].m_subPopSize.size(),
			IndexError, "subpopulation index out of range");
		DBG_FAILIF(idx > m_ancestralPops[genIdx].m_subPopSize[subPop],
			IndexError, "individual index out of range");
		size_t shift = 0;
		if (subPop > 0) {
			for (size_t i = 0; i < subPop; ++i)
				shift += m_ancestralPops[genIdx].m_subPopSize[i];
		}
		return m_ancestralPops[genIdx].m_inds[shift + idx];
	}
}


const Individual & Population::ancestor(double fidx, ssize_t gen, vspID vsp) const
{
	size_t idx = toID(fidx);

	DBG_FAILIF(fabs(fidx - idx) > 1e-8, ValueError,
		"individual index has to be integer (or a double round to full iteger).");
	DBG_FAILIF(vsp.isVirtual(), ValueError,
		"Function genotype currently does not support virtual subpopulation");

	DBG_FAILIF(static_cast<size_t>(gen) > m_ancestralPops.size(), IndexError,
		"Ancestray generation " + toStr(gen) + " does not exist");
	if (!vsp.valid()) {
		if (gen == m_curAncestralGen)
			return m_inds[idx];
		ssize_t genIdx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
		DBG_FAILIF(idx > m_ancestralPops[genIdx].m_inds.size(),
			IndexError, "individual index out of range");
		return m_ancestralPops[genIdx].m_inds[idx];
	} else {
		size_t subPop = vsp.subPop();
		if (gen == m_curAncestralGen)
			return m_inds[idx + subPopBegin(subPop)];
		ssize_t genIdx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
		DBG_FAILIF(static_cast<size_t>(subPop) > m_ancestralPops[genIdx].m_subPopSize.size(),
			IndexError, "subpopulation index out of range");
		DBG_FAILIF(idx > m_ancestralPops[genIdx].m_subPopSize[subPop],
			IndexError, "individual index out of range");
		size_t shift = 0;
		if (subPop > 0) {
			for (size_t i = 0; i < subPop; ++i)
				shift += m_ancestralPops[genIdx].m_subPopSize[i];
		}
		return m_ancestralPops[genIdx].m_inds[shift + idx];
	}
}


IndAlleleIterator Population::alleleIterator(size_t locus)
{
	CHECKRANGEABSLOCUS(locus);

	// if there is virtual subpop, use Individual based iterator
	// or
	// if requires order, but the alleles are not ordered
	// use Individual based
	size_t ct = chromType(chromLocusPair(locus).first);
	if (hasActivatedVirtualSubPop() || !indOrdered()
	    || (ct != AUTOSOME && ct != CUSTOMIZED) || isHaplodiploid())
		// this is a complex case
		return IndAlleleIterator(locus, indIterator());
	else
		// a simplere case with straightforward iterator
		return IndAlleleIterator(locus, m_genotype.begin(), m_genotype.end(),
			totNumLoci());
}


/// CPPONLY allele begin, for given subPop
IndAlleleIterator Population::alleleIterator(size_t locus, size_t subPop)
{
	CHECKRANGEABSLOCUS(locus);
	CHECKRANGESUBPOP(subPop);

	size_t ct = chromType(chromLocusPair(locus).first);
	// this is a complex case
	if (hasActivatedVirtualSubPop() || !indOrdered()
	    || (ct != AUTOSOME && ct != CUSTOMIZED) || isHaplodiploid())
		// this is a complex case
		return IndAlleleIterator(locus, indIterator(subPop));
	else
		// this is a complex case
		return IndAlleleIterator(locus,
			m_genotype.begin() + m_subPopIndex[subPop] * genoSize(),
			m_genotype.begin() + m_subPopIndex[subPop + 1] * genoSize(),
			totNumLoci());
}


PyObject * Population::genotype(vspID subPopID)
{

	DBG_WARNIF(true, "The returned object of function Population.genotype() is a special "
		             "carray object that reflects the underlying genotype of a "
		             "population. It will become invalid once the population changes. "
		             "Please use list(pop.genotype()) if you would like to keep a copy of genotypes");

	vspID vsp = subPopID.resolve(*this);

	DBG_FAILIF(vsp.isVirtual(), ValueError,
		"Function genotype currently does not support virtual subpopulation");
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	syncIndPointers();
	if (!vsp.valid()) {
		// directly expose values. Do not copy data over.
		return Allele_Vec_As_NumArray(m_genotype.begin(), m_genotype.end());
	} else {
		size_t subPop = vsp.subPop();
		CHECKRANGESUBPOP(subPop);
		// directly expose values. Do not copy data over.
		return Allele_Vec_As_NumArray(genoBegin(subPop, true), genoEnd(subPop, true));
	}
	return NULL;
}


void Population::setGenotype(const uintList & genoList, vspID subPopID)
{
	const vectoru & geno = genoList.elems();

	vspID subPop = subPopID.resolve(*this);

	syncIndPointers();
	if (!subPop.valid()) {
		GenoIterator ptr = m_genotype.begin();
		size_t sz = geno.size();
		for (size_t i = 0; i < popSize() * genoSize(); ++i)
			*(ptr++) = ToAllele(geno[i % sz]);
		return;
	}

	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	size_t sp = subPop.subPop();
	CHECKRANGESUBPOP(sp);

	size_t sz = geno.size();
	if (!subPop.isVirtual()) {
		GenoIterator ptr = genoBegin(sp, true);
		for (size_t i = 0; i < subPopSize(sp) * genoSize(); ++i)
			*(ptr++) = ToAllele(geno[i % sz]);
	} else {
		activateVirtualSubPop(subPop);
		IndIterator it = indIterator(sp);
		size_t i = 0;
		for (; it.valid(); ++it)
			for (GenoIterator git = it->genoBegin(); git != it->genoEnd(); ++git, ++i)
				*git = ToAllele(geno[i % sz]);
		deactivateVirtualSubPop(subPop.subPop());
	}
}


void Population::validate(const string & msg) const
{
#ifdef OPTIMIZED
	(void)msg;  // avoid a warning message of unused parameter in optmized module
#else
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


void Population::fitSubPopStru(const vectoru & newSubPopSizes,
                               const vectorstr & newSubPopNames)
{
	size_t newSize = accumulate(newSubPopSizes.begin(), newSubPopSizes.end(), size_t(0));

	bool needsResize = m_popSize != newSize;

	if (needsResize) {
		size_t is = infoSize();
		size_t step = genoSize();
		m_popSize = newSize;
		try {
			if (step != 0 && m_popSize > MaxIndexSize / step)
				throw RuntimeError("Population size times number of loci exceed maximum index size.");
			m_genotype.resize(m_popSize * step);
			m_info.resize(m_popSize * is);
			m_inds.resize(m_popSize);
		} catch (...) {
			throw RuntimeError("Failed to create population (popSize=" + toStr(m_popSize)
				+ ", totNumLoci*ploidy=" + toStr(step) + ", maximum population size for such a long genome=" +
				toStr(MaxIndexSize / step) + ", requested memory="
#ifdef BINARYALLELE
				+ toStr((m_popSize * step / 8 + m_popSize * is * sizeof(double) + m_popSize * sizeof(Individual)) / 1024)
#else
				+ toStr((m_popSize * step * sizeof(Allele) + m_popSize * is * sizeof(double) + m_popSize * sizeof(Individual)) / 1024)
#endif
				+ "k bytes)");
		}
		// reset individual pointers
		GenoIterator ptr = m_genotype.begin();
		InfoIterator infoPtr = m_info.begin();
		for (size_t i = 0; i < m_popSize; ++i, ptr += step, infoPtr += is) {
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


void Population::fitGenoStru(size_t stru)
{
	// set genotypic structure to a population.
	// This function will try not to change population size.
	size_t oldSize = genoSize();
	size_t oldInfoSize = infoSize();

	decGenoStruRef();   // dec ref to the old
	setGenoStruIdx(stru);
	incGenoStruRef();   // inc ref to the new
	size_t newSize = genoSize();
	size_t newInfoSize = infoSize();

	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		if (oldSize != newSize)
			m_genotype.resize(newSize * popSize());
		if (oldInfoSize != newInfoSize)
			m_info.resize(newInfoSize * popSize());
		// reset structure
		GenoIterator ptr = m_genotype.begin();
		InfoIterator infoPtr = m_info.begin();
		for (size_t i = 0; i < m_popSize; ++i, ptr += newSize, infoPtr += newInfoSize) {
			m_inds[i].setGenoStruIdx(stru);
			m_inds[i].setGenoPtr(ptr);
			m_inds[i].setInfoPtr(infoPtr);
		}
	}
}


void Population::setSubPopStru(const vectoru & newSubPopSizes,
                               const vectorstr & newSubPopNames)
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	DBG_ASSERT(accumulate(newSubPopSizes.begin(), newSubPopSizes.end(), size_t(0)) == m_popSize, ValueError,
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
	size_t i = 1;
	for (m_subPopIndex[0] = 0; i <= numSubPop(); ++i)
		m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];
}


size_t Population::subPopSize(vspID subPopID, int ancGen) const
{
	vspID subPop = subPopID.resolve(*this);

	if (!subPop.valid())
		return popSize(ancGen);

	DBG_FAILIF(ancGen > ancestralGens(),
		IndexError, "Ancestral generation " + toStr(ancGen) + " out of range of 0 ~ "
		+ toStr(ancestralGens()));

	if (ancGen < 0 || ancGen == m_curAncestralGen) {
		CHECKRANGESUBPOP(subPop.subPop());
		CHECKRANGEVIRTUALSUBPOP(subPop.virtualSubPop());
		if (subPop.isVirtual())
			return m_vspSplitter->size(*this, subPop.subPop(), subPop.virtualSubPop());
		else
			return m_subPopSize[subPop.subPop()];
	} else if (subPop.isVirtual()) {
		int curGen = m_curAncestralGen;
		const_cast<Population *>(this)->useAncestralGen(ancGen);
		CHECKRANGESUBPOP(subPop.subPop());
		CHECKRANGEVIRTUALSUBPOP(subPop.virtualSubPop());
		size_t size = m_vspSplitter->size(*this, subPop.subPop(), subPop.virtualSubPop());
		const_cast<Population *>(this)->useAncestralGen(curGen);
		return size;
	} else {
		const vectoru & sizes = m_ancestralPops[ancGen - 1].m_subPopSize;
		DBG_FAILIF(static_cast<size_t>(subPop.subPop()) >= sizes.size(), IndexError,
			"Subpopulation index " + toStr(subPop.subPop()) + " out of range of 0 ~ "
			+ toStr(sizes.size() - 1) + " at ancestral generation " + toStr(ancGen));
		return sizes[subPop.subPop()];
	}
	// avoid a warning message.
	return 0;
}


void Population::sortIndividuals(const stringList & infoList)
{
	const vectorstr & infoFields = infoList.elems();

	if (infoFields.size() == 0)
		return;
	vectoru fields(infoFields.size());
	for (size_t i = 0; i < infoFields.size(); ++i)
		fields[i] = infoIdx(infoFields[i]);
	for (size_t sp = 0; sp < numSubPop(); ++sp)
		std::sort(rawIndBegin(sp), rawIndEnd(sp), indCompare(fields));
	setIndOrdered(false);
}


void Population::setSubPopByIndInfo(const string & field)
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	size_t info = infoIdx(field);

	DBG_DO(DBG_POPULATION, cerr << "Sorting individuals." << endl);
	// sort individuals first
	std::sort(rawIndBegin(), rawIndEnd(), indCompare(info));
	setIndOrdered(false);

	// sort individuals first
	// remove individuals with negative index.
	if (indIterator()->info(info) < 0) {
		// popsize etc will be changed.
		size_t newPopSize = m_popSize;
		IndIterator it = indIterator();
		for (; it.valid(); ++it) {
			if (it->info(info) < 0)
				newPopSize-- ;
			else
				break;
		}
		// 'it' now point to the one with positive info(info)

		DBG_DO(DBG_POPULATION, cerr << "New pop size" << newPopSize << endl);

		// allocate new genotype and inds
		vectora newGenotype(genoSize() * newPopSize);
		vectorf newInfo(newPopSize * infoSize());
		vector<Individual> newInds(newPopSize);

		// assign genotype location and set structure information for individuals
		GenoIterator ptr = newGenotype.begin();
		InfoIterator infoPtr = newInfo.begin();
		size_t step = genoSize();
		size_t infoStep = infoSize();
		for (size_t i = 0; i < newPopSize; ++i, ptr += step, ++it, infoPtr += infoStep) {
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
		size_t numSubPop = static_cast<size_t>(m_inds.back().info(info)) + 1;
		m_subPopSize.resize(numSubPop);
		m_subPopIndex.resize(numSubPop + 1);

		// check subpop size
		fill(m_subPopSize.begin(), m_subPopSize.end(), 0);
		for (IndIterator it = indIterator(); it.valid(); ++it)
			m_subPopSize[ static_cast<size_t>(it->info(info)) ]++;
	}
	// rebuild index
	size_t i = 1;
	for (m_subPopIndex[0] = 0; i <= numSubPop(); ++i)
		m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];
	// subpopulation names
	if (!m_subPopNames.empty())
		m_subPopNames.resize(numSubPop(), UnnamedSubPop);
}


vectoru Population::splitSubPop(size_t subPop, const vectorf & sizes, const vectorstr & names)
{
	if (sizes.size() <= 1)
		return vectoru(1, subPop);

	double s = accumulate(sizes.begin(), sizes.end(), 0.);
	vectoru count(sizes.size());

	if (fcmp_eq(s, 1.)) {
		// proportions
		propToCount(sizes.begin(), sizes.end(), subPopSize(subPop), count);
	} else if (fcmp_eq(s, static_cast<double>(subPopSize(subPop)))) {
		for (size_t i = 0; i < sizes.size(); ++i) {
			count[i] = static_cast<size_t>(sizes[i] + 0.5);
		}
	}
	DBG_FAILIF(accumulate(count.begin(), count.end(), size_t(0)) != subPopSize(subPop), ValueError,
		"Sum of parameter sizes should be 1 or the size of subpopulation " + toStr(subPop));

	DBG_ASSERT(names.empty() || sizes.size() == names.size(), ValueError,
		"Names should be given to none or all of the split subpopulations");

	if (!names.empty() && m_subPopNames.empty())
		m_subPopNames.resize(numSubPop(), UnnamedSubPop);

	vectoru subPopSizes;
	vectorstr subPopNames;
	vectoru ret(count.size());
	for (size_t sp = 0; sp < numSubPop(); ++sp) {
		if (sp != subPop) {
			subPopSizes.push_back(subPopSize(sp));
			if (!m_subPopNames.empty())
				subPopNames.push_back(m_subPopNames[sp]);
		} else {
			subPopSizes.insert(subPopSizes.end(), count.begin(), count.end());
			if (!m_subPopNames.empty()) {
				if (names.empty()) {
					for (size_t i = 0; i < count.size(); ++i)
						subPopNames.push_back(m_subPopNames[subPop]);
				} else
					subPopNames.insert(subPopNames.end(), names.begin(), names.end());
			}
			for (size_t i = 0; i < count.size(); ++i)
				ret[i] = sp + i;
		}
	}
	setSubPopStru(subPopSizes, subPopNames);
	return ret;
}


void Population::removeSubPops(const subPopList & subPops)
{
	syncIndPointers();
	vectoru new_size;
	vectorstr new_spNames;

	size_t step = genoSize();
	size_t infoStep = infoSize();
	RawIndIterator oldInd = m_inds.begin();
	RawIndIterator newInd = m_inds.begin();
	GenoIterator oldPtr = m_genotype.begin();
	GenoIterator newPtr = m_genotype.begin();
	InfoIterator oldInfoPtr = m_info.begin();
	InfoIterator newInfoPtr = m_info.begin();

#ifndef OPTIMIZED
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	for (; it != itEnd; ++it) {
		CHECKRANGESUBPOP(it->subPop());
		CHECKRANGEVIRTUALSUBPOP(it->virtualSubPop());
	}
#endif
	for (size_t sp = 0; sp < numSubPop(); ++sp) {
		size_t spSize = subPopSize(sp);
		if (subPops.contains(sp)) {
			// complete removal, move oldPtr but not newPtr
			oldInd += spSize;
			oldPtr += step * spSize;
			oldInfoPtr += infoStep * spSize;
		} else if (subPops.overlap(sp)) {
			// partial removal
			//
			// mark for removal
			markIndividuals(sp, false);
			subPopList::const_iterator it = subPops.begin();
			subPopList::const_iterator itEnd = subPops.end();
			for (; it != itEnd; ++it)
				if (it->subPop() == static_cast<size_t>(sp))
					markIndividuals(*it, true);
			// remove individuals, but not subpopulation
			size_t newSize = 0;
			for (size_t i = 0; i < spSize; ++i) {
				// will be kept
				if (!oldInd->marked()) {
					++newSize;
					if (oldInd != newInd) {
						*newInd = *oldInd;
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
			//
			new_size.push_back(newSize);
			if (!m_subPopNames.empty())
				new_spNames.push_back(m_subPopNames[sp]);
		} else {
			// do not remove
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
			oldInd += spSize;
			oldPtr += step * spSize;
			oldInfoPtr += infoStep * spSize;
		}
	}
	//
	m_inds.erase(newInd, m_inds.end());
	m_genotype.erase(newPtr, m_genotype.end());
	m_info.erase(newInfoPtr, m_info.end());
	m_popSize = std::accumulate(new_size.begin(), new_size.end(), size_t(0));
	setSubPopStru(new_size, new_spNames);
	//
	GenoIterator ptr = m_genotype.begin();
	InfoIterator infoPtr = m_info.begin();
	for (size_t i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		m_inds[i].setGenoPtr(ptr);
		m_inds[i].setInfoPtr(infoPtr);
	}
}


void Population::removeMarkedIndividuals()
{
	syncIndPointers();
	vectoru new_size(numSubPop(), 0);

	size_t step = genoSize();
	size_t infoStep = infoSize();
	RawIndIterator oldInd = m_inds.begin();
	RawIndIterator newInd = m_inds.begin();
	GenoIterator oldPtr = m_genotype.begin();
	InfoIterator oldInfoPtr = m_info.begin();
	GenoIterator newPtr = m_genotype.begin();
	InfoIterator newInfoPtr = m_info.begin();
	//
	for (size_t sp = 0; sp < numSubPop(); ++sp) {
		size_t newSize = 0;
		size_t spSize = subPopSize(sp);
		for (size_t i = 0; i < spSize; ++i) {
			// will be kept
			if (!oldInd->marked()) {
				++newSize;
				if (oldInd != newInd) {
					*newInd = *oldInd;
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
		new_size[sp] = newSize;
	}
	//
	m_inds.erase(newInd, m_inds.end());
	m_genotype.erase(newPtr, m_genotype.end());
	m_info.erase(newInfoPtr, m_info.end());
	m_popSize = std::accumulate(new_size.begin(), new_size.end(), size_t(0));
	setSubPopStru(new_size, m_subPopNames);
	//
	GenoIterator ptr = m_genotype.begin();
	InfoIterator infoPtr = m_info.begin();
	for (size_t i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		m_inds[i].setGenoPtr(ptr);
		m_inds[i].setInfoPtr(infoPtr);
	}
}


void Population::removeIndividuals(const uintList & indexList, const floatList & IDList,
                                   const string & idField, PyObject * filter)
{
	const vectoru & indexes = indexList.elems();
	const vectorf & IDs = IDList.elems();

	if (IDs.empty() && indexes.empty() && filter == NULL)
		return;

	DBG_FAILIF(IDs.empty() + indexes.empty() + (filter == NULL) != 2, ValueError,
		"Please specify only one of parameters indexes, IDs and filter");

	if (!indexes.empty()) {
		markIndividuals(vspID(), false);
		for (size_t i = 0; i < indexes.size(); ++i) {
			DBG_FAILIF(indexes[i] >= m_popSize, IndexError,
				"individual index out of range.");
			m_inds[indexes[i]].setMarked(true);
		}
		removeMarkedIndividuals();
		return;
	}
	int curGen = m_curAncestralGen;
	if (!IDs.empty()) {
		size_t fieldIdx = infoIdx(idField);
		// remove by ID
		// first, build a map
		std::map<size_t, bool > idMap;
		for (size_t i = 0; i < IDs.size(); ++i)
			idMap[toID(IDs[i])] = true;

		for (int depth = ancestralGens(); depth >= 0; --depth) {
			useAncestralGen(depth);
			markIndividuals(vspID(), false);
			RawIndIterator it = rawIndBegin();
			RawIndIterator itEnd = rawIndEnd();
			for (; it != itEnd; ++it) {
				size_t id = toID(it->info(fieldIdx));
				if (idMap.find(id) != idMap.end())
					it->setMarked(true);
			}
			removeMarkedIndividuals();
		}
	}
	if (filter != NULL) {
		pyFunc func(filter);
		PyObject * args = PyTuple_New(func.numArgs());
		//
		vectoru pars(func.numArgs());
		for (size_t i = 0; i < func.numArgs(); ++i) {
			const string & arg = func.arg(i);
			if (arg == "ind")
				pars[i] = InvalidValue;
			else {
				DBG_FAILIF(!hasInfoField(arg), ValueError,
					"Only parameters 'ind', and names of information fields are "
					"acceptable in the filter function");
				pars[i] = infoIdx(arg);
			}
		}

		for (int depth = ancestralGens(); depth >= 0; --depth) {
			useAncestralGen(depth);
			markIndividuals(vspID(), false);
			RawIndIterator it = rawIndBegin();
			RawIndIterator itEnd = rawIndEnd();
			for (; it != itEnd; ++it) {
				for (size_t i = 0; i < func.numArgs(); ++i) {
					if (pars[i] == InvalidValue)
						PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(&*it)));
					else
						PyTuple_SET_ITEM(args, i, PyFloat_FromDouble(it->info(pars[i])));
				}
				if (func(PyObj_As_Bool, args))
					it->setMarked(true);
			}
			removeMarkedIndividuals();
		}
		Py_XDECREF(args);
	}
	useAncestralGen(curGen);
}


size_t Population::mergeSubPops(const uintList & subPops, const string & name)
{
	if (!name.empty() && m_subPopNames.empty())
		m_subPopNames.resize(numSubPop(), UnnamedSubPop);

	// merge all subpopulations
	if (subPops.allAvail()) {
		// [ popSize() ]
		vectoru sz(1, popSize());
		if (m_subPopNames.empty())
			setSubPopStru(sz, m_subPopNames);
		else
			setSubPopStru(sz, vectorstr(1, name.empty() ? m_subPopNames[0] : name));
		return 0;
	}
	if (subPops.elems().size() == 1) {
		if (!name.empty())
			m_subPopNames[0] = name;
		return subPops.elems()[0];
	}

	// are they in order?
	bool consecutive = true;
	vectoru sps = subPops.elems();
	sort(sps.begin(), sps.end());
	for (size_t i = 1; i < sps.size(); ++i)
		if (sps[i] != sps[i - 1] + 1) {
			consecutive = false;
			break;
		}
	// new subpopulation sizes and names
	vectoru new_size;
	vectorstr new_names;
	for (size_t sp = 0; sp < numSubPop(); ++sp) {
		if (find(sps.begin(), sps.end(), sp) != sps.end()) {
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
	syncIndPointers();
	// find the new subpop order
	vectoru sp_order;
	for (size_t sp = 0; sp < sps[0]; ++sp)
		sp_order.push_back(sp);
	sp_order.insert(sp_order.end(), sps.begin(), sps.end());
	for (size_t sp = sps[0]; sp < numSubPop(); ++sp)
		if (find(sps.begin(), sps.end(), sp) == sps.end())
			sp_order.push_back(sp);
	//
	DBG_ASSERT(sp_order.size() == numSubPop(), ValueError,
		"Incorrect resulting subpopulation number, maybe caused by duplicate entries in parameter subPops.");

	size_t step = genoSize();
	size_t infoStep = infoSize();
	vector<Individual> new_inds;
	vectora new_genotype;
	vectorf new_info;
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
	for (size_t i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		m_inds[i].setGenoPtr(ptr);
		m_inds[i].setInfoPtr(infoPtr);
	}
	return sps[0];
}


void Population::addChromFrom(const Population & pop)
{
	size_t numLoci1 = totNumLoci();
	size_t numLoci2 = pop.totNumLoci();

	// obtain new genotype structure and set it
	setGenoStructure(gsAddChromFromStru(pop.genoStruIdx()));
	//
	DBG_FAILIF(ancestralGens() != pop.ancestralGens(), ValueError,
		"Can not add chromosomes from a population with different number of ancestral generations");
	//
	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		const_cast<Population &>(pop).useAncestralGen(depth);
		//
		DBG_FAILIF(m_subPopSize != pop.m_subPopSize, ValueError,
			"Can not add chromosomes from a population with different subpopulation sizes");

		vectora newGenotype(genoSize() * m_popSize);

		// append pop2 chromosomes to the first one
		GenoIterator ptr = newGenotype.begin();
		size_t pEnd = ploidy();
		for (size_t i = 0; i < m_popSize; ++i) {
			// set new geno structure
			m_inds[i].setGenoStruIdx(genoStruIdx());
			GenoIterator ptr1 = m_inds[i].genoPtr();
			GenoIterator ptr2 = pop.m_inds[i].genoPtr();
			m_inds[i].setGenoPtr(ptr);
			for (size_t p = 0; p < pEnd; ++p) {
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
		syncIndPointers(true);
}


void Population::addIndFrom(const Population & pop)
{
	DBG_FAILIF(genoStruIdx() != pop.genoStruIdx(), ValueError,
		"Cannot add Individual from a population with different genotypic structure.");
	DBG_FAILIF(ancestralGens() != pop.ancestralGens(), ValueError,
		"Two populations should have the same number of ancestral generations.");
	// genotype pointers may be reset so this is needed.
	syncIndPointers();
	const_cast<Population &>(pop).syncIndPointers();
	// go to the oldest generation
	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		const_cast<Population &>(pop).useAncestralGen(depth);
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
		size_t step = genoSize();
		size_t infoStep = infoSize();
		// set pointers
		for (size_t i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
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


void Population::addLociFrom(const Population & pop)
{
	DBG_FAILIF(ancestralGens() != pop.ancestralGens(), ValueError,
		"Can not add chromosomes from a population with different number of ancestral generations");

	size_t size1 = totNumLoci();
	size_t size2 = pop.totNumLoci();
	// obtain new genotype structure and set it
	vectoru indexes1;
	vectoru indexes2;
	setGenoStructure(gsAddLociFromStru(pop.genoStruIdx(),
			indexes1, indexes2));

	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		const_cast<Population &>(pop).useAncestralGen(depth);
		//
		DBG_FAILIF(m_subPopSize != pop.m_subPopSize, ValueError,
			"Can not add chromosomes from a population with different subpopulation sizes");
		//
		vectora newGenotype(genoSize() * m_popSize);

		// merge chromosome by chromosome
		GenoIterator ptr = newGenotype.begin();
		size_t pEnd = ploidy();
		size_t newSize = totNumLoci();
		for (size_t i = 0; i < m_popSize; ++i) {
			// set new geno structure
			m_inds[i].setGenoStruIdx(genoStruIdx());
			GenoIterator ptr1 = m_inds[i].genoPtr();
			GenoIterator ptr2 = pop.m_inds[i].genoPtr();
			// new genotype
			m_inds[i].setGenoPtr(ptr);
			// copy each allele
			for (size_t p = 0; p < pEnd; ++p) {
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
	syncIndPointers(true);
}


void Population::addChrom(const floatList & lociPosList, const stringList & lociNameList,
                          const string & chromName, const stringMatrix & alleleNames,
                          size_t chromType)
{
	const vectorf & lociPos = lociPosList.elems();
	const vectorstr & lociNames = lociNameList.elems();

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
		size_t newPopGenoSize = genoSize() * m_popSize;
		vectora newGenotype(newPopGenoSize, 0);

		// copy data over
		GenoIterator newPtr = newGenotype.begin();
		size_t pEnd = ploidy();
		size_t gap = totNumLoci() - oldNumLoci;
		for (size_t i = 0; i < m_popSize; ++i) {
			// set new geno structure
			m_inds[i].setGenoStruIdx(genoStruIdx());
			GenoIterator oldPtr = m_inds[i].genoPtr();
			// new genotype
			m_inds[i].setGenoPtr(newPtr);
			// copy each chromosome
			for (size_t p = 0; p < pEnd; ++p) {
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
	syncIndPointers(true);
}


vectoru Population::addLoci(const uintList & chromList, const floatList & posList,
                            const stringList & lociNameList, const stringMatrix & alleleNamesMatrix)
{
	const vectoru & chrom = chromList.elems();
	const vectorf & pos = posList.elems();
	const matrixstr & alleleNames = alleleNamesMatrix.elems();
	const vectorstr & lociNames = lociNameList.elems();

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
		size_t newPopGenoSize = genoSize() * m_popSize;
		vectora newGenotype(newPopGenoSize, 0);

		// copy data over
		GenoIterator newPtr = newGenotype.begin();
		size_t pEnd = ploidy();
		for (size_t i = 0; i < m_popSize; ++i) {
			// set new geno structure
			m_inds[i].setGenoStruIdx(genoStruIdx());
			GenoIterator oldPtr = m_inds[i].genoPtr();
			// new genotype
			m_inds[i].setGenoPtr(newPtr);
			// copy each chromosome
			for (size_t p = 0; p < pEnd; ++p) {
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
	syncIndPointers(true);
	return newIndex;
}


void Population::resize(const uintList & sizeList, bool propagate)
{
	const vectoru & newSubPopSizes = sizeList.elems();

	DBG_FAILIF(newSubPopSizes.size() != numSubPop(), ValueError,
		"Resize should give subpopulation size for each subpopulation");

	size_t newPopSize = accumulate(newSubPopSizes.begin(), newSubPopSizes.end(), size_t(0));

	// prepare new Population
	vector<Individual> newInds(newPopSize);
	vectora newGenotype(genoSize() * newPopSize);
	vectorf newInfo(newPopSize * infoSize());
	// iterators ready
	GenoIterator ptr = newGenotype.begin();
	InfoIterator infoPtr = newInfo.begin();
	size_t step = genoSize();
	size_t infoStep = infoSize();
	// set pointers
	for (size_t i = 0; i < newPopSize; ++i, ptr += step, infoPtr += infoStep) {
		newInds[i].setGenoStruIdx(genoStruIdx());
		newInds[i].setGenoPtr(ptr);
		newInds[i].setInfoPtr(infoPtr);
	}
	// copy stuff over
	size_t startSP = 0;
	for (size_t sp = 0; sp < numSubPop(); ++sp) {
		size_t spSize = subPopSize(sp);
		for (size_t i = 0, j = 0; i < newSubPopSizes[sp]; ++j, ++i) {
			// repeating?
			if (spSize == 0 || ((j / spSize) > 0 && !propagate))
				break;
			newInds[startSP + i].copyFrom(m_inds[j % spSize + subPopBegin(sp)]);
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


Population & Population::extractSubPops(const subPopList & subPops, bool rearrange) const
{
#ifndef OPTIMIZED
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	for (; it != itEnd; ++it) {
		CHECKRANGESUBPOP(it->subPop());
		CHECKRANGEVIRTUALSUBPOP(it->virtualSubPop());
	}
#endif
	Population & pop = *new Population();

	pop.setGenoStruIdx(genoStruIdx());
	incGenoStruRef();
	pop.setVirtualSplitter(virtualSplitter());

	syncIndPointers();
	vectoru new_size;
	vectorstr new_spNames;

	size_t step = genoSize();
	size_t infoStep = infoSize();

	vector<Individual> new_inds;
	vectora new_genotype;
	vectorf new_info;

	if (rearrange) {
		size_t sz = 0;
		subPopList::const_iterator it = subPops.begin();
		subPopList::const_iterator itEnd = subPops.end();
		for (; it != itEnd; ++it) {
			new_size.push_back(subPopSize(*it));
			sz += new_size.back();
			if (!m_subPopNames.empty())
				new_spNames.push_back(m_subPopNames[it->subPop()]);
		}
		//
		new_inds.resize(sz);
		new_genotype.resize(sz * step);
		new_info.resize(sz * infoStep);
		//
		RawIndIterator newInd = new_inds.begin();
		GenoIterator newPtr = new_genotype.begin();
		InfoIterator newInfoPtr = new_info.begin();

		it = subPops.begin();
		for (; it != itEnd; ++it) {
			activateVirtualSubPop(*it);
			ConstIndIterator oldInd = indIterator(it->subPop());
			for (; oldInd.valid(); ++oldInd) {
				*newInd = *oldInd;
				copy(oldInd->genoBegin(), oldInd->genoEnd(), newPtr);
				copy(oldInd->infoBegin(), oldInd->infoEnd(), newInfoPtr);
				++newInd;
				newPtr += step;
				newInfoPtr += infoStep;
			}
			deactivateVirtualSubPop(it->subPop());
		}
	} else {
		size_t sz = 0;
		// there is overestimate here
		for (size_t sp = 0; sp < numSubPop(); ++sp) {
			if (subPops.contains(sp) || subPops.overlap(sp))
				sz += subPopSize(sp);
		}
		new_inds.resize(sz);
		new_genotype.resize(sz * step);
		new_info.resize(sz * infoStep);
		//
		RawIndIterator newInd = new_inds.begin();
		GenoIterator newPtr = new_genotype.begin();
		InfoIterator newInfoPtr = new_info.begin();

		//
		ConstRawIndIterator oldInd = m_inds.begin();
		ConstGenoIterator oldPtr = m_genotype.begin();
		ConstInfoIterator oldInfoPtr = m_info.begin();
		for (size_t sp = 0; sp < numSubPop(); ++sp) {
			size_t spSize = subPopSize(sp);
			if (subPops.contains(sp)) {
				// complete copy
				new_size.push_back(spSize);
				if (!m_subPopNames.empty())
					new_spNames.push_back(m_subPopNames[sp]);
				//
				copy(oldInd, oldInd + spSize, newInd);
				copy(oldPtr, oldPtr + step * spSize, newPtr);
				copy(oldInfoPtr, oldInfoPtr + infoStep * spSize, newInfoPtr);

				oldInd += spSize;
				oldPtr += step * spSize;
				oldInfoPtr += infoStep * spSize;
				newInd += spSize;
				newPtr += step * spSize;
				newInfoPtr += infoStep * spSize;
			} else if (subPops.overlap(sp)) {
				// partial copy
				//
				// mark for copy
				markIndividuals(sp, false);
				subPopList::const_iterator it = subPops.begin();
				subPopList::const_iterator itEnd = subPops.end();
				for (; it != itEnd; ++it)
					if (it->subPop() == static_cast<size_t>(sp))
						markIndividuals(*it, true);
				//
				size_t newSize = 0;
				for (size_t i = 0; i < spSize; ++i) {
					// will be kept
					if (oldInd->marked()) {
						++newSize;
						*newInd = *oldInd;
						copy(oldPtr, oldPtr + step, newPtr);
						copy(oldInfoPtr, oldInfoPtr + infoStep, newInfoPtr);
						++newInd;
						newPtr += step;
						newInfoPtr += infoStep;
					}
					++oldInd;
					oldPtr += step;
					oldInfoPtr += infoStep;
				}
				//
				new_size.push_back(newSize);
				if (!m_subPopNames.empty())
					new_spNames.push_back(m_subPopNames[sp]);
			} else {
				// do not copy
				oldInd += spSize;
				oldPtr += step * spSize;
				oldInfoPtr += infoStep * spSize;
			}
		}
	}
	size_t sz = std::accumulate(new_size.begin(), new_size.end(), size_t(0));
	new_inds.resize(sz);
	new_genotype.resize(sz * step);
	new_info.resize(sz * infoStep);
	//
	pop.m_inds.swap(new_inds);
	pop.m_genotype.swap(new_genotype);
	pop.m_info.swap(new_info);
	pop.m_popSize = sz;
	pop.setSubPopStru(new_size, new_spNames);
	//
	GenoIterator ptr = pop.m_genotype.begin();
	InfoIterator infoPtr = pop.m_info.begin();
	for (size_t i = 0; i < pop.m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		pop.m_inds[i].setGenoPtr(ptr);
		pop.m_inds[i].setInfoPtr(infoPtr);
	}

	return pop;
}


Population & Population::extractMarkedIndividuals() const
{
	Population & pop = *new Population();

	pop.setGenoStruIdx(genoStruIdx());
	incGenoStruRef();
	pop.setVirtualSplitter(virtualSplitter());

	syncIndPointers();
	vectoru new_size;

	size_t step = genoSize();
	size_t infoStep = infoSize();
	ConstRawIndIterator oldInd = m_inds.begin();
	ConstGenoIterator oldPtr = m_genotype.begin();
	ConstInfoIterator oldInfoPtr = m_info.begin();

	size_t sz = 0;
	ConstRawIndIterator it = rawIndBegin();
	ConstRawIndIterator itEnd = rawIndEnd();
	for (; it != itEnd; ++it)
		if (it->marked())
			++sz;

	vector<Individual> new_inds(sz);
	vectora new_genotype(sz * step);
	vectorf new_info(sz * infoStep);

	RawIndIterator newInd = new_inds.begin();
	GenoIterator newPtr = new_genotype.begin();
	InfoIterator newInfoPtr = new_info.begin();
	//
	for (size_t sp = 0; sp < numSubPop(); ++sp) {
		size_t newSize = 0;
		size_t spSize = subPopSize(sp);
		for (size_t i = 0; i < spSize; ++i) {
			// will be kept
			if (oldInd->marked()) {
				++newSize;
				*newInd = *oldInd;
				copy(oldPtr, oldPtr + step, newPtr);
				copy(oldInfoPtr, oldInfoPtr + infoStep, newInfoPtr);
				++newInd;
				newPtr += step;
				newInfoPtr += infoStep;
			}
			++oldInd;
			oldPtr += step;
			oldInfoPtr += infoStep;
		}
		new_size.push_back(newSize);
	}
	//
	pop.m_inds.swap(new_inds);
	pop.m_genotype.swap(new_genotype);
	pop.m_info.swap(new_info);
	pop.m_popSize = std::accumulate(new_size.begin(), new_size.end(), size_t(0));
	pop.setSubPopStru(new_size, m_subPopNames);
	//
	GenoIterator ptr = pop.m_genotype.begin();
	InfoIterator infoPtr = pop.m_info.begin();
	for (size_t i = 0; i < pop.m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		pop.m_inds[i].setGenoPtr(ptr);
		pop.m_inds[i].setInfoPtr(infoPtr);
	}
	return pop;
}


Population & Population::extractIndividuals(const uintList & indexList,
                                            const floatList & IDList, const string & idField,
                                            PyObject * filter) const
{
	const vectoru & indexes = indexList.elems();
	const vectorf & IDs = IDList.elems();

	if (IDs.empty() && indexes.empty() && filter == NULL) {
		// extract no individuals
		Population & pop = *new Population();
		pop.setGenoStruIdx(genoStruIdx());
		incGenoStruRef();
		pop.setVirtualSplitter(virtualSplitter());
		pop.setSubPopStru(vectoru(numSubPop(), 0), m_subPopNames);
		return pop;
	}

	DBG_FAILIF(IDs.empty() + indexes.empty() + (filter == NULL) != 2, ValueError,
		"Please specify only one of parameters indexes, IDs and filter");

	if (!indexes.empty()) {
		markIndividuals(vspID(), false);
		for (size_t i = 0; i < indexes.size(); ++i) {
			DBG_FAILIF(indexes[i] >= m_popSize, IndexError,
				"individual index " + toStr(indexes[i]) + " out of range of 0 ~ " + toStr(m_popSize) + ".");
			m_inds[indexes[i]].setMarked(true);
		}
		return extractMarkedIndividuals();
	}
	int curGen = m_curAncestralGen;
	if (!IDs.empty()) {
		size_t fieldIdx = infoIdx(idField);
		// remove by ID
		// first, build a map
		std::map<size_t, bool > idMap;
		for (size_t i = 0; i < IDs.size(); ++i)
			idMap[toID(IDs[i])] = true;

		for (int depth = ancestralGens(); depth >= 0; --depth) {
			const_cast<Population *>(this)->useAncestralGen(depth);
			markIndividuals(vspID(), false);
			ConstRawIndIterator it = rawIndBegin();
			ConstRawIndIterator itEnd = rawIndEnd();
			for (; it != itEnd; ++it) {
				size_t id = toID(it->info(fieldIdx));
				if (idMap.find(id) != idMap.end())
					it->setMarked(true);
			}
		}
	}
	if (filter != NULL) {
		pyFunc func(filter);
		PyObject * args = PyTuple_New(func.numArgs());
		vectori pars(func.numArgs());
		for (size_t i = 0; i < func.numArgs(); ++i) {
			const string & arg = func.arg(i);
			if (arg == "ind")
				pars[i] = -1;
			else {
				DBG_FAILIF(!hasInfoField(arg), ValueError,
					"Only parameters 'ind', and names of information fields are "
					"acceptable in the filter function");
				pars[i] = static_cast<int>(infoIdx(arg));
			}
		}

		for (int depth = ancestralGens(); depth >= 0; --depth) {
			const_cast<Population *>(this)->useAncestralGen(depth);
			markIndividuals(vspID(), false);
			ConstRawIndIterator it = rawIndBegin();
			ConstRawIndIterator itEnd = rawIndEnd();
			for (; it != itEnd; ++it) {
				for (size_t i = 0; i < func.numArgs(); ++i) {
					if (pars[i] < 0)
						PyTuple_SET_ITEM(args, i, pyIndObj(
								static_cast<void *>(const_cast<Individual *>(&*it))));
					else
						PyTuple_SET_ITEM(args, i, PyFloat_FromDouble(it->info(pars[i])));
				}
				if (func(PyObj_As_Bool, args))
					it->setMarked(true);
			}
		}
		Py_XDECREF(args);
	}
	Population * allPop = NULL;
	for (int depth = ancestralGens(); depth >= 0; --depth) {
		const_cast<Population *>(this)->useAncestralGen(depth);
		Population & pop = extractMarkedIndividuals();
		if (allPop == NULL) {
			allPop = &pop;
			allPop->setAncestralDepth(ancestralGens());
		} else
			allPop->push(pop);
	}
	const_cast<Population *>(this)->useAncestralGen(curGen);
	return *allPop;
}


Population & Population::extract(const lociList & extractedLoci, const stringList & infoFieldList,
                                 const subPopList & _subPops, const uintList & ancGens) const
{
	Population & pop = *new Population();

	// the usual whole population, easy case.
	subPopList subPops = _subPops.expandFrom(*this);

	bool removeInd = !_subPops.allAvail();
	bool removeLoci = !extractedLoci.allAvail();
	const vectoru & loci = extractedLoci.elems(this);
	bool removeInfo = !infoFieldList.allAvail();
	const vectorstr & infoFields = infoFieldList.elems();

	DBG_DO(DBG_POPULATION, cerr << "Remove ind: " << removeInd
		                        << "\nRemove loci: " << removeLoci
		                        << "\nRemove info: " << removeInfo << endl);

	// will keep a sorted version of loci
	vectoru new_loci = loci;

	vectorstr keptInfoFields = removeLoci ? infoFields : this->infoFields();

	// Population strcture.
	if (!removeLoci && !removeInfo) {
		pop.setGenoStruIdx(genoStruIdx());
		incGenoStruRef();
	} else if (!removeLoci) {
		// only change information fields
		pop.setGenoStructure(ploidy(), numLoci(), chromTypes(), isHaplodiploid(),
			lociPos(), chromNames(), allAlleleNames(), lociNames(), keptInfoFields);
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
			size_t last_ch = chromLocusPair(*it).first;
			// create the first chromosome
			new_numLoci.push_back(0);
			new_chromTypes.push_back(chromType(last_ch));
			new_chromNames.push_back(chromName(last_ch));
			for (; it != it_end; ++it) {
				DBG_FAILIF(*it >= totNumLoci(), IndexError,
					"Locus index " + toStr(*it) + " out of range.");
				// if new chromosome
				size_t ch = chromLocusPair(*it).first;
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
			                        << "\ninfoFields: " << keptInfoFields
			                        << endl);
		pop.setGenoStructure(ploidy(), new_numLoci, new_chromTypes, isHaplodiploid(),
			new_lociPos, new_chromNames, new_alleleNames, new_lociNames, keptInfoFields);
	}
	size_t step = pop.genoSize();
	size_t infoStep = pop.infoSize();
	// will be used to copy from *this
	size_t pStep = totNumLoci();
	size_t pEnd = ploidy() * totNumLoci();
	vectoru::iterator lociPtr = new_loci.begin();
	vectoru::iterator lociEnd = new_loci.end();
	//
	vectoru infoList;
	vectorstr::const_iterator iit = keptInfoFields.begin();
	vectorstr::const_iterator iit_end = keptInfoFields.end();
	for (; iit != iit_end; ++iit)
		infoList.push_back(infoIdx(*iit));
	vectoru::iterator infoPtr = infoList.begin();
	vectoru::iterator infoEnd = infoList.end();
	//
	vectoru gens = ancGens.elems();
	if (ancGens.allAvail())
		for (int gen = 0; gen <= ancestralGens(); ++gen)
			gens.push_back(gen);
	else if (ancGens.unspecified())
		gens.push_back(m_curAncestralGen);
	std::sort(gens.begin(), gens.end());

	// ancestral depth can be -1
	pop.setAncestralDepth(m_ancestralGens);
	for (ssize_t genIdx = gens.size() - 1; genIdx >= 0; --genIdx) {
		ssize_t depth = gens[genIdx];
		const_cast<Population *>(this)->useAncestralGen(depth);
		syncIndPointers();
		// determine the number of individuals
		vectoru spSizes(numSubPop());
		vector<vectoru> indIdx(numSubPop());
		size_t size;
		if (!removeInd) {
			spSizes = subPopSizes();
			size = popSize();
		} else {
			// mark individuals for extraction
			// unmark all
			markIndividuals(vspID(), false);
			for (size_t i = 0; i < subPops.size(); ++i)
				markIndividuals(subPops[i], true);
			// select Individuals
			for (size_t sp = 0; sp < numSubPop(); ++sp) {
				size_t spBegin = subPopBegin(sp);
				for (size_t i = 0; i < subPopSize(sp); ++i) {
					if (!m_inds[i + subPopBegin(sp)].marked())
						continue;
					++spSizes[sp];
					indIdx[sp].push_back(spBegin + i);
				}
			}
			size = accumulate(spSizes.begin(), spSizes.end(), size_t(0));
			DBG_DO(DBG_POPULATION, cerr << "New subpopulation size " << spSizes << endl);
		}

		vector<Individual> new_inds;
		vectora new_genotype;
		vectorf new_info;

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
					for (size_t p = 0; p < pEnd; p += pStep) {
						for (lociPtr = new_loci.begin(); lociPtr != lociEnd; ++lociPtr)
							new_genotype.push_back(*(ptr + *lociPtr + p));
					}
				}
			}
			// handle information fields
			if (!removeInfo)
				new_info.insert(new_info.end(), m_info.begin(), m_info.end());
			else {
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
						for (size_t p = 0; p < pEnd; p += pStep) {
							for (lociPtr = new_loci.begin(); lociPtr != lociEnd; ++lociPtr)
								new_genotype.push_back(*(ptr + *lociPtr + p));
						}
					}
					if (!removeInfo)
						new_info.insert(new_info.end(), m_inds[*it].infoBegin(),
							m_inds[*it].infoEnd());
					else {
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
		vectorf::iterator infoPtr = new_info.begin();
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
		if (genIdx == 0) { // current generation
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
	// remove empty subpopulations
	subPopList emptySubPops;
	for (size_t i = 0; i < pop.numSubPop(); ++i)
		if (pop.subPopSize(i) == 0)
			emptySubPops.push_back(i);
	if (!emptySubPops.empty())
		pop.removeSubPops(emptySubPops);
	return pop;
}


void Population::removeLoci(const lociList & removeList, const lociList & keepList)
{
	if (removeList.unspecified() && keepList.unspecified())
		return;

	DBG_FAILIF(removeList.unspecified() + keepList.unspecified() != 1, ValueError,
		"Please specify only one of parameters loci and keep");

	vectoru kept = keepList.elems(this);
	if (!keepList.unspecified()) {
		// use keep list
		if (keepList.allAvail())
			for (size_t i = 0; i < totNumLoci(); ++i)
				kept.push_back(i);
		// kept must be in order so that genotypes could be copied correctly
		std::sort(kept.begin(), kept.end());
	} else {
		// use remove list
		if (removeList.allAvail())
			kept.clear();
		else {
			vectoru loci = removeList.elems(this);
			for (size_t loc = 0; loc < totNumLoci(); ++loc) {
				if (find(loci.begin(), loci.end(), loc) == loci.end())
					kept.push_back(loc);
			}
		}
	}
	setGenoStructure(gsRemoveLoci(kept));

	size_t oldTotNumLoci = totNumLoci();

	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);
		//
		vectora newGenotype(genoSize() * m_popSize);
		// copy data over
		GenoIterator newPtr = newGenotype.begin();
		size_t pEnd = ploidy();
		for (size_t i = 0; i < m_popSize; ++i) {
			// set new geno structure
			m_inds[i].setGenoStruIdx(genoStruIdx());
			GenoIterator oldPtr = m_inds[i].genoPtr();
			// new genotype
			m_inds[i].setGenoPtr(newPtr);
			for (size_t p = 0; p < pEnd; ++p) {
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


void Population::recodeAlleles(const uintListFunc & newAlleles, const lociList & loci_,
                               const stringMatrix & alleleNamesMatrix)
{
	DBG_FAILIF(newAlleles.empty() && !newAlleles.func().isValid(), ValueError,
		"Please specify new alleles or a conversion function");

	const matrixstr & alleleNames = alleleNamesMatrix.elems();

	const vectoru & loci = loci_.elems(this);
	if (!loci_.allAvail() && loci.empty())
		return;

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

	//
	// unordered map does not appear to work here, perhaps std::pair<Allele, size_t>
	// is not hashable.
	//
	//#if TR1_SUPPORT == 0
	typedef std::map<Allele, Allele> AlleleMap;
	typedef std::map<std::pair<Allele, size_t>, Allele> AlleleLocusMap;
	//#else
	//typedef std::tr1::unordered_map<Allele, Allele> AlleleMap;
	//typedef std::tr1::unordered_map<std::pair<Allele, size_t>, Allele> AlleleLocusMap;
	//#endif

	AlleleMap alleleMap;
	AlleleLocusMap alleleLocusMap;

	size_t oldGen = curAncestralGen();
	for (int depth = ancestralGens(); depth >= 0; --depth) {
		useAncestralGen(depth);

		GenoIterator ptr = m_genotype.begin();
		GenoIterator ptrEnd = m_genotype.end();
		if (!newAlleles.empty()) {
			const vectoru & map = newAlleles.elems();
			if (loci_.allAvail()) {
				for (; ptr != ptrEnd; ++ptr) {
					if (static_cast<size_t>(*ptr) >= map.size()) {
						DBG_WARNIF(true,
							"Allele " + toStr(static_cast<size_t>(*ptr)) + " can not be recoded");
						continue;
					}
					*ptr = ToAllele(map[*ptr]);
				}
			} else {
				size_t numLoci = totNumLoci();
				size_t iEnd = loci.size();
				for (; ptr != ptrEnd; ptr += numLoci) {
					for (size_t i = 0; i < iEnd; ++i) {
						DBG_FAILIF(loci[i] >= numLoci, IndexError, "Loci index out of range");
						GenoIterator allele = ptr + loci[i];
						DBG_FAILIF(static_cast<size_t>(*allele) >= map.size(),
							ValueError, "Allele " + toStr(static_cast<size_t>(*allele)) + " can not be recoded");
						*allele = ToAllele(map[*allele]);
					}
				}
			}
		} else {
			pyFunc func = newAlleles.func();
			PyObject * args = PyTuple_New(func.numArgs());
			size_t alleleIndex = InvalidValue;
			size_t locusIndex = InvalidValue;
			for (size_t i = 0; i < func.numArgs(); ++i) {
				const string & arg = func.arg(i);
				if (arg == "allele")
					alleleIndex = i;
				else if (arg == "locus")
					locusIndex = i;
				else {
					DBG_FAILIF(true, ValueError,
						"Only parameters 'allele' and 'locus' are acceptable in a user-provided recode function.");
				}
			}
			size_t numLoci = totNumLoci();
			size_t iEnd = loci.size();
			for (; ptr != ptrEnd; ptr += numLoci) {
				if (loci_.allAvail()) {
					for (size_t i = 0; i < numLoci; ++i) {
						if (alleleIndex != InvalidValue)
							PyTuple_SET_ITEM(args, alleleIndex, PyInt_FromLong(static_cast<int>(*(ptr + i))));
						if (locusIndex != InvalidValue)
							PyTuple_SET_ITEM(args, locusIndex, PyInt_FromLong(static_cast<int>(i)));
						if (locusIndex != InvalidValue) {
							std::pair<Allele, size_t> key(*(ptr + i), i);
							AlleleLocusMap::iterator it = alleleLocusMap.find(key);
							if (it != alleleLocusMap.end())
								*(ptr + i) = it->second;
							else {
								*(ptr + i) = ToAllele(func(PyObj_As_Int, args));
								alleleLocusMap[key] = *(ptr + i);
							}
						} else {
							AlleleMap::iterator it = alleleMap.find(*(ptr + i));
							if (it != alleleMap.end())
								*(ptr + i) = it->second;
							else {
								Allele oldAllele = *(ptr + i);
								*(ptr + i) = ToAllele(func(PyObj_As_Int, args));
								alleleMap[oldAllele] = *(ptr + i);
							}
						}
					}
				} else {
					for (size_t i = 0; i < iEnd; ++i) {
						DBG_FAILIF(loci[i] >= numLoci, IndexError, "Loci index out of range");
						if (alleleIndex != InvalidValue)
							PyTuple_SET_ITEM(args, alleleIndex, PyInt_FromLong(static_cast<int>(*(ptr + loci[i]))));
						if (locusIndex != InvalidValue)
							PyTuple_SET_ITEM(args, locusIndex, PyInt_FromLong(static_cast<int>(loci[i])));
						if (locusIndex != InvalidValue) {
							std::pair<Allele, size_t> key(*(ptr + loci[i]), loci[i]);
							AlleleLocusMap::iterator it = alleleLocusMap.find(key);
							if (it != alleleLocusMap.end())
								*(ptr + loci[i]) = it->second;
							else {
								*(ptr + loci[i]) = ToAllele(func(PyObj_As_Int, args));
								alleleLocusMap[key] = *(ptr + loci[i]);
							}
						} else {
							AlleleMap::iterator it = alleleMap.find(*(ptr + loci[i]));
							if (it != alleleMap.end())
								*(ptr + loci[i]) = it->second;
							else {
								Allele oldAllele = *(ptr + loci[i]);
								*(ptr + loci[i]) = ToAllele(func(PyObj_As_Int, args));
								alleleMap[oldAllele] = *(ptr + loci[i]);
							}
						}
					}
				}
			}
			Py_DECREF(args);
		}
	}
	useAncestralGen(oldGen);
}


void Population::push(Population & rhs)
{
	DBG_ASSERT(rhs.genoStruIdx() == genoStruIdx(), ValueError,
		"Evolution can not continue because the new generation has different \n"
		"genotypic structure.\n");

	DBG_FAILIF(this == &rhs, ValueError,
		"Passed population is a reference of current population, population.push failed.");

	// front -1 pop, -2 pop, .... end
	//
	if (m_ancestralGens > 0
	    && ancestralGens() == m_ancestralGens)
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
	// can not use Population::swap because it swaps too much data
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
	rhs.validate("Outside Population after push and discard:");
}


void Population::addInfoFields(const stringList & fieldList, double init)
{
	const vectorstr & fields = fieldList.elems();

	DBG_ASSERT(m_info.size() == infoSize() * popSize(), SystemError,
		"Info size is wrong");

	vectorstr newfields;

	// oldsize, this is valid for rank 0
	size_t os = infoSize();
	for (vectorstr::const_iterator it = fields.begin(); it != fields.end(); ++it) {
		try {
			// has field
			size_t idx = infoIdx(*it);
			// only needs to initialize
			int oldAncPop = m_curAncestralGen;
			for (size_t anc = 0; anc <= m_ancestralPops.size(); anc++) {
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
		size_t is = infoSize();
		int oldAncPop = m_curAncestralGen;
		for (size_t anc = 0; anc <= m_ancestralPops.size(); anc++) {
			useAncestralGen(anc);
			vectorf newInfo(is * popSize(), 0.);
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


void Population::setInfoFields(const stringList & fieldList, double init)
{
	const vectorstr & fields = fieldList.elems();

	setGenoStructure(gsSetInfoFields(fields));
	// reset info vector
	int oldAncPop = m_curAncestralGen;
	size_t is = infoSize();
	for (size_t anc = 0; anc <= m_ancestralPops.size(); anc++) {
		useAncestralGen(anc);
		vectorf newInfo(is * popSize(), init);
		InfoIterator ptr = newInfo.begin();
		for (IndIterator ind = indIterator(); ind.valid(); ++ind, ptr += is) {
			ind->setInfoPtr(ptr);
			ind->setGenoStruIdx(genoStruIdx());
		}
		m_info.swap(newInfo);
	}
	useAncestralGen(oldAncPop);
}


void Population::removeInfoFields(const stringList & fieldList)
{
	const vectorstr & fields = fieldList.elems();

	if (fields.size() == 0)
		return;

	DBG_ASSERT(m_info.size() == infoSize() * popSize(), SystemError,
		"Info size is wrong");

	vectorstr newfields;
	vectoru oldIdx;
	for (size_t idx = 0; idx < infoSize(); ++idx) {
		string field = infoField(idx);
		if (find(fields.begin(), fields.end(), field) == fields.end()) {
			oldIdx.push_back(idx);
			newfields.push_back(field);
		}
	}

	setGenoStructure(gsSetInfoFields(newfields));

	int oldAncPop = m_curAncestralGen;
	size_t sz = infoSize();
	for (size_t anc = 0; anc <= m_ancestralPops.size(); anc++) {
		useAncestralGen(anc);
		vectorf newInfo(sz * popSize(), 0.);
		// copy the old stuff in
		InfoIterator ptr = newInfo.begin();

		for (IndIterator ind = indIterator(); ind.valid(); ++ind) {
			InfoIterator oldptr = ind->infoPtr();
			ind->setInfoPtr(ptr);
			ind->setGenoStruIdx(genoStruIdx());
			for (size_t i = 0; i < sz; ++i)
				*(ptr++) = *(oldptr + oldIdx[i]);
		}
		m_info.swap(newInfo);
	}
	useAncestralGen(oldAncPop);
}


void Population::updateInfoFieldsFrom(const stringList & fieldList, const Population & pop,
                                      const stringList & fromFieldList, const uintList & ancGens)
{
	const vectorstr & fields = fieldList.elems();
	const vectorstr & fromFields = fromFieldList.elems();

	int oldGen = m_curAncestralGen;
	vectoru gens = ancGens.elems();

	if (ancGens.allAvail())
		for (int gen = 0; gen <= ancestralGens(); ++gen)
			gens.push_back(gen);
	else if (ancGens.unspecified())
		gens.push_back(m_curAncestralGen);
	//
	for (size_t genIdx = 0; genIdx < gens.size(); ++genIdx) {
		ssize_t depth = gens[genIdx];
		useAncestralGen(depth);
		const_cast<Population &>(pop).useAncestralGen(depth);
		DBG_FAILIF(subPopSizes() != pop.subPopSizes(), ValueError,
			"Two populations should have the same population structure.");
		for (size_t i = 0; i < fields.size(); ++i) {
			size_t fromIdx = fromFields.empty() ? pop.infoIdx(fields[i]) : pop.infoIdx(fromFields[i]);
			size_t toIdx = pop.infoIdx(fields[i]);
			// indInfo is supposed to be const, but it is troublesome to change that.
			setIndInfo(const_cast<Population &>(pop).indInfo(fromIdx), toIdx);
		}
	}
	useAncestralGen(oldGen);
}


void Population::setIndInfo(const floatList & valueList, const uintString &
                            field, vspID subPopID)
{
	vspID subPop = subPopID.resolve(*this);

	DBG_FAILIF(subPop.valid() && hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	size_t idx = field.empty() ? field.value() : infoIdx(field.name());

	CHECKRANGEINFO(idx);
	const vectorf & values = valueList.elems();
	size_t valueSize = values.size();
	if (subPop.valid()) {
		activateVirtualSubPop(subPop);
		IndInfoIterator ptr = infoBegin(idx, subPop);
		for (size_t i = 0; ptr != infoEnd(idx, subPop); ++ptr, ++i)
			*ptr = static_cast<double>(values[i % valueSize]);
		deactivateVirtualSubPop(subPop.subPop());
	} else {
		IndInfoIterator ptr = infoBegin(idx);
		for (size_t i = 0; ptr != infoEnd(idx); ++ptr, ++i)
			*ptr = static_cast<double>(values[i % valueSize]);
	}
}


void Population::markIndividuals(vspID subPop, bool mark) const
{
	if (subPop.valid()) {
		activateVirtualSubPop(subPop);
		ConstIndIterator it = indIterator(subPop.subPop());
		for (; it.valid(); ++it)
			it->setMarked(mark);
		deactivateVirtualSubPop(subPop.subPop());
	} else {
		ConstIndIterator it = indIterator();
		for (; it.valid(); ++it)
			it->setMarked(mark);
	}
}


// set ancestral depth, can be -1
void Population::setAncestralDepth(int depth)
{
	// just to make sure.
	useAncestralGen(0);
	//
	if (depth >= 0 && m_ancestralPops.size() > static_cast<size_t>(depth)) {
		ssize_t numRemove = m_ancestralPops.size() - depth;
		while (numRemove-- > 0)
			m_ancestralPops.pop_back();
	}
	DBG_ASSERT(depth < 0 || m_ancestralPops.size() <= static_cast<size_t>(depth), SystemError,
		"Failed to change ancestral Depth");

	m_ancestralGens = depth;
}


void Population::keepAncestralGens(const uintList & ancGens)
{
	if (ancGens.allAvail())
		return;

	useAncestralGen(0);
	vectoru gens = ancGens.elems();
	std::sort(gens.begin(), gens.end());
	for (size_t genIdx = 0; genIdx < gens.size(); ++genIdx) {
		size_t depth = gens[genIdx];
		if (genIdx == 0) {
			// move to current generation
			if (depth != 0) { // move that gen to current
				popData & pd = m_ancestralPops[depth - 1];
				pd.swap(*this);
				m_popSize = m_inds.size();
				setSubPopStru(m_subPopSize, m_subPopNames);
			}
		} else {
			// switch around in popData
			if (depth != genIdx) {
				// depth is the existing place
				// genIdx is the new location
				popData & pd = m_ancestralPops[depth - 1];
				popData & pd1 = m_ancestralPops[genIdx - 1];
				pd1.m_subPopSize.swap(pd.m_subPopSize);
				pd1.m_subPopNames.swap(pd.m_subPopNames);
				pd1.m_genotype.swap(pd.m_genotype);
				pd1.m_info.swap(pd.m_info);
				pd1.m_inds.swap(pd.m_inds);
				std::swap(pd1.m_indOrdered, pd.m_indOrdered);
			}
		}
	}
	for (size_t genIdx = gens.size(); genIdx <= m_ancestralPops.size(); ++genIdx)
		m_ancestralPops.pop_back();
	m_curAncestralGen = 0;
}


void Population::useAncestralGen(ssize_t idx)
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), RuntimeError, "Can not switch ancestral generation with an activated virtual subpopulation");

	if (idx == m_curAncestralGen)
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
	DBG_ASSERT(static_cast<size_t>(idx) <= m_ancestralPops.size(),
		ValueError, "Ancestry generation " + toStr(idx) + " does not exist.");

	// now idx should be at least 1
	m_curAncestralGen = static_cast<int>(idx);
	// swap  1 ==> 0, 2 ==> 1

	popData & pd = m_ancestralPops[m_curAncestralGen - 1];
	pd.swap(*this);
	m_popSize = m_inds.size();
	setSubPopStru(m_subPopSize, m_subPopNames);
}


//template<class Archive>
void Population::save(boost::archive::text_oarchive & ar, const unsigned int /* version */) const
{
	// deep adjustment: everyone in order
	const_cast<Population *>(this)->syncIndPointers();

	ar & ModuleMaxAllele;

	DBG_DO(DBG_POPULATION, cerr << "Handling geno structure" << endl);
	// GenoStructure genoStru = this->genoStru();
	ar & genoStru();

	ar & m_subPopSize;
	ar & m_subPopNames;
	DBG_DO(DBG_POPULATION, cerr << "Handling genotype" << endl);
#ifdef BINARYALLELE
	size_t size = m_genotype.size();
	ar & size;
	ConstGenoIterator ptr = m_genotype.begin();
	WORDTYPE data = 0;
	for (size_t i = 0; i < size; ++i) {
		data |= (*ptr++) << (i % 32);
		// end of block of end of data
		if (i % 32 == 31 || i == size - 1) {
			// on 64 systems, the upper 32bit of the data might not be 0
			// so we need to clear the upper 32bit so that the genotype saved
			// on a 64bit system can be loaded on 32bit systems.
			data &= 0xFFFFFFFF;
			ar & data;
			data = 0;
		}
	}
#else
	ar & m_genotype;
#endif
	DBG_DO(DBG_POPULATION, cerr << "Handling information" << endl);
	ar & m_info;
	DBG_DO(DBG_POPULATION, cerr << "Handling Individuals" << endl);
	ar & m_inds;
	DBG_DO(DBG_POPULATION, cerr << "Handling ancestral populations" << endl);
	ar & m_ancestralGens;
	size_t sz = m_ancestralPops.size();
	ar & sz;
	for (size_t i = 0; i < m_ancestralPops.size(); ++i) {
		const_cast<Population *>(this)->useAncestralGen(i + 1);
		// need to make sure ancestral pop also in order
		const_cast<Population *>(this)->syncIndPointers();
		ar & m_subPopSize;
		ar & m_subPopNames;
#ifdef BINARYALLELE
		size_t size = m_genotype.size();
		ar & size;
		ptr = m_genotype.begin();
		WORDTYPE data = 0;
		for (size_t i = 0; i < size; ++i) {
			data |= (*ptr++) << (i % 32);
			// end of block of end of data
			if (i % 32 == 31 || i == size - 1) {
				data &= 0xFFFFFFFF;
				ar & data;
				data = 0;
			}
		}
#else
		ar & m_genotype;
#endif
		ar & m_info;
		ar & m_inds;
	}
	const_cast<Population *>(this)->useAncestralGen(0);

	// save shared variables as string.
	// note that many format are not supported.
	DBG_DO(DBG_POPULATION, cerr << "Handling shared variables" << endl);
	string vars = varsAsString();
	ar & vars;
}


//template<class Archive>
void Population::load(boost::archive::text_iarchive & ar, const unsigned int /* version */)
{
	size_t ma;
	ar & ma;

	DBG_WARNIF(ma > ModuleMaxAllele, "Warning: The population is saved in library with more allele states. \n"
		                             "Unless all alleles are less than " + toStr(ModuleMaxAllele) +
		", you should use the modules used to save this file. (c.f. simuOpt.setOptions()\n");

	GenoStructure stru;
	DBG_DO(DBG_POPULATION, cerr << "Handling geno structure" << endl);
	ar & stru;
	ar & m_subPopSize;
	ar & m_subPopNames;
	DBG_DO(DBG_POPULATION, cerr << "Handling genotype" << endl);

#ifdef BINARYALLELE
	// binary from binary
	if (ma == 1) {
		size_t size;
		ar & size;
		m_genotype.resize(size);
		GenoIterator ptr = m_genotype.begin();
		WORDTYPE data = 0;
		for (size_t i = 0; i < size; ++i) {
			if (i % 32 == 0)
				ar & data;
			*ptr++ = (data & (1UL << (i % 32))) != 0;
		}
	}
	// binary from others (long types)
	else {
		DBG_DO(DBG_POPULATION, cerr << "Load bin from long. " << endl);
		vectoru tmpgeno;
		ar & tmpgeno;
		m_genotype.resize(tmpgeno.size());
		for (size_t i = 0; i < tmpgeno.size(); ++i)
			m_genotype[i] = ToAllele(tmpgeno[i]);
	}
#else
	// long from binary
	if (ma == 1) {
		// for version 2 and higher, archive in 32bit blocks.
		size_t size;
		ar & size;
		m_genotype.resize(size);
		GenoIterator ptr = m_genotype.begin();
		WORDTYPE data = 0;
		for (size_t i = 0; i < size; ++i) {
			if (i % 32 == 0)
				ar & data;
			*ptr++ = (data & (1UL << (i % 32))) != 0;
		}
	}                                                                                   // if ma == 1
	else {                                                                              // for non-binary types, ...
		DBG_DO(DBG_POPULATION, cerr << "Load long from long. " << endl);
		// long from long
		ar & m_genotype;
	}
#endif

	DBG_DO(DBG_POPULATION, cerr << "Handling info" << endl);
	ar & m_info;

	DBG_DO(DBG_POPULATION, cerr << "Handling Individuals" << endl);
	ar & m_inds;

	// set genostructure, check duplication
	// we can not use setGenoStruIdx since stru may be new.
	this->setGenoStructure(stru);

	m_popSize = accumulate(m_subPopSize.begin(), m_subPopSize.end(), size_t(0));

	DBG_FAILIF(m_info.size() != m_popSize * infoSize(), ValueError, "Wgong size of info vector");

	if (m_popSize != m_inds.size()) {
		throw ValueError("Number of individuals does not match population size.\n"
			             "Please use the same (binary, short or long) module to save and load files.");
	}

	DBG_DO(DBG_POPULATION, cerr << "Reconstruct individual genotype" << endl);
	m_subPopIndex.resize(m_subPopSize.size() + 1);
	size_t i = 1;
	for (m_subPopIndex[0] = 0; i <= m_subPopSize.size(); ++i)
		m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];

	// assign genotype location and set structure information for individuals
	GenoIterator ptr = m_genotype.begin();
	size_t step = genoSize();
	InfoIterator infoPtr = m_info.begin();
	size_t infoStep = infoSize();
	for (size_t i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		m_inds[i].setGenoStruIdx(genoStruIdx());
		m_inds[i].setGenoPtr(ptr);
		m_inds[i].setInfoPtr(infoPtr);
	}
	m_ancestralGens = 0;
	m_ancestralPops.clear();

	// ancestry populations
	DBG_DO(DBG_POPULATION, cerr << "Handling ancestral populations" << endl);
	ar & m_ancestralGens;
	size_t na;
	ar & na;
	for (size_t ap = 0; ap < na; ++ap) {
		popData pd;
		ar & pd.m_subPopSize;
		ar & pd.m_subPopNames;
#ifdef BINARYALLELE
		// binary from binary
		if (ma == 1) {
			DBG_DO(DBG_POPULATION, cerr << "Load bin from bin. " << endl);
			size_t size;
			ar & size;
			pd.m_genotype.resize(size);
			ptr = pd.m_genotype.begin();
			WORDTYPE data = 0;
			for (size_t i = 0; i < size; ++i) {
				if (i % 32 == 0)
					ar & data;
				*ptr++ = (data & (1UL << i % 32)) != 0;
			}
		} else {
			DBG_DO(DBG_POPULATION, cerr << "Load bin from long. " << endl);
			// binary from long types
			vector<unsigned char> tmpgeno;
			ar & tmpgeno;
			pd.m_genotype.resize(tmpgeno.size());
			for (size_t i = 0; i < tmpgeno.size(); ++i)
				pd.m_genotype[i] = ToAllele(tmpgeno[i]);
		}
#else
		if (ma == 1) {
			// long type from binary
			size_t size;
			ar & size;
			pd.m_genotype.resize(size);
			ptr = pd.m_genotype.begin();
			WORDTYPE data = 0;
			for (size_t i = 0; i < size; ++i) {
				if (i % 32 == 0)
					ar & data;
				*ptr++ = (data & (1UL << i % 32)) != 0;
			}
		} else {
			DBG_DO(DBG_POPULATION, cerr << "Load long from long. " << endl);
			// long type from long type.
			ar & pd.m_genotype;
		}
#endif
		ar & pd.m_info;
		ar & pd.m_inds;
		// set pointer after copy this thing again (push_back)
		m_ancestralPops.push_back(pd);
		// now set pointers
		popData & p = m_ancestralPops.back();
		// set pointers
		vector<Individual> & inds = p.m_inds;
		size_t ps = inds.size();
		ptr = p.m_genotype.begin();
		infoPtr = p.m_info.begin();

		for (size_t i = 0; i < ps; ++i, ptr += step, infoPtr += infoStep) {
			inds[i].setGenoPtr(ptr);
			inds[i].setInfoPtr(infoPtr);
			// set new genoStructure
			inds[i].setGenoStruIdx(genoStruIdx());
		}
	}

	// load vars from string
	DBG_DO(DBG_POPULATION, cerr << "Handling shared variables" << endl);
	string vars;
	ar & vars;
	varsFromString(vars);

	setIndOrdered(true);
}


void Population::save(const string & filename) const
{
	boost::iostreams::filtering_ostream ofs;

	ofs.push(boost::iostreams::gzip_compressor());
	ofs.push(boost::iostreams::file_sink(filename, std::ios::binary));

	if (!ofs)
		throw ValueError("Can not open file " + filename);

	boost::archive::text_oarchive oa(ofs);
	oa << *this;
}


void Population::load(const string & filename)
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
		throw ValueError("Failed to load Population " + filename + ".\n");
	}
}


PyObject * Population::vars(vspID vsp)
{
	if (!vsp.valid()) {
		Py_INCREF(m_vars.dict());
		return m_vars.dict();
	}

	DBG_ASSERT(static_cast<size_t>(vsp.subPop()) < numSubPop(),
		IndexError, "Subpop index out of range of 0 ~ " + toStr(numSubPop() - 1));

	if (!m_vars.hasVar("subPop"))
		throw ValueError("Population local namespace does not have key 'subPop'. "
			             "You may forgot to call the Stat operator, or use the 'vars' parameter "
			             "to generate subpopulation-specific statistics.");

	PyObject * spObj = m_vars.getVar("subPop");
	// vsp? A tube with (sp, vsp)
	PyObject * key = NULL;
	if (vsp.isVirtual())
		key = Py_BuildValue("(ii)", vsp.subPop(), vsp.virtualSubPop());
	else
		key = PyInt_FromLong(static_cast<int>(vsp.subPop()));

	spObj = PyDict_GetItem(spObj, key);

	if (spObj == NULL)
		throw ValueError("Statistics for specified (virtual) subpopulation does not exist.");

	Py_INCREF(spObj);
	return spObj;
}


// The same as vars(), but without increasing
// reference count.
PyObject * Population::dict(vspID vsp)
{
	if (!vsp.valid())
		return m_vars.dict();

	DBG_ASSERT(static_cast<size_t>(vsp.subPop()) < numSubPop(),
		IndexError, "Subpop index out of range of 0 ~ " + toStr(numSubPop() - 1));

	if (!m_vars.hasVar("subPop"))
		throw ValueError("Population local namespace does not have key 'subPop'. "
			             "You may forgot to call the Stat operator, or use the 'vars' parameter "
			             "to generate subpopulation-specific statistics.");

	PyObject * spObj = m_vars.getVar("subPop");
	// vsp? A tube with (sp, vsp)
	PyObject * key = NULL;
	if (vsp.isVirtual())
		key = Py_BuildValue("(ii)", vsp.subPop(), vsp.virtualSubPop());
	else
		key = PyInt_FromLong(static_cast<int>(vsp.subPop()));

	spObj = PyDict_GetItem(spObj, key);

	if (spObj == NULL)
		throw ValueError("Statistics for specified (virtual) subpopulation does not exist.");

	return spObj;
}


void Population::syncIndPointers(bool infoOnly) const
{
	if (indOrdered())
		return;

	if (infoOnly) {
		DBG_DO(DBG_POPULATION, cerr << "Adjust info position " << endl);
		size_t is = infoSize();
		if (is == 0) {
			setIndOrdered(true);
			return;
		}
		vectorf tmpInfo(m_popSize * is);
		vectorf::iterator infoPtr = tmpInfo.begin();

		IndIterator ind = const_cast<Population *>(this)->indIterator();
		for (; ind.valid(); ++ind) {
			copy(ind->infoBegin(), ind->infoEnd(), infoPtr);
			ind->setInfoPtr(infoPtr);
			infoPtr += is;
		}
		const_cast<Population *>(this)->m_info.swap(tmpInfo);
	} else {
		DBG_DO(DBG_POPULATION, cerr << "Adjust geno and info position " << endl);

		size_t sz = genoSize();
		size_t is = infoSize();
		vectora tmpGenotype(m_popSize * genoSize());
		vectorf tmpInfo(m_popSize * infoSize());
		vectora::iterator it = tmpGenotype.begin();
		vectorf::iterator infoPtr = tmpInfo.begin();

		IndIterator ind = const_cast<Population *>(this)->indIterator();
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
		const_cast<Population *>(this)->m_genotype.swap(tmpGenotype);
		const_cast<Population *>(this)->m_info.swap(tmpInfo);
	}
	setIndOrdered(true);
}


Population & loadPopulation(const string & file)
{
	Population * p = new Population();

	p->load(file);
	return *p;
}


}


