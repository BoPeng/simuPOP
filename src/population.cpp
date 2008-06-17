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

#include "population.h"
#include "virtualSubPop.h"

// for file compression
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

namespace io = boost::iostreams;

namespace simuPOP {

population::population(const vectorlu & size,
                       float ploidy,
                       const vectoru & loci,
                       bool sexChrom,
                       const vectorf & lociPos,
                       int ancestralDepth,
                       const vectorstr & chromNames,
                       const vectorstr & alleleNames,
                       const vectorstr & lociNames,
                       UINT maxAllele,
                       const vectorstr & infoFields)
	:
	GenoStruTrait(),
	m_popSize(0),
	m_numSubPop(size.size()),
	m_subPopSize(size),
	m_subPopIndex(size.size() + 1),
	m_vspSplitter(NULL),
	m_genotype(0),                                                                          // resize later
	m_info(0),
	m_inds(0),                                                                              // default constructor will be called.
	m_ancestralDepth(ancestralDepth),
	m_vars(NULL, true),                                                                     // invalid shared variables initially
	m_ancestralPops(0),                                                                     // no history first
	m_rep(-1),
	m_grp(-1),
	m_gen(0),
	m_curAncestralGen(0),
	m_indOrdered(true),
	m_selectionFlags()
{
	DBG_FAILIF(maxAllele > ModuleMaxAllele, ValueError,
		"maxAllele is bigger than maximum allowed allele state of this library (" + toStr(ModuleMaxAllele) +
		")\nPlease use simuOpt.setOptions(alleleType='long') to use the long allele version of simuPOP.");

	DBG_FAILIF(maxAllele == 0, ValueError,
		"maxAllele should be at least 1 (0,1 two states). ");

	DBG_DO(DBG_POPULATION, cout << "Constructor of population is called\n");

	DBG_FAILIF(m_subPopSize.size() > MaxSubPopID, ValueError,
		"Number of subpopulations exceed maximum allowed subpopulation numbers");

	// get a GenoStructure with parameters. GenoStructure may be shared by some populations
	// a whole set of functions ploidy() etc in GenoStruTriat can be used after this step.
	DBG_FAILIF(static_cast<UINT>(ploidy) * 1.0 != ploidy && fcmp_ne(ploidy, Haplodiploid),
		ValueError, "Only integer ploidy number or Haplodiploid can be specified");

	setGenoStructure(fcmp_eq(ploidy, Haplodiploid) ? 2 : static_cast<UINT>(ploidy),
		loci, sexChrom, fcmp_eq(ploidy, Haplodiploid), lociPos, chromNames, alleleNames,
		lociNames, maxAllele, infoFields);

	DBG_DO(DBG_DEVEL, cout << "individual size is " << sizeof(individual) << '+'
		                   << sizeof(Allele) << '*' << genoSize() << endl
		                   << ", infoPtr: " << sizeof(double *)
		                   << ", GenoPtr: " << sizeof(Allele *) << ", Flag: " << sizeof(unsigned char)
		                   << ", plus genoStru"
		                   << "\ngenoSize " << genoSize()
		                   << endl);

	// m_popSize will be defined in fitSubPopStru
	if (m_subPopSize.empty())
		m_subPopSize.resize(1, 0);
	fitSubPopStru(m_subPopSize);
	// set local variable
	setRep(-1);
	setGrp(-1);
}


/// destroy a population
population::~population()
{
	if (m_vspSplitter)
		delete m_vspSplitter;

	DBG_DO(DBG_POPULATION,
		cout << "Destructor of population is called" << endl);
}


population::population(const population & rhs) :
	GenoStruTrait(rhs),
	m_popSize(rhs.m_popSize),
	m_numSubPop(rhs.m_numSubPop),
	m_subPopSize(rhs.m_subPopSize),
	m_subPopIndex(rhs.m_subPopIndex),
	m_vspSplitter(NULL),
	m_genotype(0),
	m_info(0),
	m_inds(0),
	m_ancestralDepth(rhs.m_ancestralDepth),
	m_vars(rhs.m_vars),                                                                     // variables will be copied
	m_rep(-1),                                                                              // rep is set to -1 for new pop (until simulator really set them
	m_grp(-1),
	m_gen(0),
	m_curAncestralGen(rhs.m_curAncestralGen),
	m_indOrdered(true),
	m_selectionFlags()
{
	DBG_DO(DBG_POPULATION,
		cout << "Copy constructor of population is called" << endl);

	try {
		m_inds.resize(rhs.m_popSize);
		m_genotype.resize(m_popSize * genoSize());
		// have 0 length for mpi/non-head node
		m_info.resize(rhs.m_popSize * infoSize());
	} catch (...) {
		throw OutOfMemory("Memory allocation fail");
	}

	// individuals will always have the correct genostructure
	// by using their copied pointer
	// population, however, need to set this pointer correctly
	//
	setGenoStruIdx(rhs.genoStruIdx());

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
		cout << "Unable to copy ancestral populations. "
		     << "The popolation size may be too big." << endl
		     << "The population will still be usable but without any ancestral population stored." << endl;
		m_ancestralDepth = 0;
		m_ancestralPops.clear();
	}

	// copy virtual subpop splitters
	setVirtualSplitter(rhs.virtualSplitter());

	// set local variable
	setRep(-1);
	setGrp(-1);
}


void population::popData::swap(population & pop)
{
	pop.m_subPopSize.swap(m_subPopSize);
	pop.m_genotype.swap(m_genotype);
	pop.m_info.swap(m_info);
	pop.m_inds.swap(m_inds);
	std::swap(pop.m_indOrdered, m_indOrdered);
}


///
population * population::clone(int keepAncestralPops) const
{
	population * p = new population(*this);
	int oldDepth = m_ancestralDepth;

	if (keepAncestralPops >= 0)
		// try to remove excessive ancestra generations.
		p->setAncestralDepth(keepAncestralPops);
	p->setAncestralDepth(oldDepth);
	return p;
}


ULONG population::virtualSubPopSize(SubPopID subPop, SubPopID virtualSubPop) const
{
	CHECKRANGESUBPOP(subPop);
	// if there is no virtual subpopulation,
	if (hasActivatedVirtualSubPop() || virtualSubPop != InvalidSubPopID)
		return m_vspSplitter->size(*this, subPop, virtualSubPop);
	return subPopSize(subPop);
}


string population::virtualSubPopName(SubPopID subPop, SubPopID virtualSubPop) const
{
	CHECKRANGESUBPOP(subPop);
	DBG_ASSERT(virtualSubPop != InvalidSubPopID, ValueError,
		"Subpopulation id is not virtual");
	DBG_ASSERT(hasVirtualSubPop(), ValueError,
		"There is no virtual subpopulation in subpop " + toStr(subPop));
	return m_vspSplitter->name(virtualSubPop);
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


void population::activateVirtualSubPop(SubPopID subPop, SubPopID virtualSubPop,
                                       vspSplitter::activateType type)
{
	CHECKRANGESUBPOP(subPop);
	DBG_ASSERT(virtualSubPop != InvalidSubPopID, ValueError, "Given virtual subpopulation ID is wrong");
	DBG_ASSERT(hasVirtualSubPop(), ValueError,
		"Subpopulation " + toStr(subPop) + " has no virtual subpopulations");
	m_vspSplitter->activate(*this, subPop, virtualSubPop, type);
	DBG_ASSERT(type !=  vspSplitter::Visible || 
		m_vspSplitter->activatedSubPop() == subPop, SystemError,
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
		DBG_DO(DBG_POPULATION, cout << "Genotype structures are different" << endl);
		return 1;
	}

	if (popSize() != rhs.popSize() ) {
		DBG_DO(DBG_POPULATION, cout << "Population sizes are different" << endl);
		return 1;
	}

	for (ULONG i = 0, iEnd = popSize(); i < iEnd; ++i)
		if (m_inds[i] != rhs.m_inds[i]) {
			DBG_DO(DBG_POPULATION, cout << "Individuals are different" << endl);
			return 1;
		}

	return 0;
}

individual & population::ancestor(ULONG ind, UINT gen)
{
	DBG_FAILIF(gen > m_ancestralPops.size(), IndexError,
		"Ancestray generation " + toStr(gen) + " does not exist");
	if (gen == m_curAncestralGen)
		return this->ind(ind);
	UINT idx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
	DBG_FAILIF(ind > m_ancestralPops[idx].m_inds.size(),
		IndexError, "Individual index out of range");
	return m_ancestralPops[idx].m_inds[ind];
}


const individual & population::ancestor(ULONG ind, UINT gen) const
{
	DBG_FAILIF(gen > m_ancestralPops.size(), IndexError,
		"Ancestray generation " + toStr(gen) + " does not exist");
	if (gen == m_curAncestralGen)
		return this->ind(ind);
	UINT idx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
	DBG_FAILIF(ind > m_ancestralPops[idx].m_inds.size(),
		IndexError, "Individual index out of range");
	return m_ancestralPops[idx].m_inds[ind];
}


individual & population::ancestor(ULONG ind, UINT subPop, UINT gen)
{
	DBG_FAILIF(gen > m_ancestralPops.size(), IndexError,
		"Ancestray generation " + toStr(gen) + " does not exist");
	if (gen == m_curAncestralGen)
		return this->ind(ind, subPop);
	UINT idx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
	DBG_FAILIF(subPop > m_ancestralPops[idx].m_subPopSize.size(),
		IndexError, "subpopulation index out of range");
	DBG_FAILIF(ind > m_ancestralPops[idx].m_subPopSize[subPop],
		IndexError, "Individual index out of range");
	ULONG shift = 0;
	if (subPop > 0) {
		for (size_t i = 0; i < subPop; ++i)
			shift += m_ancestralPops[idx].m_subPopSize[i];
	}
	return m_ancestralPops[idx].m_inds[shift + ind];
}


const individual & population::ancestor(ULONG ind, UINT subPop, UINT gen) const
{
	DBG_FAILIF(gen > m_ancestralPops.size(), IndexError,
		"Ancestray generation " + toStr(gen) + " does not exist");
	if (gen == m_curAncestralGen)
		return this->ind(ind, subPop);
	UINT idx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
	DBG_FAILIF(subPop > m_ancestralPops[idx].m_subPopSize.size(),
		IndexError, "subpopulation index out of range");
	DBG_FAILIF(ind > m_ancestralPops[idx].m_subPopSize[subPop],
		IndexError, "Individual index out of range");
	ULONG shift = 0;
	if (subPop > 0) {
		for (size_t i = 0; i < subPop; ++i)
			shift += m_ancestralPops[idx].m_subPopSize[i];
	}
	return m_ancestralPops[idx].m_inds[shift + ind];
}


PyObject * population::arrGenotype(bool order)
{
	if (order)
		sortIndividuals();
	// directly expose values. Do not copy data over.
	return Allele_Vec_As_NumArray(m_genotype.begin(), m_genotype.end());
}


void population::setIndSubPopID(const vectori & id, bool ancestralPops)
{
	UINT oldGen = ancestralGen();
	size_t sz = id.size();
	for (UINT anc = 0; anc <= ancestralDepth(); ++anc) {
		if (!ancestralPops && anc != oldGen)
			continue;
		useAncestralPop(anc);
		for (ULONG it = 0; it < m_popSize; ++it)
			ind(it).setSubPopID(static_cast<SubPopID>(id[it % sz]));
	}
	useAncestralPop(oldGen);
}


void population::setIndSubPopIDWithID(bool ancestralPops)
{
	UINT oldGen = ancestralGen();
	for (UINT anc = 0; anc <= ancestralDepth(); ++anc) {
		if (!ancestralPops && anc != oldGen)
			continue;
		useAncestralPop(anc);
		for (UINT i = 0, iEnd = numSubPop(); i < iEnd;  ++i)
			for (IndIterator it = indBegin(i); it.valid();  ++it)
				it->setSubPopID(i);
	}
	useAncestralPop(oldGen);
}


// get the whole genotype.
// individuals will be in order before exposing
// their genotypes.
//
// if order: keep order
// otherwise: respect subpop structure
PyObject * population::arrGenotype(UINT subPop, bool order)
{
	CHECKRANGESUBPOP(subPop);
	sortIndividuals();
	return Allele_Vec_As_NumArray(genoBegin(subPop, order), genoEnd(subPop, order));
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
		for (ConstIndIterator it = indBegin(); it.valid(); ++it) {
			DBG_ASSERT(it->genoPtr() >= gb && it->genoPtr() < ge, SystemError,
				msg + "Wrong genotype pointer");
		}
	}
	if (infoSize() > 0) {
		for (ConstIndIterator it = indBegin(); it.valid(); ++it) {
			DBG_ASSERT(it->infoPtr() >= ib && it->infoPtr() < ie, SystemError,
				msg + "Wrong information field pointer. (number of information fields: "
				+ toStr(infoSize()) + ")");
		}
	}
#endif
}


void population::fitSubPopStru(const vectorlu & newSubPopSizes)
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
			throw OutOfMemory("Memory allocation fail. (popSize=" + toStr(m_popSize) + ")");
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

	setSubPopStru(newSubPopSizes);
}


void population::setSubPopStru(const vectorlu & newSubPopSizes)
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	DBG_ASSERT(accumulate(newSubPopSizes.begin(), newSubPopSizes.end(), 0UL) == m_popSize, ValueError,
		"Overall population size should not be changed in setSubPopStru.");

	if (newSubPopSizes.empty()) {
		m_numSubPop = 1;
		m_subPopSize = vectorlu(1, 0);
	} else {
		m_numSubPop = newSubPopSizes.size();
		m_subPopSize = newSubPopSizes;
	}
	m_subPopIndex.resize(m_numSubPop + 1);

	// build subPop index
	UINT i = 1;
	for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
		m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];
}


void population::setSubPopByIndID(vectori id)
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	if (!id.empty()) {
		DBG_ASSERT(id.size() == m_popSize, ValueError,
			"Info should have the same length as pop size");
		for (ULONG it = 0; it < m_popSize; ++it)
			ind(it).setSubPopID(id[it]);
	}

	DBG_DO(DBG_POPULATION, cout << "Sorting individuals." << endl);
	// sort individuals first
	std::sort(indBegin(), indEnd());
	setIndOrdered(false);

	// sort individuals first
	// remove individuals with negative index.
	if (indBegin()->subPopID() < 0) {
		// popsize etc will be changed.
		ULONG newPopSize = m_popSize;
		IndIterator it = indBegin();
		for (; it.valid();  ++it) {
			if (it->subPopID() < 0)
				newPopSize-- ;
			else
				break;
		}
		// 'it' now point to the one with positive subPopID()

		DBG_DO(DBG_POPULATION, cout << "New pop size" << newPopSize << endl);

		// allocate new genotype and inds
		vectora newGenotype(genoSize() * newPopSize);
		vectorinfo newInfo(newPopSize * infoSize());
		vector<individual> newInds(newPopSize);

		DBG_ASSERT(indEnd() == it + newPopSize, SystemError,
			"Pointer misplaced. ");

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
		m_numSubPop = 1;
		m_subPopSize.resize(1, 0);
		m_subPopIndex.resize(2);
	} else {
		// reset indexes etc.
		m_numSubPop = static_cast<UINT>(m_inds.back().subPopID()) + 1;
		m_subPopSize.resize(m_numSubPop);
		m_subPopIndex.resize(m_numSubPop + 1);

		// check subpop size
		fill(m_subPopSize.begin(), m_subPopSize.end(), 0);
		for (IndIterator it = indBegin(); it.valid();  ++it)
			m_subPopSize[ static_cast<UINT>(it->subPopID()) ]++;
	}
	// rebuild index
	size_t i = 1;
	for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
		m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];
}


void population::splitSubPop(UINT which, vectorlu sizes, vectoru subPopID)
{
	DBG_ASSERT(accumulate(sizes.begin(), sizes.end(), 0UL) == subPopSize(which),
		ValueError,
		"Sum of subpopulation sizes does not equal to the size of subpopulation to be splitted.");

	DBG_FAILIF(!subPopID.empty() && subPopID.size() != sizes.size(), ValueError,
		"If subPopID is given, it should have the same length as subPOP");

	if (sizes.size() == 1)
		return;

	// set initial info
	setIndSubPopIDWithID();

	UINT spID;
	if (subPopID.empty())  // starting sp number
		spID = which;
	else {
		spID = subPopID[0];
		DBG_WARNING(spID != which && spID < numSubPop(),
			"new subpop ID is already used. You are effectively merging two subpopulations")
	}
	ULONG sz = 0;                                                                     // idx within subpop
	size_t newSPIdx = 0;
	for (IndIterator ind = indBegin(which); ind.valid(); ++ind) {
		if (sz == sizes[newSPIdx]) {
			sz = 0;
			newSPIdx++;
			if (subPopID.empty())
				spID = numSubPop() + newSPIdx - 1;
			else {
				DBG_WARNING(subPopID[newSPIdx] != which && subPopID[newSPIdx] < numSubPop(),
					"new subpop ID is already used. You are effectively merging two subpopulations")
				spID = subPopID[newSPIdx];
			}
		}
		ind->setSubPopID(spID);
		sz++;
	}
	setSubPopByIndID();
}


void population::splitSubPopByProportion(UINT which, vectorf proportions, vectoru subPopID)
{
	DBG_ASSERT(fcmp_eq(accumulate(proportions.begin(), proportions.end(), 0.), 1.), ValueError,
		"Proportions do not add up to one.");

	if (proportions.size() == 1)
		return;

	ULONG spSize = subPopSize(which);
	vectorlu subPop(proportions.size());
	for (size_t i = 0; i < proportions.size() - 1; ++i)
		subPop[i] = static_cast<ULONG>(floor(spSize * proportions[i]));
	// to avoid round off problem, calculate the last subpopulation
	subPop[ subPop.size() - 1] = spSize - accumulate(subPop.begin(), subPop.end() - 1, 0L);
	splitSubPop(which, subPop, subPopID);
}


void population::removeEmptySubPops()
{
	// if remove empty subpops
	UINT newSPNum = m_numSubPop;
	vectorlu newSPSize;

	for (size_t sp = 0; sp < m_numSubPop; ++sp) {
		if (m_subPopSize[sp] == 0)
			newSPNum--;
		else
			newSPSize.push_back(m_subPopSize[sp]);
	}
	m_numSubPop = newSPNum;
	m_subPopSize.swap(newSPSize);
	m_subPopIndex.resize(m_numSubPop + 1);
	// rebuild index
	size_t i = 1;
	for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
		m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];
}


void population::removeSubPops(const vectoru & subPops, bool shiftSubPopID, bool removeEmptySubPops)
{
#ifndef OPTIMIZED
	// check if subPops are valid
	for (vectoru::const_iterator sp = subPops.begin(); sp < subPops.end(); ++sp) {
		DBG_WARNING(*sp >= m_numSubPop, "Subpopulation " + toStr(*sp) + " does not exist.");
	}
#endif
	setIndSubPopIDWithID();

	int shift = 0;
	for (size_t sp = 0; sp < m_numSubPop; ++sp) {
		if (find(subPops.begin(), subPops.end(), sp) != subPops.end()) {
			shift++;
			RawIndIterator ind = rawIndBegin(sp);
			RawIndIterator ind_end = rawIndEnd(sp);
			for (; ind != ind_end; ++ind)
				ind->setSubPopID(-1); // remove
		} else {
			// other subpop shift left
			if (shiftSubPopID) {
				RawIndIterator ind = rawIndBegin(sp);
				RawIndIterator ind_end = rawIndEnd(sp);
				for (; ind != ind_end; ++ind)
					ind->setSubPopID(sp - shift); // shift left
			}
		}
	}

	UINT pendingEmptySubPops = 0;
	for (UINT i = m_numSubPop - 1; i >= 0 && (subPopSize(i) == 0
	                                          || find(subPops.begin(), subPops.end(), i) != subPops.end()); --i, ++pendingEmptySubPops) ;
	setSubPopByIndID();
	// what to do with pending empty subpops?
	if (pendingEmptySubPops != 0 && !removeEmptySubPops) {
		vectorlu spSizes = subPopSizes();
		for (UINT i = 0; i < pendingEmptySubPops; ++i)
			spSizes.push_back(0);
		setSubPopStru(spSizes);
	}
	if (removeEmptySubPops)
		this->removeEmptySubPops();
}


void population::removeIndividuals(const vectoru & inds, int subPop, bool removeEmptySubPops)
{
	setIndSubPopIDWithID();
	if (subPop == -1) {
		for (size_t i = 0; i < inds.size(); ++i)
			ind(inds[i]).setSubPopID(-1); // remove
	} else {
		for (size_t i = 0; i < inds.size(); ++i)
			// remove
			ind(inds[i], subPop).setSubPopID(-1);
	}

	int oldNumSP = numSubPop();
	setSubPopByIndID();
	int pendingEmptySubPops = oldNumSP - numSubPop();
	// what to do with pending empty subpops?
	if (pendingEmptySubPops != 0 && !removeEmptySubPops) {
		vectorlu spSizes = subPopSizes();
		for (int i = 0; i < pendingEmptySubPops; ++i)
			spSizes.push_back(0);
		setSubPopStru(spSizes);
	}
	if (removeEmptySubPops)
		this->removeEmptySubPops();
}


void population::mergeSubPops(vectoru subPops, bool removeEmptySubPops)
{
	// set initial info
	setIndSubPopIDWithID();

	// merge all subpopulations
	if (subPops.empty()) {
		// [ popSize() ]
		vectorlu sz(1, popSize());
		setSubPopStru(sz);
		return;
	}

	UINT id = subPops[0];
	for (UINT sp = 0; sp < numSubPop(); ++sp) {
		if (find(subPops.begin(), subPops.end(), sp) != subPops.end())
			for (IndIterator ind = indBegin(sp); ind.valid(); ++ind)
				ind->setSubPopID(id);
	}
	int oldNumSP = numSubPop();
	setSubPopByIndID();
	int pendingEmptySubPops = oldNumSP - numSubPop();
	// what to do with pending empty subpops?
	if (pendingEmptySubPops != 0 && !removeEmptySubPops) {
		vectorlu spSizes = subPopSizes();
		for (int i = 0; i < pendingEmptySubPops; ++i)
			spSizes.push_back(0);
		setSubPopStru(spSizes);
	}

	if (removeEmptySubPops)
		this->removeEmptySubPops();
}


void population::mergePopulationPerGen(const population & pop, const vectorlu & newSubPopSizes)
{
	DBG_FAILIF(genoStruIdx() != pop.genoStruIdx(), ValueError,
		"Merged population should have the same genotype structure");
	// calculate new population size
	vectorlu newSS;
	newSS.insert(newSS.end(), m_subPopSize.begin(), m_subPopSize.end());
	newSS.insert(newSS.end(), pop.m_subPopSize.begin(), pop.m_subPopSize.end());
	// new population size
	ULONG newPopSize = accumulate(newSS.begin(), newSS.end(), 0UL);
	DBG_FAILIF(!newSubPopSizes.empty() && accumulate(newSubPopSizes.begin(), newSubPopSizes.end(), 0UL) != newPopSize,
		ValueError, "newSubPopSizes should not change overall population size");

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
	for (ULONG i = 0; i < popSize(); ++i)
		newInds[i].copyFrom(ind(i));
	ULONG start = popSize();
	for (ULONG i = 0; i < pop.popSize(); ++i)
		newInds[i + start].copyFrom(pop.ind(i));
	// now, switch!
	m_genotype.swap(newGenotype);
	m_info.swap(newInfo);
	m_inds.swap(newInds);
	m_popSize = newPopSize;
	setIndOrdered(true);
	if (newSubPopSizes.empty())
		m_subPopSize = newSS;
	else
		m_subPopSize = newSubPopSizes;
	m_numSubPop = m_subPopSize.size();
	// rebuild index
	m_subPopIndex.resize(m_numSubPop + 1);
	size_t j = 1;
	for (m_subPopIndex[0] = 0; j <= m_numSubPop; ++j)
		m_subPopIndex[j] = m_subPopIndex[j - 1] + m_subPopSize[j - 1];
}


void population::mergePopulation(const population & pop, const vectorlu & newSubPopSizes,
                                 int keepAncestralPops)
{
	DBG_FAILIF(ancestralDepth() != pop.ancestralDepth(), ValueError,
		"Merged populations should have the same number of ancestral generations");
	UINT topGen;
	if (keepAncestralPops < 0 || static_cast<UINT>(keepAncestralPops) >= ancestralDepth())
		topGen = ancestralDepth();
	else
		topGen = keepAncestralPops;
	// go to the oldest generation
	useAncestralPop(topGen);
	const_cast<population &>(pop).useAncestralPop(topGen);
	mergePopulationPerGen(pop, newSubPopSizes);
	if (topGen > 0) {
		for (int depth = topGen - 1; depth >= 0; --depth) {
			useAncestralPop(depth);
			const_cast<population &>(pop).useAncestralPop(depth);
			mergePopulationPerGen(pop, newSubPopSizes);
		}
	}
	useAncestralPop(0);
	const_cast<population &>(pop).useAncestralPop(0);
}


void population::mergePopulationByLoci(const population & pop,
                                       const vectoru & newNumLoci, const vectorf & newLociPos, bool byChromosome)
{
	DBG_FAILIF(subPopSizes() != pop.subPopSizes(), ValueError,
		"Merged population should have the same number of individuals in each subpopulation");

	UINT numloci1 = totNumLoci();
	UINT numloci2 = pop.totNumLoci();

	DBG_FAILIF(!newNumLoci.empty() && accumulate(newNumLoci.begin(), newNumLoci.end(), 0U) != numloci1 + numloci2,
		ValueError, "Sum of newNumLoci should equal to " + toStr(numloci1 + numloci2));

	DBG_FAILIF(!newLociPos.empty() && newLociPos.size() != numloci1 + numloci2,
		ValueError, "newLociPos should have the length of combined total number of loci");

	// make copy of old genotype structure
	GenoStructure * gs1 = new GenoStructure(genoStru());
	GenoStructure * gs2 = new GenoStructure(pop.genoStru());

	// obtain new genotype structure and set it
	setGenoStructure(mergeGenoStru(pop.genoStruIdx(), byChromosome));
	// if byChromosome, we have to figure out how to copy loci
	vectoru source;
	if (byChromosome) {
		DBG_DO(DBG_POPULATION, cout << "New pos " << lociPos() << endl);

		for (size_t ch = 0; ch < numChrom(); ++ch) {
			size_t idx1 = 0;
			size_t idx2 = 0;
			for (size_t loc = 0; loc < numLoci(ch); ++loc) {
				double pos = locusPos(absLocusIndex(ch, loc));
				if (pos == gs1->locusPos(gs1->chromIndex(ch) + idx1)) {
					source.push_back(0);
					idx1++;
				} else if (pos == gs2->locusPos(gs2->chromIndex(ch) + idx2)) {
					source.push_back(1);
					idx2++;
				} else
					throw SystemError("Unable to determine position. Something wrong");
			}
		}
	}
	// remove copied genostructure
	delete gs1;
	delete gs2;

	DBG_FAILIF(ancestralDepth() != pop.ancestralDepth(), ValueError,
		"Merged populations should have the same number of ancestral generations");
	for (int depth = ancestralDepth(); depth >= 0; --depth) {
		useAncestralPop(depth);
		const_cast<population &>(pop).useAncestralPop(depth);
		//
		ULONG newPopGenoSize = genoSize() * m_popSize;
		vectora newGenotype(newPopGenoSize);

		if (byChromosome) {
			// merge chromosome by chromosome
			// copy data over
			GenoIterator ptr = newGenotype.begin();
			UINT pEnd = ploidy();
			for (ULONG i = 0; i < m_popSize; ++i) {
				// set new geno structure
				m_inds[i].setGenoStruIdx(genoStruIdx());
				GenoIterator ptr1 = m_inds[i].genoPtr();
				GenoIterator ptr2 = pop.m_inds[i].genoPtr();
				// new genotype
				m_inds[i].setGenoPtr(ptr);
				// copy each allele
				for (UINT p = 0; p < pEnd; ++p) {
					for (size_t j = 0; j < numloci1 + numloci2; ++j)
						if (source[j] == 0)
							*(ptr++) = *(ptr1++);
						else
							*(ptr++) = *(ptr2++);
				}
			}
		} else {
			// append pop2 chromosomes to the first one
			// copy data over
			GenoIterator ptr = newGenotype.begin();
			UINT pEnd = ploidy();
			for (ULONG i = 0; i < m_popSize; ++i) {
				// set new geno structure
				m_inds[i].setGenoStruIdx(genoStruIdx());
				GenoIterator ptr1 = m_inds[i].genoPtr();
				GenoIterator ptr2 = pop.m_inds[i].genoPtr();
				// new genotype
				m_inds[i].setGenoPtr(ptr);
				// copy each allele
				for (UINT p = 0; p < pEnd; ++p) {
					for (size_t j = 0; j < numloci1; ++j)
						*(ptr++) = *(ptr1++);
					for (size_t j = 0; j < numloci2; ++j)
						*(ptr++) = *(ptr2++);
				}
			}
		}
		m_genotype.swap(newGenotype);
	}

	// reset genoStructure again
	if (!newNumLoci.empty() || !newLociPos.empty())
		rearrangeLoci(newNumLoci, newLociPos);

	if (!indOrdered())
		// sort information only
		sortIndividuals(true);
}


void population::insertBeforeLocus(UINT idx, double pos, const string & name)
{
	vectoru v_idx(1, idx);
	vectorf v_pos(1, pos);

	if (name == string())
		return insertBeforeLoci(v_idx, v_pos, vectorstr());
	else
		return insertBeforeLoci(v_idx, v_pos, vectorstr(1, name));
}


void population::insertBeforeLoci(const vectoru & idx, const vectorf & pos, const vectorstr & names)
{
	// use loci to keep the position of old loci in the new structure
	vectoru loci(totNumLoci());

	// 0, 1, 2, 3, 4, 5, 6, ...
	for (size_t i = 0; i < totNumLoci(); ++i)
		loci[i] = i;
	// can insert multiple loci at the same location
	for (size_t i = 0; i < idx.size(); ++i) {
		CHECKRANGEABSLOCUS(idx[i]);
		// 0, 1, 2, 3, 4, 5, 6, ...
		// if idx[i] = 3 (before index 3)
		// 0, 1, 2, 4, 5, 6, 7, ...
		for (size_t j = idx[i]; j < totNumLoci(); ++j)
			loci[j]++;
	}

	// obtain new genotype structure and set it
	setGenoStructure(insertBeforeLociToGenoStru(idx, pos, names));

	for (int depth = ancestralDepth(); depth >= 0; --depth) {
		useAncestralPop(depth);
		//
		ULONG newPopGenoSize = genoSize() * m_popSize;
		vectora newGenotype(newPopGenoSize);

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
				for (vectoru::iterator loc = loci.begin();
				     loc != loci.end(); ++loc) {
					newPtr[*loc] = *(oldPtr++);
				}
				newPtr += totNumLoci();                   // next ploidy
			}
		}
		m_genotype.swap(newGenotype);
	}
	if (!indOrdered())
		sortIndividuals(true);
}


void population::insertAfterLocus(UINT idx, double pos, const string & name)
{
	vectoru v_idx(1, idx);
	vectorf v_pos(1, pos);

	if (name == string())
		return insertAfterLoci(v_idx, v_pos, vectorstr());
	else
		return insertAfterLoci(v_idx, v_pos, vectorstr(1, name));
}


void population::insertAfterLoci(const vectoru & idx, const vectorf & pos, const vectorstr & names)
{
	// use loci to keep the position of old loci in the new structure
	vectoru loci(totNumLoci());

	// 0, 1, 2, 3, 4, 5, 6, ...
	for (size_t i = 0; i < totNumLoci(); ++i)
		loci[i] = i;
	// can insert multiple loci at the same location
	for (size_t i = 0; i < idx.size(); ++i) {
		CHECKRANGEABSLOCUS(idx[i]);
		// 0, 1, 2, 3, 4, 5, 6, ...
		// if idx[i] = 3 (before index 3)
		// 0, 1, 2, 3, 5, 6, 7, ...
		if (idx[i] != totNumLoci() - 1)
			for (size_t j = idx[i] + 1; j < totNumLoci(); ++j)
				loci[j]++;
	}

	// obtain new genotype structure and set it
	setGenoStructure(insertAfterLociToGenoStru(idx, pos, names));

	for (int depth = ancestralDepth(); depth >= 0; --depth) {
		useAncestralPop(depth);
		//
		ULONG newPopGenoSize = genoSize() * m_popSize;
		vectora newGenotype(newPopGenoSize);

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
				for (vectoru::iterator loc = loci.begin();
				     loc != loci.end(); ++loc) {
					newPtr[*loc] = *(oldPtr++);
				}
				newPtr += totNumLoci();                   // next ploidy
			}
		}
		m_genotype.swap(newGenotype);
	}
	if (!indOrdered())
		sortIndividuals(true);
}


void population::resize(const vectorlu & newSubPopSizes, bool propagate)
{
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
			if ((j / spSize) > 0 && !propagate)
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
	for (m_subPopIndex[0] = 0; idx <= m_numSubPop; ++idx)
		m_subPopIndex[idx] = m_subPopIndex[idx - 1] + m_subPopSize[idx - 1];
}


void population::reorderSubPops(const vectoru & order, const vectoru & rank,
                                bool removeEmptySubPops)
{
	DBG_FAILIF(order.empty() && rank.empty(), ValueError,
		"Please specify one of order or rank.");

	DBG_FAILIF(!order.empty() && !rank.empty(), ValueError,
		"You can specify only one of order or rank.");

	if (removeEmptySubPops)
		this->removeEmptySubPops();

	if ( (!order.empty() && order.size() != m_numSubPop)
	    || (!rank.empty() && rank.size() != m_numSubPop))
		cout << "Warning: Given order or rank does not have the length of number of subpop." << endl;

	if (!order.empty()) {
		// alow order[i] > numSubPop(). In a special case, I have last empty subpop...
		for (size_t i = 0; i < order.size(); ++i) {
			if (order[i] >= numSubPop())
				continue;
			for (IndIterator ind = indBegin(order[i]); ind.valid(); ++ind)
				ind->setSubPopID(i);
		}
	} else {
		for (size_t i = 0; i < rank.size(); ++i) {
			if (i >= numSubPop())
				continue;
			for (IndIterator ind = indBegin(i); ind.valid(); ++ind)
				ind->setSubPopID(rank[i]);
		}
	}
	// reset ...
	setSubPopByIndID();
}


population & population::newPopByIndIDPerGen(const vectori & id, bool removeEmptySubPops)
{
	// determine the size of needed individuals
	vectorlu sz;

	if (!id.empty()) {
		DBG_ASSERT(id.size() == popSize(), ValueError, "Please assign id for each individual");
		for (ULONG i = 0; i != id.size(); ++i) {
			if (id[i] < 0)
				continue;
			if (static_cast<UINT>(id[i]) >= sz.size())
				sz.resize(id[i] + 1);
			sz[id[i]]++;
		}
	} else {
		for (UINT sp = 0; sp < numSubPop(); ++sp) {
			for (IndIterator it = indBegin(sp); it.valid(); ++it) {
				int indID = it->subPopID();
				if (indID < 0)
					continue;
				if (static_cast<UINT>(indID) >= sz.size())
					sz.resize(indID + 1);
				sz[indID]++;
			}
		}
	}
	DBG_DO(DBG_POPULATION, cout << "newPopByIndIDPerGen: New population size: " << sz << endl);

	// create a population with this size
	population * pop = new population(sz, ploidy(), numLoci(), sexChrom(), lociPos(), 0,
		chromNames(), alleleNames(), lociNames(), maxAllele(), infoFields());
	// copy individuals over
	IndIterator from = indBegin();
	vector<IndIterator> to;
	for (UINT sp = 0; sp < sz.size(); ++sp)
		to.push_back(pop->indBegin(sp));
	if (!id.empty()) {
		for (ULONG i = 0; i != id.size(); ++i) {
			if (id[i] >= 0) {
				to[id[i]]->copyFrom(ind(i));
				++to[id[i]];
			}
		}
	} else {
		for (; from.valid(); ++from) {
			int indID = from->subPopID();
			if (indID >= 0) {
				to[indID]->copyFrom(*from);
				++to[indID];
			}
		}
	}
	if (removeEmptySubPops)
		pop->removeEmptySubPops();
	return *pop;
}


/** form a new population according to info, info can be given directly */
population & population::newPopByIndID(int keepAncestralPops,
                                       const vectori & id, bool removeEmptySubPops)
{
	UINT topGen;

	if (keepAncestralPops < 0 || static_cast<UINT>(keepAncestralPops) >= ancestralDepth())
		topGen = ancestralDepth();
	else
		topGen = keepAncestralPops;
	// go to the oldest generation
	useAncestralPop(topGen);
	population & ret = newPopByIndIDPerGen(id, removeEmptySubPops);
	// prepare for push and discard
	ret.setAncestralDepth(topGen);
	if (topGen > 0) {
		for (int depth = topGen - 1; depth >= 0; --depth) {
			useAncestralPop(depth);
			ret.pushAndDiscard(newPopByIndIDPerGen(id, removeEmptySubPops));
		}
	}
	return ret;
}


void population::removeLoci(const vectoru & remove, const vectoru & keep)
{
	DBG_FAILIF(!keep.empty() && !remove.empty(), ValueError,
		"Please specify one and only one of keep or remove.");

	if (keep.empty() && remove.empty() )
		return;

	vectoru loci;
	if (!keep.empty())
		loci = keep;
	else {
		for (size_t loc = 0; loc < this->totNumLoci(); ++loc)
			// if not removed
			if (find(remove.begin(), remove.end(), loc) == remove.end())
				loci.push_back(loc);
	}

#ifndef OPTIMIZED
	for (size_t i = 0; i < loci.size(); ++i) {
		DBG_FAILIF(loci[i] >= this->totNumLoci(), ValueError,
			"Given loci " + toStr(loci[i]) + " exceed max number of loci.");
		DBG_FAILIF(i > 0 && loci[i] <= loci[i - 1], ValueError,
			"Given loci should be in order.");
	}
#endif
	// adjust order before doing anything
	UINT oldTotNumLoci = totNumLoci();

	// prepare data
	//
	// keep m_popSize;
	// keep m_numSubPop;
	// keep m_subPopSize;
	// keep m_subPopIndex;

	// genotype
	// allocate new genotype and inds
	// new geno structure is in effective now!
	setGenoStructure(removeLociFromGenoStru(vectoru(), loci));

	for (int depth = ancestralDepth(); depth >= 0; --depth) {
		useAncestralPop(depth);

		//
		ULONG newPopGenoSize = genoSize() * m_popSize;
		vectora newGenotype(newPopGenoSize);
		// keep newInds();

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
				for (vectoru::iterator loc = loci.begin();
				     loc != loci.end(); ++loc) {
					*(newPtr++) = oldPtr[*loc];
				}
				oldPtr += oldTotNumLoci;                  // next ploidy
			}
		}
		m_genotype.swap(newGenotype);
	}
	setIndOrdered(true);
}


/** get a new population with selected loci */
population & population::newPopWithPartialLoci(
                                               const vectoru & remove,
                                               const vectoru & keep)
{
	// copy the population over (info is also copied)
	population * pop = new population(*this);

	pop->removeLoci(remove, keep);
	return *pop;
}


void population::rearrangeLoci(const vectoru & newNumLoci, const vectorf & newLociPos)
{
	/// total number of loci can not change
	DBG_FAILIF(std::accumulate(newNumLoci.begin(), newNumLoci.end(), 0U) != totNumLoci(), ValueError,
		"Re-arrange loci must keep the same total number of loci");
	setGenoStructure(ploidy(), newNumLoci.empty() ? numLoci() : newNumLoci,
		sexChrom(), haplodiploid(), newLociPos.empty() ? lociPos() : newLociPos,
		// chromosome names are discarded
		vectorstr(), alleleNames(), lociNames(), maxAllele(), infoFields());
	for (int depth = ancestralDepth(); depth >= 0; --depth) {
		useAncestralPop(depth);

		// now set geno structure
		for (ULONG i = 0; i < m_popSize; ++i)
			// set new geno structure
			m_inds[i].setGenoStruIdx(genoStruIdx());
	}
}


void population::pushAndDiscard(population & rhs, bool force)
{
	// time consuming!
	DBG_ASSERT(rhs.genoStruIdx() == genoStruIdx(), ValueError,
		"Evolution can not continue because the new generation has different \n"
		"genotypic structure. Note that genetypic structure of a population \n"
		"might be changed unexpectedly, e.g. when a sample is drawn from a \n"
		"population.\n");

	DBG_ASSERT(m_genotype.begin() != rhs.m_genotype.begin(), ValueError,
		"Passed population is a reference of current population, swapPop failed.");

	// front -1 pop, -2 pop, .... end
	//
	if (!force && m_ancestralDepth > 0
	    && ancestralDepth() == static_cast<size_t>(m_ancestralDepth) )
		m_ancestralPops.pop_back();

	// save current population
	if (force || m_ancestralDepth != 0) {
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
	m_numSubPop = rhs.m_numSubPop;
	m_subPopSize.swap(rhs.m_subPopSize);
	m_subPopIndex.swap(rhs.m_subPopIndex);
	std::swap(m_vspSplitter, rhs.m_vspSplitter);
	m_genotype.swap(rhs.m_genotype);
	m_info.swap(rhs.m_info);
	m_inds.swap(rhs.m_inds);
	// current population should be working well
	// (with all datamember copied form rhs
	// rhs may not be working well since m_genotype etc
	// may be from ancestral pops
	if (rhs.m_popSize != rhs.m_inds.size()) {
		// keep size if pop size is OK.
		// remove all supopulation structure of rhs
		rhs.m_popSize = rhs.m_inds.size();
		rhs.setSubPopStru(rhs.m_subPopSize);
	}
	validate("Current population after push and discard:");
	rhs.validate("Outside population after push and discard:");
}


// add field
void population::addInfoField(const string field, double init)
{
	DBG_ASSERT(m_info.size() == infoSize() * popSize(), SystemError,
		"Info size is wrong");

	vectorstr newfields;
	UINT os = infoSize();
	UINT idx;
	// if this field exists, return directly
	try {
		idx = infoIdx(field);
		// only needs to initialize
		int oldAncPop = m_curAncestralGen;
		for (UINT anc = 0; anc <= m_ancestralPops.size(); anc++) {
			useAncestralPop(anc);
			for (IndIterator ind = indBegin(); ind.valid(); ++ind)
				ind->setInfo(init, idx);
		}
		useAncestralPop(oldAncPop);
		return;
	} catch (IndexError &) {
		newfields.push_back(field);
	}

	// adjust information size.
	if (!newfields.empty()) {
		setGenoStructure(struAddInfoFields(newfields));
		UINT is = infoSize();
		int oldAncPop = m_curAncestralGen;
		for (UINT anc = 0; anc <= m_ancestralPops.size(); anc++) {
			useAncestralPop(anc);
			vectorinfo newInfo(is * popSize());
			// copy the old stuff in
			InfoIterator ptr = newInfo.begin();
			for (IndIterator ind = indBegin(); ind.valid(); ++ind) {
				copy(ind->infoBegin(), ind->infoBegin() + is - 1, ptr);
				ind->setInfoPtr(ptr);
				ind->setGenoStruIdx(genoStruIdx());
				fill(ind->infoBegin() + os, ind->infoEnd(), init);
				ptr += is;
			}
			m_info.swap(newInfo);
		}
		useAncestralPop(oldAncPop);
	}
	return;
}


void population::addInfoFields(const vectorstr & fields, double init)
{
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
				useAncestralPop(anc);

				for (IndIterator ind = indBegin(); ind.valid(); ++ind)
					ind->setInfo(init, idx);
			}
			useAncestralPop(oldAncPop);
		} catch (IndexError &) {
			newfields.push_back(*it);
		}
	}

	// add these fields
	if (!newfields.empty()) {
		setGenoStructure(struAddInfoFields(newfields));

		// adjust information size.
		UINT is = infoSize();
		int oldAncPop = m_curAncestralGen;
		for (UINT anc = 0; anc <= m_ancestralPops.size(); anc++) {
			useAncestralPop(anc);
			vectorinfo newInfo(is * popSize(), 0.);
			// copy the old stuff in
			InfoIterator ptr = newInfo.begin();
			for (IndIterator ind = indBegin(); ind.valid(); ++ind) {
				copy(ind->infoBegin(), ind->infoBegin() + os, ptr);
				ind->setInfoPtr(ptr);
				ind->setGenoStruIdx(genoStruIdx());
				fill(ind->infoBegin() + os, ind->infoEnd(), init);
				ptr += is;
			}
			m_info.swap(newInfo);
		}
		useAncestralPop(oldAncPop);
	}
}


/// set fields
void population::setInfoFields(const vectorstr & fields, double init)
{
	setGenoStructure(struSetInfoFields(fields));
	// reset info vector
	int oldAncPop = m_curAncestralGen;
	UINT is = infoSize();
	for (UINT anc = 0; anc <= m_ancestralPops.size(); anc++) {
		useAncestralPop(anc);
		vectorinfo newInfo(is * popSize(), init);
		InfoIterator ptr = newInfo.begin();
		for (IndIterator ind = indBegin(); ind.valid(); ++ind, ptr += is) {
			ind->setInfoPtr(ptr);
			ind->setGenoStruIdx(genoStruIdx());
		}
		m_info.swap(newInfo);
	}
	useAncestralPop(oldAncPop);
}


void population::locateRelatives(RelativeType type, const vectorstr & infoFields,
	const vectorstr & parentFields)
{
	if (type == REL_Self) {
		DBG_ASSERT(infoFields.size() == 1, ValueError,
			"Please provide one information field to store REL_Self individuals");
		UINT fieldIdx = infoIdx(infoFields[0]);
			
		for (size_t ans = 0; ans <= ancestralDepth(); ++ans) {
			useAncestralPop(ans);
			for (size_t idx = 0; idx < popSize(); ++idx)
				ind(idx).setInfo(idx, fieldIdx);
		}
		useAncestralPop(0);
	} else if (type == REL_Spouse) {
		DBG_ASSERT(parentFields.size() == 2, ValueError,
			"This relative only exists when there are two parents for each indidivual");

		DBG_ASSERT(infoFields.size() >= 1, ValueError,
			"Please provide at least one information field to store REL_Self individuals");

		UINT maxSpouse = infoFields.size();

		vectori spouseIdx(maxSpouse);
		for (size_t i = 0; i < maxSpouse; ++i) {
			spouseIdx[i] = infoIdx(infoFields[i]);
			// clear these fields for the last generation
			for (IndInfoIterator ptr = infoBegin(spouseIdx[i]);
				ptr != infoEnd(spouseIdx[i]); ++ptr)
				*ptr = static_cast<InfoType>(-1);
		}

		// start from the parental generation
		for (unsigned ans = 1; ans <= ancestralDepth(); ++ans) {
			vectoru numSpouse;
			// go to offspring generation
			useAncestralPop(ans - 1);
			vectorf father = indInfo(parentFields[0]);
			vectorf mother = indInfo(parentFields[1]);
			//
			useAncestralPop(ans);
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
					bool found = false;
					for (size_t s = 0; s < numSpouse[p]; ++s)
						if (ind(p).info(spouseIdx[s]) == m) {
							found = true;
							break;
						}
					if (!found) {
						ind(p).setInfo(m, spouseIdx[numSpouse[p]]);
						++numSpouse[p];
					}
				}
				if (numSpouse[m] < maxSpouse) {
					bool found = false;
					for (size_t s = 0; s < numSpouse[m]; ++s)
						if (ind(p).info(spouseIdx[s]) == p) {
							found = true;
							break;
						}
					if (!found) {
						ind(m).setInfo(p, spouseIdx[numSpouse[m]]);
						++numSpouse[m];
					}
				}				                                                                           // idx
			}                                                                                               // ancestal generations
			// set the rest of the field to -1
			for (size_t idx = 0; idx < popSize(); ++idx) {
				for (size_t ns = numSpouse[idx]; ns < maxSpouse; ++ns)
					ind(idx).setInfo(-1, spouseIdx[ns]);
			}
		}
		useAncestralPop(0);
	} else if (type == REL_Offspring) {
		DBG_ASSERT(infoFields.size() >= 1, ValueError,
			"Please provide at least one information field to store offspring");

		UINT maxOffspring = infoFields.size();

		vectori offspringIdx(maxOffspring);
		for (size_t i = 0; i < maxOffspring; ++i) {
			offspringIdx[i] = infoIdx(infoFields[i]);
			// clear these fields for the last generation
			for (IndInfoIterator ptr = infoBegin(offspringIdx[i]);
				ptr != infoEnd(offspringIdx[i]); ++ptr)
				*ptr = static_cast<InfoType>(-1);
		}

		// start from the parental generation
		for (unsigned ans = 1; ans <= ancestralDepth(); ++ans) {
			vectoru numOffspring;
			// for each type of parental relationship
			for (vectorstr::const_iterator field = parentFields.begin();
				field != parentFields.end(); ++field) {
				useAncestralPop(ans - 1);
				vectorf parent = indInfo(*field);
				//
				useAncestralPop(ans);
				if (numOffspring.empty())
					numOffspring.resize(popSize(), 0);
				//
				for (size_t idx = 0; idx < parent.size(); ++idx) {
					DBG_ASSERT(fcmp_ne(parent[idx], -1), ValueError, "Invalid parental index (-1)");
					ULONG p = static_cast<ULONG>(parent[idx]);
					DBG_ASSERT(p < popSize(), IndexError, "Parental index out of range of 0 ~ " + toStr(popSize() - 1));
					if (numOffspring[p] < maxOffspring) {
						ind(p).setInfo(idx, offspringIdx[numOffspring[p]]);
						++numOffspring[p];
					}
				}                                                                                           // idx
			}                                                                                               // ancestal generations
			// set the rest of the field to -1
			for (size_t idx = 0; idx < popSize(); ++idx) {
				for (size_t no = numOffspring[idx]; no < maxOffspring; ++no)
					ind(idx).setInfo(-1, offspringIdx[no]);
			}
		}
		useAncestralPop(0);
	} else if (type == REL_Sibling ) {
	}
}


// set ancestral depth, can be -1
void population::setAncestralDepth(int depth)
{
	// just to make sure.
	useAncestralPop(0);
	//
	if (depth >= 0 && m_ancestralPops.size() > static_cast<size_t>(depth)) {
		int numRemove = m_ancestralPops.size() - depth;
		while (numRemove-- > 0)
			m_ancestralPops.pop_back();
	}
	DBG_ASSERT(depth < 0 || m_ancestralPops.size() <= static_cast<size_t>(depth), SystemError,
		"Failed to change ancestral Depth");

	m_ancestralDepth = depth;
}


void population::useAncestralPop(UINT idx)
{
	if (m_curAncestralGen >= 0 && idx == static_cast<UINT>(m_curAncestralGen))
		return;

	DBG_DO(DBG_POPULATION, cout << "Use ancestralPop: " << idx <<
		"Curidx: " << m_curAncestralGen << endl);

	if (idx == 0 || m_curAncestralGen != 0) {         // recover pop.
		popData & pd = m_ancestralPops[ m_curAncestralGen - 1];
		pd.swap(*this);
		m_curAncestralGen = 0;
		if (idx == 0) {                                               // restore key parameters from data
			m_popSize = m_inds.size();
			setSubPopStru(m_subPopSize);
			return;
		}
	}

	// now m_curAncestralGen is zero.
	DBG_ASSERT(idx <= m_ancestralPops.size(),
		ValueError, "Ancestry population " + toStr(idx) + " does not exist.");

	// now idx should be at least 1
	m_curAncestralGen = idx;
	// swap  1 ==> 0, 2 ==> 1

	popData & pd = m_ancestralPops[m_curAncestralGen - 1];
	pd.swap(*this);
	m_popSize = m_inds.size();
	setSubPopStru(m_subPopSize);
}


void population::savePopulation(const string & filename, const string & format, bool compress) const
{
	io::filtering_ostream ofs;
	// get file extension
	string ext = fileExtension(filename);

#ifndef DISABLE_COMPRESSION
	if (compress || (ext.size() > 3 && ext.substr(ext.size() - 3, 3) == ".gz"))
		ofs.push(io::gzip_compressor());
#endif
	ofs.push(io::file_sink(filename));

	if (!ofs)
		throw ValueError("Can not open file " + filename);

	if (format == "text" || (format == "auto" && (ext == "txt" || ext == "txt.gz"))) {
		boost::archive::text_oarchive oa(ofs);
		oa << *this;
	} else if (format == "xml" || (format == "auto" && (ext == "xml" || ext == "xml.gz"))) {
		boost::archive::xml_oarchive oa(ofs);
		oa << boost::serialization::make_nvp("population", *this);
	} else if (format == "bin" || (format == "auto" && (ext == "bin" || ext == "bin.gz"))) {
		boost::archive::binary_oarchive oa(ofs);
		oa << *this;
	} else
		throw ValueError("Wrong format type. Use one of text, xml, bin or appropriate extension txt, xml or bin");
}


void population::loadPopulation(const string & filename, const string & format)
{
	io::filtering_istream ifs;
	bool gzipped = isGzipped(filename);

	if (gzipped)
#ifdef DISABLE_COMPRESSION
		throw ValueError("This version of simuPOP can not handle compressed file");
#else
		ifs.push(io::gzip_decompressor());
#endif
	ifs.push(io::file_source(filename));
	// do not have to test again.
	if (!ifs)
		throw ValueError("Can not open file " + filename);

	// get file extension
	string ext = fileExtension(filename);

	// try to load the file, according to file extension.
	try {
		if (format == "text" || (format == "auto" && (ext == "txt" || ext == "txt.gz") )) {
			boost::archive::text_iarchive ia(ifs);
			ia >> *this;
		} else if (format == "xml" || (format == "auto" && (ext == "xml" || ext == "xml.gz") )) {
			boost::archive::xml_iarchive ia(ifs);
			ia >> boost::serialization::make_nvp("population", *this);
		} else if (format == "bin" || (format == "auto" && (ext == "bin" || ext == "bin.gz") )) {
			boost::archive::binary_iarchive ia(ifs);
			ia >> *this;
		} else                                                                              // need special handling
			throw;
	} catch (...) {                                                                         // if any error happens, or can not determine format, try different methods
		// first close the file handle.

		DBG_DO(DBG_POPULATION,
			cout << "Can not determine file type, or file type is wrong. Trying different ways." << endl);

		// open a fresh ifstream
		io::filtering_istream ifbin;
		if (gzipped)
			ifbin.push(io::gzip_decompressor());
		ifbin.push(io::file_source(filename));

		// try to load the file using different iarchives.
		try                                                                               // binary?
		{
			boost::archive::binary_iarchive ia(ifbin);
			ia >> *this;
		} catch (...) {                                                             // not binary, text?
			io::filtering_istream iftxt;
			if (gzipped)
				iftxt.push(io::gzip_decompressor());
			iftxt.push(io::file_source(filename));
			try {
				boost::archive::text_iarchive ia(iftxt);
				ia >> *this;
			} catch (...) {                                                     // then xml?
				io::filtering_istream ifxml;
				if (gzipped)
					ifxml.push(io::gzip_decompressor());
				ifxml.push(io::file_source(filename));
				try {
					boost::archive::xml_iarchive ia(ifxml);
					ia >> boost::serialization::make_nvp("population", *this);
				} catch (...) {
					throw ValueError("Failed to load population " + filename + " in " + format + " format.\n");
				}
			}                                                                                       // try xml
		}                                                                                           // try text
	}                                                                                               // try bin
}


PyObject * population::vars(int subPop)
{
	if (subPop < 0) {
		Py_INCREF(m_vars.dict());
		return m_vars.dict();
	} else {
		DBG_ASSERT(static_cast<UINT>(subPop) < numSubPop(),
			IndexError, "Subpop index out of range of 0 ~ " + toStr(numSubPop() - 1) );

		DBG_ASSERT(hasVar("subPop"), ValueError,
			"subPop statistics does not exist yet.");

		PyObject * spObj = getVar("subPop");
		spObj = PyList_GetItem(spObj, subPop);

		DBG_ASSERT(spObj != NULL, SystemError,
			"Something is wrong about the length of subPop list. ");

		Py_INCREF(spObj);
		return spObj;
	}
}


// CPPONLY
// The same as vars(), but without increasing
// reference count.
PyObject * population::dict(int subPop)
{
	if (subPop < 0)
		return m_vars.dict();
	else {
		DBG_ASSERT(static_cast<UINT>(subPop) < numSubPop(),
			IndexError, "Subpop index out of range of 0 ~ " + toStr(numSubPop() - 1) );

		DBG_ASSERT(hasVar("subPop"), ValueError,
			"subPop statistics does not exist yet.");

		PyObject * spObj = getVar("subPop");
		spObj = PyList_GetItem(spObj, subPop);

		DBG_ASSERT(spObj != NULL, SystemError,
			"Something is wrong about the length of subPop list. ");

		return spObj;
	}
}


/// CPPONLY
void population::sortIndividuals(bool infoOnly)
{
	if (indOrdered())
		return;

	if (infoOnly) {
		DBG_DO(DBG_POPULATION, cout << "Adjust info position " << endl);
		UINT is = infoSize();
		vectorinfo tmpInfo(m_popSize * is);
		vectorinfo::iterator infoPtr = tmpInfo.begin();

		for (IndIterator ind = indBegin(); ind.valid(); ++ind) {
			copy(ind->infoBegin(), ind->infoEnd(), infoPtr);
			ind->setInfoPtr(infoPtr);
			infoPtr += is;
		}
		m_info.swap(tmpInfo);
	} else {
		DBG_DO(DBG_POPULATION, cout << "Adjust geno and info position " << endl);

		size_t sz = genoSize();
		UINT is = infoSize();
		vectora tmpGenotype(m_popSize * genoSize());
		vectorinfo tmpInfo(m_popSize * infoSize());
		vectora::iterator it = tmpGenotype.begin();
		vectorinfo::iterator infoPtr = tmpInfo.begin();

		for (IndIterator ind = indBegin(); ind.valid(); ++ind) {
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
		m_genotype.swap(tmpGenotype);
		m_info.swap(tmpInfo);
	}
	setIndOrdered(true);
}


population & LoadPopulation(const string & file, const string & format)
{
#ifndef _NO_SERIALIZATION_
	population * p = new population();
	p->loadPopulation(file, format);
	return *p;
#else
	cout << "This feature is not supported in this platform" << endl;
	return *new population(1);
#endif
}


vectorf testGetinfoFromInd(population & pop)
{
	vectorf a(pop.popSize());
	size_t i = 0;

	for (IndIterator ind = pop.indBegin(); ind.valid(); ++ind)
		a[i++] = ind->info(0);
	return a;
}


vectorf testGetinfoFromPop(population & pop, bool order)
{
	vectorf a(pop.popSize());
	size_t i = 0;

	IndInfoIterator it = pop.infoBegin(0);
	IndInfoIterator it_end = pop.infoEnd(0);

	for (; it != it_end; ++it)
		a[i++] = *it;
	return a;
}


}


