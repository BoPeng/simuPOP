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


namespace simuPOP {

population::population(const vectorlu & size,
	float ploidy,
	const vectoru & loci,
	const vectoru & chromTypes,
	const vectorf & lociPos,
	int ancestralGens,
	const vectorstr & chromNames,
	const vectorstr & alleleNames,
	const vectorstr & lociNames,
	const vectorstr & infoFields)
	:
	GenoStruTrait(),
	m_popSize(0),
	m_subPopSize(size),
	m_subPopIndex(size.size() + 1),
	m_vspSplitter(NULL),
	m_genotype(0),                                                                          // resize later
	m_info(0),
	m_inds(0),                                                                              // default constructor will be called.
	m_ancestralGens(ancestralGens),
	m_vars(NULL, true),                                                                     // invalid shared variables initially
	m_ancestralPops(0),                                                                     // no history first
	m_rep(-1),
	m_gen(0),
	m_curAncestralGen(0),
	m_indOrdered(true),
	m_selectionFlags()
{
	DBG_DO(DBG_POPULATION, cout << "Constructor of population is called\n");

	DBG_FAILIF(m_subPopSize.size() > MaxSubPopID, ValueError,
		"Number of subpopulations exceed maximum allowed subpopulation numbers");

	// get a GenoStructure with parameters. GenoStructure may be shared by some populations
	// a whole set of functions ploidy() etc in GenoStruTriat can be used after this step.
	DBG_FAILIF(static_cast<UINT>(ploidy) * 1.0 != ploidy && fcmp_ne(ploidy, Haplodiploid),
		ValueError, "Only integer ploidy number or Haplodiploid can be specified");

	setGenoStructure(fcmp_eq(ploidy, Haplodiploid) ? 2 : static_cast<UINT>(ploidy),
		loci, chromTypes, fcmp_eq(ploidy, Haplodiploid), lociPos, chromNames, alleleNames,
		lociNames, infoFields);

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
}


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
	m_subPopSize(rhs.m_subPopSize),
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
		m_ancestralGens = 0;
		m_ancestralPops.clear();
	}

	// copy virtual subpop splitters
	setVirtualSplitter(rhs.virtualSplitter());

	// set local variable
	setRep(-1);
}


void population::popData::swap(population & pop)
{
	pop.m_subPopSize.swap(m_subPopSize);
	pop.m_genotype.swap(m_genotype);
	pop.m_info.swap(m_info);
	pop.m_inds.swap(m_inds);
	std::swap(pop.m_indOrdered, m_indOrdered);
}


population * population::clone(int ancGen) const
{
	population * p = new population(*this);
	int oldDepth = m_ancestralGens;

	if (ancGen >= 0)
		// try to remove excessive ancestra generations.
		p->setAncestralDepth(ancGen);
	p->setAncestralDepth(oldDepth);
	return p;
}


string population::virtualSubPopName(vspID subPop) const
{
	DBG_ASSERT(hasVirtualSubPop(), ValueError,
		"No virtual subpopulation is defined for this population.");
	// if a single number is given, it will be passed as (sp, None),
	// but we will treat sp as vsp here.
	if (!subPop.isVirtual())
		return m_vspSplitter->name(subPop.subPop());
	else
		return m_vspSplitter->name(subPop.virtualSubPop());
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
                                       IterationType type)
{
	CHECKRANGESUBPOP(subPop);
	DBG_ASSERT(virtualSubPop != InvalidSubPopID, ValueError, "Given virtual subpopulation ID is wrong");
	DBG_ASSERT(hasVirtualSubPop(), ValueError,
		"Subpopulation " + toStr(subPop) + " has no virtual subpopulations");
	m_vspSplitter->activate(*this, subPop, virtualSubPop, type);
	DBG_ASSERT(type != VisibleInds ||
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


individual & population::ancestor(ULONG idx, UINT gen)
{
	DBG_FAILIF(gen > m_ancestralPops.size(), IndexError,
		"Ancestray generation " + toStr(gen) + " does not exist");
	if (gen == m_curAncestralGen)
		return this->ind(idx);
	UINT genIdx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
	DBG_FAILIF(idx > m_ancestralPops[genIdx].m_inds.size(),
		IndexError, "Individual index out of range");
	return m_ancestralPops[genIdx].m_inds[idx];
}


const individual & population::ancestor(ULONG idx, UINT gen) const
{
	DBG_FAILIF(gen > m_ancestralPops.size(), IndexError,
		"Ancestray generation " + toStr(gen) + " does not exist");
	if (gen == m_curAncestralGen)
		return this->ind(idx);
	UINT genIdx = gen == 0 ? m_curAncestralGen - 1 : gen - 1;
	DBG_FAILIF(idx > m_ancestralPops[genIdx].m_inds.size(),
		IndexError, "Individual index out of range");
	return m_ancestralPops[genIdx].m_inds[idx];
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


PyObject * population::genotype()
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	sortIndividuals();
	// directly expose values. Do not copy data over.
	return Allele_Vec_As_NumArray(m_genotype.begin(), m_genotype.end());
}


PyObject * population::genotype(SubPopID subPop)
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	CHECKRANGESUBPOP(subPop);
	sortIndividuals();
	// directly expose values. Do not copy data over.
	return Allele_Vec_As_NumArray(genoBegin(subPop, true), genoEnd(subPop, true));
}


void population::setGenotype(vectora geno)
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	sortIndividuals();
	GenoIterator ptr = m_genotype.begin();
	ULONG sz = geno.size();
	for (ULONG i = 0; i < popSize() * genoSize(); ++i)
		*(ptr++) = geno[i % sz];
}


void population::setGenotype(vectora geno, SubPopID subPop)
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	CHECKRANGESUBPOP(subPop);
	sortIndividuals();

	GenoIterator ptr = genoBegin(subPop, true);
	ULONG sz = geno.size();
	for (ULONG i = 0; i < subPopSize(subPop) * genoSize(); ++i)
		*(ptr++) = geno[i % sz];
}


void population::setIndSubPopID(const vectori & id, bool ancestralPops)
{
	UINT oldGen = curAncestralGen();
	size_t sz = id.size();

	for (UINT anc = 0; anc <= ancestralGens(); ++anc) {
		if (!ancestralPops && anc != oldGen)
			continue;
		useAncestralGen(anc);
		for (ULONG it = 0; it < m_popSize; ++it)
			ind(it).setSubPopID(static_cast<SubPopID>(id[it % sz]));
	}
	useAncestralGen(oldGen);
}


void population::setIndSubPopIDWithID(bool ancestralPops)
{
	UINT oldGen = curAncestralGen();

	for (UINT anc = 0; anc <= ancestralGens(); ++anc) {
		if (!ancestralPops && anc != oldGen)
			continue;
		useAncestralGen(anc);
		for (UINT i = 0, iEnd = numSubPop(); i < iEnd;  ++i)
			for (IndIterator it = indBegin(i); it.valid();  ++it)
				it->setSubPopID(i);
	}
	useAncestralGen(oldGen);
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

	if (newSubPopSizes.empty())
		m_subPopSize = vectorlu(1, 0);
	else
		m_subPopSize = newSubPopSizes;
	m_subPopIndex.resize(numSubPop() + 1);

	// build subPop index
	UINT i = 1;
	for (m_subPopIndex[0] = 0; i <= numSubPop(); ++i)
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
		m_subPopSize.resize(1, 0);
		m_subPopIndex.resize(2);
	} else {
		// reset indexes etc.
		UINT numSubPop = static_cast<UINT>(m_inds.back().subPopID()) + 1;
		m_subPopSize.resize(numSubPop);
		m_subPopIndex.resize(numSubPop + 1);

		// check subpop size
		fill(m_subPopSize.begin(), m_subPopSize.end(), 0);
		for (IndIterator it = indBegin(); it.valid();  ++it)
			m_subPopSize[ static_cast<UINT>(it->subPopID()) ]++;
	}
	// rebuild index
	size_t i = 1;
	for (m_subPopIndex[0] = 0; i <= numSubPop(); ++i)
		m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];
}


void population::splitSubPop(UINT subPop, vectorf sizes)
{
	if (sizes.size() <= 1)
		return;

	double all = accumulate(sizes.begin(), sizes.end(), 0.);
	DBG_ASSERT(static_cast<ULONG>(all) == subPopSize(subPop) || fcmp_eq(all, 1.),
		ValueError,
		"Sum of parameter sizes should be 1 (proportions) or the size of subpopulation subPop.");


	ULONG spSize = subPopSize(subPop);
	vectorlu newSizes(sizes.size());
	for (size_t i = 0; i < sizes.size() - 1; ++i) {
		if (fcmp_eq(all, 1.))
			newSizes[i] = static_cast<ULONG>(floor(spSize * sizes[i]));
		else
			newSizes[i] = static_cast<ULONG>(sizes[i]);
	}
	// to avoid round off problem, calculate the last subpopulation
	newSizes[sizes.size() - 1] = spSize - accumulate(newSizes.begin(), newSizes.end() - 1, 0L);

	vectorlu subPopSizes;
	for (size_t sp = 0; sp < numSubPop(); ++sp)
		if (sp != subPop)
			subPopSizes.push_back(subPopSize(sp));
		else
			subPopSizes.insert(subPopSizes.end(), newSizes.begin(), newSizes.end());
	setSubPopStru(subPopSizes);
	return;
}


void population::removeEmptySubPops()
{
	// if remove empty subpops
	UINT newSPNum = numSubPop();
	vectorlu newSPSize;

	for (size_t sp = 0; sp < numSubPop(); ++sp) {
		if (m_subPopSize[sp] == 0)
			newSPNum--;
		else
			newSPSize.push_back(m_subPopSize[sp]);
	}
	m_subPopSize.swap(newSPSize);
	m_subPopIndex.resize(numSubPop() + 1);
	// rebuild index
	size_t i = 1;
	for (m_subPopIndex[0] = 0; i <= numSubPop(); ++i)
		m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];
}


void population::removeSubPops(const vectoru & subPops)
{
#ifndef OPTIMIZED
	// check if subPops are valid
	for (vectoru::const_iterator sp = subPops.begin(); sp < subPops.end(); ++sp) {
		DBG_WARNING(*sp >= numSubPop(), "Subpopulation " + toStr(*sp) + " does not exist.");
	}
#endif
	sortIndividuals();
	vectorlu new_size;

	UINT step = genoSize();
	UINT infoStep = infoSize();
	vector<individual>::iterator oldInd = m_inds.begin();
	vector<individual>::iterator newInd = m_inds.begin();
	GenoIterator oldPtr = m_genotype.begin();
	InfoIterator oldInfoPtr = m_info.begin();
	GenoIterator newPtr = m_genotype.begin();
	InfoIterator newInfoPtr = m_info.begin();

	for (size_t sp = 0; sp < numSubPop(); ++sp) {
		ULONG spSize = subPopSize(sp);
		if (find(subPops.begin(), subPops.end(), sp) == subPops.end()) {
			new_size.push_back(spSize);
			// do not remove.
			if (oldPtr != newPtr) {
				copy(oldInd, oldInd + spSize, newInd);
				copy(oldPtr, oldPtr + step * spSize, newPtr);
				copy(oldInfoPtr, oldInfoPtr + infoStep * spSize, newInfoPtr);
			}
			newInd += spSize;
			newPtr += step * spSize;
			newPtr += infoStep * spSize;
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
	setSubPopStru(new_size);
	//
	GenoIterator ptr = m_genotype.begin();
	InfoIterator infoPtr = m_info.begin();
	for (ULONG i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		m_inds[i].setGenoPtr(ptr);
		m_inds[i].setInfoPtr(infoPtr);
	}
}


void population::removeIndividuals(const vectoru & inds)
{
	sortIndividuals();
	vectorlu new_size = m_subPopSize;

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
	for (; it != it_end; ++it)
		removed[*it] = true;
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
				newPtr += infoStep;
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
	setSubPopStru(new_size);
	//
	GenoIterator ptr = m_genotype.begin();
	InfoIterator infoPtr = m_info.begin();
	for (ULONG i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		m_inds[i].setGenoPtr(ptr);
		m_inds[i].setInfoPtr(infoPtr);
	}
}


void population::mergeSubPops(const vectoru & subPops)
{
	// merge all subpopulations
	if (subPops.empty()) {
		// [ popSize() ]
		vectorlu sz(1, popSize());
		setSubPopStru(sz);
		return;
	}
	if (subPops.size() == 1)
		return;

	// are they in order?
	bool consecutive = true;
	vectoru sps = subPops;
	sort(sps.begin(), sps.end());
	for (size_t i = 1; i < sps.size(); ++i)
		if (subPops[i] != subPops[i - 1] + 1) {
			consecutive = false;
			break;
		}
	// new subpopulation sizes
	vectorlu new_size;
	for (UINT sp = 0; sp < numSubPop(); ++sp) {
		if (find(subPops.begin(), subPops.end(), sp) != subPops.end()) {
			if (new_size.size() <= sps[0])
				new_size.push_back(subPopSize(sp));
			else
				new_size[sps[0]] += subPopSize(sp);
		} else
			new_size.push_back(subPopSize(sp));
	}
	// if consecutive, no need to move anyone
	if (consecutive) {
		setSubPopStru(new_size);
		return;
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
	setSubPopStru(new_size);
	//
	GenoIterator ptr = m_genotype.begin();
	InfoIterator infoPtr = m_info.begin();
	for (ULONG i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
		m_inds[i].setGenoPtr(ptr);
		m_inds[i].setInfoPtr(infoPtr);
	}
}


void population::addChromFromPop(const population & pop)
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


void population::addIndFromPop(const population & pop)
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
}


void population::addLociFromPop(const population & pop)
{
	DBG_FAILIF(ancestralGens() != pop.ancestralGens(), ValueError,
		"Can not add chromosomes from a population with different number of ancestral generations");

	vectorstr lociNames1 = lociNames();
	vectorstr lociNames2 = pop.lociNames();
	// obtain new genotype structure and set it
	setGenoStructure(gsAddLociFromStru(pop.genoStruIdx()));
	vectoru indexes1 = lociByNames(lociNames1);
	vectoru indexes2 = lociByNames(lociNames2);

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
		UINT size1 = lociNames1.size();
		UINT size2 = lociNames2.size();
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
                          const string & chromName, UINT chromType)
{
	DBG_ASSERT(lociNames.empty() || lociPos.size() == lociNames.size(), ValueError,
		"Please specifiy locus name for all inserted loci.");

	size_t oldNumLoci = totNumLoci();
	// obtain new genotype structure and set it
	setGenoStructure(gsAddChrom(lociPos, lociNames, chromName, chromType));

	DBG_FAILIF(totNumLoci() - oldNumLoci == lociPos.size(), SystemError,
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


vectoru population::addLoci(const vectoru & chrom, const vectorf & pos,
                            const vectorstr & names)
{
	DBG_ASSERT(chrom.size() == pos.size(), ValueError,
		"Chromosome and position lists should have the same length");
	DBG_ASSERT(names.empty() || pos.size() == names.size(), ValueError,
		"Please specifiy locus name for all inserted loci.");

	vectoru newIndex;
	vectoru loci(totNumLoci());
	// obtain new genotype structure and set it
	setGenoStructure(gsAddLoci(chrom, pos, names, newIndex));
	// use loci to keep the position of old loci in the new structure
	for (size_t i = 0, j = 0; j < totNumLoci(); ++j) {
		// i is the index to loci before insertion.
		// j is the index to loci after insertion.
		if (find(newIndex.begin(), newIndex.end(), i) == newIndex.end()) {
			loci[i] = j;
			++i;
		}
	}

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
	for (m_subPopIndex[0] = 0; idx <= numSubPop(); ++idx)
		m_subPopIndex[idx] = m_subPopIndex[idx - 1] + m_subPopSize[idx - 1];
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
	population * pop = new population(sz, ploidy(), numLoci(), chromTypes(), lociPos(), 0,
		chromNames(), alleleNames(), lociNames(), infoFields());
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
population & population::newPopByIndID(int ancGen,
                                       const vectori & id, bool removeEmptySubPops)
{
	UINT topGen;

	if (ancGen < 0 || static_cast<UINT>(ancGen) >= ancestralGens())
		topGen = ancestralGens();
	else
		topGen = ancGen;
	// go to the oldest generation
	useAncestralGen(topGen);
	population & ret = newPopByIndIDPerGen(id, removeEmptySubPops);
	// prepare for push and discard
	ret.setAncestralDepth(topGen);
	if (topGen > 0) {
		for (int depth = topGen - 1; depth >= 0; --depth) {
			useAncestralGen(depth);
			ret.pushAndDiscard(newPopByIndIDPerGen(id, removeEmptySubPops));
		}
	}
	return ret;
}


void population::removeLoci(const vectoru & loci, const vectoru & keep)
{
	DBG_FAILIF(!loci.empty() && !keep.empty(), ValueError,
		"Please specify only one of parameters loci and keep");

	if (loci.empty() && keep.empty())
		return;

	vectoru kept = keep;
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
					*(newPtr++) = oldPtr[*loc];
				oldPtr += oldTotNumLoci;
			}
		}
		m_genotype.swap(newGenotype);
	}
	setIndOrdered(true);
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
	if (!force && m_ancestralGens > 0
	    && ancestralGens() == static_cast<size_t>(m_ancestralGens) )
		m_ancestralPops.pop_back();

	// save current population
	if (force || m_ancestralGens != 0) {
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
			useAncestralGen(anc);
			for (IndIterator ind = indBegin(); ind.valid(); ++ind)
				ind->setInfo(init, idx);
		}
		useAncestralGen(oldAncPop);
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
			useAncestralGen(anc);
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
		useAncestralGen(oldAncPop);
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
				useAncestralGen(anc);

				for (IndIterator ind = indBegin(); ind.valid(); ++ind)
					ind->setInfo(init, idx);
			}
			useAncestralGen(oldAncPop);
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
			useAncestralGen(anc);
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
		useAncestralGen(oldAncPop);
	}
}


void population::setInfoFields(const vectorstr & fields, double init)
{
	setGenoStructure(struSetInfoFields(fields));
	// reset info vector
	int oldAncPop = m_curAncestralGen;
	UINT is = infoSize();
	for (UINT anc = 0; anc <= m_ancestralPops.size(); anc++) {
		useAncestralGen(anc);
		vectorinfo newInfo(is * popSize(), init);
		InfoIterator ptr = newInfo.begin();
		for (IndIterator ind = indBegin(); ind.valid(); ++ind, ptr += is) {
			ind->setInfoPtr(ptr);
			ind->setGenoStruIdx(genoStruIdx());
		}
		m_info.swap(newInfo);
	}
	useAncestralGen(oldAncPop);
}


void population::setIndInfo(const vectorinfo & values, UINT idx)
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	CHECKRANGEINFO(idx);
	size_t valueSize = values.size();
	IndInfoIterator ptr = infoBegin(idx);
	for (size_t i = 0; ptr != infoEnd(idx); ++ptr, ++i)
		*ptr = static_cast<InfoType>(values[i % valueSize]);
}


void population::setIndInfo(const vectorinfo & values, UINT idx, vspID subPop)
{
	DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
		"This operation is not allowed when there is an activated virtual subpopulation");

	CHECKRANGEINFO(idx);
	size_t valueSize = values.size();
	IndInfoIterator ptr = infoBegin(idx, subPop);
	for (size_t i = 0; ptr != infoEnd(idx, subPop); ++ptr, ++i)
		*ptr = static_cast<InfoType>(values[i % valueSize]);
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


void population::sortIndividuals(bool infoOnly)
{
	if (indOrdered())
		return;

	if (infoOnly) {
		DBG_DO(DBG_POPULATION, cout << "Adjust info position " << endl);
		UINT is = infoSize();
		if (is == 0) {
			setIndOrdered(true);
			return;
		}
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


void population::scramble()
{
	UINT step = genoSize();
	UINT infoStep = infoSize();

	vectorlu pointers(m_popSize);

	for (ULONG i = 0; i < m_popSize; ++i)
		pointers[i] = i;
	random_shuffle(pointers.begin(), pointers.end());

	vectora newGenotype(m_popSize * genoSize());
	vectorinfo newInfo(m_popSize * infoSize());
	vectora::iterator newGenoPtr = newGenotype.begin();
	vectorinfo::iterator newInfoPtr = newInfo.begin();

	GenoIterator ptr = m_genotype.begin();
	InfoIterator infoPtr = m_info.begin();
	for (ULONG i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
#ifdef BINARYALLELE
		copyGenotype(m_inds[i].genoBegin(), newGenoPtr + pointers[i] * step, step);
#else
		copy(m_inds[i].genoBegin(), m_inds[i].genoEnd(), newGenoPtr + pointers[i] * step);
#endif
		m_inds[i].setGenoPtr(newGenoPtr + pointers[i] * step);

		copy(m_inds[i].infoBegin(), m_inds[i].infoEnd(), newInfoPtr + pointers[i] * infoStep);
		m_inds[i].setInfoPtr(newInfoPtr + pointers[i] * infoStep);
	}
	m_genotype.swap(newGenotype);
	m_info.swap(newInfo);
	setIndOrdered(false);
}


population & LoadPopulation(const string & file)
{
	population * p = new population();

	p->load(file);
	return *p;
}


}


