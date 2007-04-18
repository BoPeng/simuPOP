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

// for file compression
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#ifdef SIMUMPI
#include <boost/parallel/mpi.hpp>
namespace mpi = boost::parallel::mpi;

// define action codes
#include "slave.h"
#endif

namespace io = boost::iostreams;

namespace simuPOP
{
	individual& individualIterator::next()
	{
		if(m_index == m_end)
			throw StopIteration("");
		else
			return m_population->ind(m_index++);
	}

	population::population( ULONG size,
		UINT ploidy,
		const vectoru& loci,
		bool sexChrom,
		const vectorf& lociPos,
		const vectorlu& subPop,
		int ancestralDepth,
		const vectorstr& alleleNames,
		const vectorstr& lociNames,
		UINT maxAllele,
		const vectorstr& infoFields,
		const vectori& chromMap)
		:
	GenoStruTrait(),
		m_popSize(size),
		m_numSubPop(subPop.size()),
		m_subPopSize(subPop),
		m_subPopIndex(subPop.size()+1),
		m_genotype(0),							  // resize later
		m_info(0),
		m_inds(0),								  // default constructor will be called.
		m_ancestralDepth(ancestralDepth),
		m_vars(NULL, true),						  // invalid shared variables initially
		m_ancestralPops(0),						  // no history first
		m_rep(-1),
		m_grp(-1),
		m_gen(0),
		m_curAncestralPop(0),
		m_shallowCopied(false),
		m_infoOrdered(true)
	{
		DBG_FAILIF(maxAllele > MaxAllele, ValueError,
			"maxAllele is bigger than maximum allowed allele state of this library (" + toStr(MaxAllele) +
			")\nPlease use simuOpt.setOptions(alleleType='long') to use the long allele version of simuPOP.");

		DBG_FAILIF(maxAllele == 0, ValueError,
			"maxAllele should be at least 1 (0,1 two states). ");

		DBG_DO( DBG_POPULATION, cout << "Constructor of population is called\n");

		DBG_FAILIF(m_subPopSize.size() > MaxSubPopID, ValueError,
			"Number of subpopulations exceed maximum allowed subpopulation numbers");

		// if specify subPop but not m_popSize
		if( !subPop.empty() )
		{
			if( size == 0 )
				m_popSize = accumulate(subPop.begin(), subPop.end(), 0UL);
			else
				DBG_ASSERT( m_popSize == accumulate(subPop.begin(), subPop.end(), 0UL),
					ValueError, "If both size and subPop are specified, size should equal to sum(subPop)");
		}

#ifdef SIMUMPI
		UINT nodes = chromMap.empty()?(loci.size()+1):(chromMap.size() + 1);
		if(mpiSize() < nodes)
		{
			if (mpiRank() == 0)
			{
				cout << "Available nodes: " << mpiSize() << endl;
				cout << "Number of chromosomes: " << loci.size() << endl;
				if (!chromMap.empty())
					cout << "Length of chromosome map: " << chromMap.size() << endl;
			}
			throw SystemError("Insufficient nodes. At least " + toStr(nodes) + " (1 + number of chromosomes "
				"or 1 + length of chromMap ) number of nodes are required\n");
		}
		m_popID = uniqueID();
		// create the population on other nodes by sending other nodes the command and parameters
		if (mpiRank() == 0)
		{
			for(size_t node = 1; node < nodes; ++node)
			{
				int action = SLAVE_POPULATION_CREATE;
				mpiComm().send(node, 0, action);
				mpiComm().send(node, 1, m_popID);
				mpiComm().send(node, 2, size);
				mpiComm().send(node, 3, ploidy);
				mpiComm().send(node, 4, loci);
				mpiComm().send(node, 5, sexChrom);
				mpiComm().send(node, 6, lociPos);
				mpiComm().send(node, 7, subPop);
				mpiComm().send(node, 8, ancestralDepth);
				mpiComm().send(node, 9, alleleNames);
				mpiComm().send(node, 10, lociNames);
				mpiComm().send(node, 11, maxAllele);
				mpiComm().send(node, 12, infoFields);
				mpiComm().send(node, 13, chromMap);
			}
		}
#endif
		// get a GenoStructure with parameters. GenoStructure may be shared by some populations
		// a whole set of functions ploidy() etc in GenoStruTriat can be used after this step.
		this->setGenoStructure(ploidy, loci, sexChrom, lociPos, alleleNames,
			lociNames, maxAllele, infoFields, chromMap);

		DBG_DO(DBG_DEVEL, cout << "individual size is " << sizeof(individual) << '+'
			<< sizeof(Allele) << '*' << genoSize() << endl
			<< ", infoPtr: " << sizeof(double*)
			<< ", GenoPtr: " << sizeof(Allele*) << ", Flag: " << sizeof(unsigned char)
			<< ", plus genoStru" << endl );

		try
		{
			// allocate memory here (not in function definition)
			m_inds.resize(m_popSize);

			// create genotype vector holding alleles for all individuals.
#ifdef SIMUMPI
			m_genotype.resize(m_popSize*localGenoSize());
			/// only head node allocate info
			size_t is = localInfoSize();
			DBG_FAILIF(mpiRank() == 0 && is != infoFields.size(), SystemError, "Wrong geno structure");
			m_info.resize(m_popSize*is);
#else
			m_genotype.resize(m_popSize*genoSize());
			/// allocate info
			size_t is = infoSize();
			DBG_ASSERT(is == infoFields.size(), SystemError, "Wrong geno structure");
			m_info.resize(m_popSize*is);
#endif

			// set subpopulation indexes, do not allow popsize change
			setSubPopStru(subPop, false);

			// set individual pointers
			// reset individual pointers
			GenoIterator ptr = m_genotype.begin();
			InfoIterator infoPtr = m_info.begin();
#ifdef SIMUMPI
			UINT step = localGenoSize();
#else
			UINT step = genoSize();
#endif
			for(ULONG i=0; i< m_popSize; ++i, ptr+=step, infoPtr+=is)
			{
				m_inds[i].setGenoPtr(ptr);
				m_inds[i].setGenoStruIdx(genoStruIdx());
				m_inds[i].setShallowCopied(false);
				m_inds[i].setInfoPtr(infoPtr);
			}
		}
		catch(...)
		{
			cout << "Memory allocation fail. A population of size 1 is created." << endl;
			*this = population(0);
			throw OutOfMemory("Memory allocation fail");
		}
		// set local variable
		setRep(-1);
		setGrp(-1);
	}

	population::population(const population& rhs):
	GenoStruTrait(rhs),
		m_popSize(rhs.m_popSize),
		m_numSubPop(rhs.m_numSubPop),
		m_subPopSize(rhs.m_subPopSize),
		m_subPopIndex(rhs.m_subPopIndex),
		m_genotype(0),
		m_info(0),
		m_inds(0),
		m_ancestralDepth(rhs.m_ancestralDepth),
		m_vars(rhs.m_vars),						  // variables will be copied
		m_rep(-1),								  // rep is set to -1 for new pop (until simulator really set them
		m_grp(-1),
		m_gen(0),
		m_curAncestralPop(rhs.m_curAncestralPop),
		m_shallowCopied(false),
		m_infoOrdered(true)
	{
		DBG_DO(DBG_POPULATION,
			cout << "Copy constructor of population is called\n" << endl);

		try
		{
			m_inds.resize(rhs.m_popSize);
#ifdef SIMUMPI
			m_genotype.resize(m_popSize*localGenoSize());
			// have 0 length for mpi/non-head node
			m_info.resize(rhs.m_popSize*localInfoSize());
#else
			m_genotype.resize(m_popSize*genoSize());
			// have 0 length for mpi/non-head node
			m_info.resize(rhs.m_popSize*infoSize());
#endif
		}
		catch(...)
		{
			cout << "Memory allocation fail. A population of size 1 is created." << endl;
			*this = population(0);
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
#ifdef SIMUMPI
		UINT step = localGenoSize();
		UINT infoStep = localInfoSize();
#else
		UINT step = this->genoSize();
		UINT infoStep = this->infoSize();
#endif
		for(ULONG i=0; i< m_popSize; ++i, ptr+=step, infoPtr+=infoStep)
		{
			m_inds[i].setGenoPtr(ptr);
			m_inds[i].setInfoPtr(infoPtr);
			m_inds[i].copyFrom(rhs.m_inds[i]);
		}

		// copy ancestral populations
		try
		{
			// copy all. individual will be shallow copied
			m_ancestralPops = rhs.m_ancestralPops;
			// need to setGenoPtr
			for( size_t ap = 0; ap < m_ancestralPops.size(); ++ap)
			{
				popData& lp = m_ancestralPops[ap];
				const popData& rp = rhs.m_ancestralPops[ap];

				vector<individual>& linds = lp.m_inds;
				const vector<individual>& rinds = rp.m_inds;

				GenoIterator lg = lp.m_genotype.begin();
				constGenoIterator rg = rp.m_genotype.begin();

				InfoIterator li = lp.m_info.begin();
				InfoConstIterator ri = rp.m_info.begin();

				ULONG ps = rinds.size();

				for(ULONG i=0; i<ps; ++i)
				{
					linds[i].setGenoPtr( rinds[i].genoPtr() - rg + lg );
					linds[i].setInfoPtr( rinds[i].infoPtr() - ri + li );
				}
			}
		}
		catch(...)
		{
			cout << "Unable to copy ancestral populations. "
				<< "The popolation size may be too big." << endl
				<< "The population will still be usable but without any ancestral population stored." << endl;
			m_ancestralDepth = 0;
			m_ancestralPops.clear();
		}

		// set local variable
		setRep(-1);
		setGrp(-1);
	}

	///
	population * population::clone(int keepAncestralPops) const
	{
		population * p = new population(*this);
		int oldDepth = m_ancestralDepth;
		if(keepAncestralPops >= 0)
			// try to remove excessive ancestra generations.
			p->setAncestralDepth(keepAncestralPops);
		p->setAncestralDepth(oldDepth);
		return p;
	}

#ifdef SIMUMPI
	/// CPPONLY
	/// info iterator
	/// if order=true, keep order,
	/// if flase, do not respect pop structure
	GappedInfoIterator population::infoBegin(UINT idx, bool order)
	{
		CHECKRANGELOCALINFO(idx);
		if(order && !infoOrdered())
			adjustInfoPosition(true);
		return GappedInfoIterator(m_info.begin()+idx, localInfoSize());
	}

	/// CPPONLY
	GappedInfoIterator population::infoEnd(UINT idx, bool order)
	{
		CHECKRANGELOCALINFO(idx);
		if(order && !infoOrdered())
			adjustInfoPosition(true);
		return GappedInfoIterator(m_info.begin()+idx+m_info.size(), localInfoSize());
	}

	/// info iterator
	/// oder = true: keep order
	/// otherwise, respect subpop structure
	GappedInfoIterator population::infoBegin(UINT index, UINT subPop, bool order)
	{
		CHECKRANGELOCALINFO(index);
		CHECKRANGESUBPOP(subPop);

		if(!infoOrdered())
			adjustInfoPosition(order);

		return GappedInfoIterator(m_info.begin()+index+m_subPopIndex[subPop]*localInfoSize(), localInfoSize());
	}

	///
	GappedInfoIterator population::infoEnd(UINT index, UINT subPop, bool order)
	{
		CHECKRANGELOCALINFO(index);
		CHECKRANGESUBPOP(subPop);

		if(!infoOrdered())
			adjustInfoPosition(order);

		return GappedInfoIterator(m_info.begin()+index+m_subPopIndex[subPop+1]*localInfoSize(), localInfoSize());
	}

	vectorinfo population::indInfo(UINT idx, bool order)
	{
		vectorinfo info;
		if(mpiRank() == 0)
			vectorinfo info = vectorinfo(infoBegin(idx, order), infoEnd(idx, order));
		broadcast(mpiComm(), info, 0);
		return info;
	}

	vectorinfo population::indInfo(const string& name, bool order)
	{
		vectorinfo info;
		if(mpiRank() == 0)
		{
			UINT idx = localInfoIdx(name);
			info = vectorinfo(infoBegin(idx, order), infoEnd(idx, order));
		}
		broadcast(mpiComm(), info, 0);
		return info;
	}

	vectorinfo population::indInfo(UINT idx, UINT subPop, bool order)
	{
		vectorinfo info;
		if(mpiRank() == 0)
			info = vectorinfo(infoBegin(idx, subPop, order), infoEnd(idx, subPop, order));
		broadcast(mpiComm(), info, 0);
		return info;
	}

	vectorinfo population::indInfo(const string& name, UINT subPop, bool order)
	{
		vectorinfo info;
		if(mpiRank() == 0)
		{
			UINT idx = localInfoIdx(name);
			info = vectorinfo(infoBegin(idx, subPop, order), infoEnd(idx, subPop, order));
		}
		broadcast(mpiComm(), info, 0);
		return info;
	}

	/// if order: keep order
	/// otherwise: do not respect subpop info
	PyObject* population::arrIndInfo(bool order)
	{
		if(mpiRank() == 0)
		{
			if(order && !infoOrdered())
				adjustInfoPosition(true);

			return Info_Vec_As_NumArray(m_info.begin(), m_info.end());
		}
		else
			DBG_ASSERT(false, ValueError,
				"Only head node has info information");
	}

	/// if order: keep order
	/// otherwise: respect subpop info
	PyObject* population::arrIndInfo(UINT subPop, bool order)
	{
		CHECKRANGESUBPOP(subPop);

		if(mpiRank() == 0)
		{
			if(!infoOrdered())
				adjustInfoPosition(order);

			return Info_Vec_As_NumArray(m_info.begin() + m_subPopIndex[subPop]*localInfoSize(),
				m_info.begin() + m_subPopIndex[subPop+1]*localInfoSize());
		}
		else
			DBG_ASSERT(false, ValueError,
				"Only head node has info information");
	}
#endif

	int population::__cmp__(const population& rhs) const
	{
#ifdef SIMUMPI
		bool res = 0;
		if( genoStruIdx() != rhs.genoStruIdx() )
		{
			DBG_DO(DBG_POPULATION, cout << "Genotype structures are different" << endl);
			res = 1;
		}

		if(res == 0 && popSize() != rhs.popSize() )
		{
			DBG_DO(DBG_POPULATION, cout << "Population sizes are different" << endl);
			res = 1;
		}

		if (res == 0)
		{
			for( ULONG i=0, iEnd = popSize(); i < iEnd; ++i)
				if( m_inds[i] != rhs.m_inds[i])
			{
				DBG_DO(DBG_POPULATION, cout << "Individuals are different" << endl);
				res = 1;
				break;
			}
		}
		// now collect result for all nodes
		bool allRes = 0;
		if (mpiRank() == 0)
			reduce(mpiComm(), res, allRes, std::logical_or<bool>(), 0);
		else
			reduce(mpiComm(), res, std::logical_or<bool>(), 0);
		broadcast(mpiComm(), !!allRes, 0);

		return allRes;
#else
		if( genoStruIdx() != rhs.genoStruIdx() )
		{
			DBG_DO(DBG_POPULATION, cout << "Genotype structures are different" << endl);
			return 1;
		}

		if(popSize() != rhs.popSize() )
		{
			DBG_DO(DBG_POPULATION, cout << "Population sizes are different" << endl);
			return 1;
		}

		for( ULONG i=0, iEnd = popSize(); i < iEnd; ++i)
			if( m_inds[i] != rhs.m_inds[i])
		{
			DBG_DO(DBG_POPULATION, cout << "Individuals are different" << endl);
			return 1;
		}

		return 0;
#endif
	}

	PyObject* population::arrGenotype(bool order)
	{
		if(shallowCopied() && order)
			// adjust position. deep=true
			adjustGenoPosition(true);
#ifdef SIMUMPI
		// shift (starting point), size (total size)
		// trunk size, pieces map
		return Allele_Vec_As_NumArray(0, genoSize()*popSize(),
			totNumLoci(), locusMap());
#else
		// directly expose values. Do not copy data over.
		return Allele_Vec_As_NumArray(m_genotype.begin(), m_genotype.end());
#endif
	}

	/// get the whole genotype.
	/// individuals will be in order before exposing
	/// their genotypes.
	///
	/// if order: keep order
	/// otherwise: respect subpop structure
	PyObject* population::arrGenotype(UINT subPop, bool order)
	{
		CHECKRANGESUBPOP(subPop);
		if(shallowCopied())
			adjustGenoPosition(order);
#ifdef SIMUMPI
		return Allele_Vec_As_NumArray(m_subPopIndex[subPop]*genoSize(),
			genoSize()*subPopSize(subPop), totNumLoci(), locusMap());
#else
		return Allele_Vec_As_NumArray( genoBegin(subPop, order), genoEnd(subPop, order));
#endif
	}

	void population::setSubPopStru(const vectorlu& newSubPopSizes, bool allowPopSizeChange)
	{
		// case 1: remove all subpopulation structure
		// do not change population size
		// individuals are valid....
		if ( newSubPopSizes.empty() )
		{
			m_numSubPop = 1;
			m_subPopSize.resize(1, m_popSize);
			m_subPopIndex.resize(2);
		}
		else									  // may change populaiton size
		{
			m_numSubPop = newSubPopSizes.size();
			m_subPopSize = newSubPopSizes;
			m_subPopIndex.resize( m_numSubPop+1);

			ULONG totSize = accumulate(newSubPopSizes.begin(), newSubPopSizes.end(), 0UL);

			// usually, totSize == m_popSize, individuals are valid
			if( totSize != m_popSize)
			{
				DBG_DO( DBG_POPULATION, "Populaiton size changed to " + toStr(totSize) +
					" Genotype information may be lost");

#ifndef OPTIMIZED
				if( !allowPopSizeChange)
				{
					DBG_DO(DBG_POPULATION, cout << "Total new size " << totSize << endl);
					DBG_DO(DBG_POPULATION, cout << "Attempted subpop " << newSubPopSizes << endl);
					DBG_DO(DBG_POPULATION, cout << "Total current " << m_popSize << endl);
					DBG_DO(DBG_POPULATION, cout << "Current subpop size " <<
						this->subPopSizes() << endl);
					throw ValueError("Populaiton size is fixed (by allowPopSizeChange parameter).\n"
						" Subpop sizes should add up to popsize");
				}
#endif

				// change populaiton size
				// genotype and individual info will be kept
				// but pointers need to be recalibrated.
				m_popSize = totSize;

				try
				{
#ifdef SIMUMPI
					m_genotype.resize(m_popSize*localGenoSize());
					m_info.resize(m_popSize*localInfoSize());
#else
					m_genotype.resize(m_popSize*genoSize());
					m_info.resize(m_popSize*infoSize());
#endif
					m_inds.resize(m_popSize);

				}
				catch(...)
				{
					throw OutOfMemory("Memory allocation fail");
				}
				// reset individual pointers
				GenoIterator ptr = m_genotype.begin();
				InfoIterator infoPtr = m_info.begin();
#ifdef SIMUMPI
				UINT step = localGenoSize();
				UINT is = localInfoSize();
#else
				UINT step = genoSize();
				UINT is = infoSize();
#endif
				for(ULONG i=0; i< m_popSize; ++i, ptr+=step, infoPtr+=is)
				{
					m_inds[i].setGenoPtr( ptr );
					m_inds[i].setInfoPtr( infoPtr );
					m_inds[i].setGenoStruIdx(genoStruIdx());
					m_inds[i].setShallowCopied(false);
				}
				m_shallowCopied = false;
				m_infoOrdered = true;
			}
		}

		// build subPop index
		UINT i = 1;
		for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
			m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];
	}

	void population::setSubPopByIndID(vectori id)
	{
		if( !id.empty())
		{
			DBG_ASSERT( id.size() == m_popSize, ValueError,
				"Info should have the same length as pop size");
			for(ULONG it=0; it < m_popSize; ++it)
				ind(it).setSubPopID( id[it] );
		}

		DBG_DO(DBG_POPULATION, cout << "Sorting individuals."<< endl);
		// sort individuals first
		std::sort(indBegin(), indEnd());
		setShallowCopied(true);
		setInfoOrdered(false);

		// sort individuals first
		// remove individuals with negative index.
		if( indBegin()->subPopID() < 0 )
		{
			// popsize etc will be changed.
			ULONG newPopSize = m_popSize;
			IndIterator it=indBegin();
			for(; it != indEnd();  ++it)
			{
				if( it->subPopID() < 0 )
					newPopSize -- ;
				else
					break;
			}
			// 'it' now point to the one with positive subPopID()

			DBG_DO(DBG_POPULATION, cout << "New pop size" << newPopSize << endl);

			// allocate new genotype and inds
#ifdef SIMUMPI
			vectora newGenotype(localGenoSize() * newPopSize);
			vectorinfo newInfo(newPopSize*localInfoSize());
#else
			vectora newGenotype(genoSize() * newPopSize);
			vectorinfo newInfo(newPopSize*infoSize());
#endif
			vector<individual> newInds(newPopSize);

			DBG_ASSERT( indEnd()== newPopSize+it, SystemError,
				"Pointer misplaced. ");

			// assign genotype location and set structure information for individuals
			GenoIterator ptr = newGenotype.begin();
			InfoIterator infoPtr = newInfo.begin();
#ifdef SIMUMPI
			UINT step = localGenoSize();
			UINT infoStep = localInfoSize();
#else
			UINT step = genoSize();
			UINT infoStep = infoSize();
#endif
			for(ULONG i=0; i< newPopSize; ++i, ptr+=step, ++it, infoPtr+=infoStep)
			{
				newInds[i].setGenoStruIdx(genoStruIdx());
				newInds[i].setGenoPtr( ptr );
				newInds[i].setInfoPtr( infoPtr );
				newInds[i].copyFrom(*it);		  // copy everything, with info value
			}

			// now, switch!
			m_genotype.swap(newGenotype);
			m_info.swap(newInfo);
			m_inds.swap(newInds);

			m_popSize = newPopSize;
			setShallowCopied(false);
			setInfoOrdered(true);
		}

		if( m_inds.empty())
		{
			m_numSubPop = 1;
			m_subPopSize.resize(1,0);
			m_subPopIndex.resize(2);
		}
		else
		{
			// reset indices etc.
			m_numSubPop = static_cast<UINT>(m_inds.back().subPopID())+1;
			m_subPopSize.resize(m_numSubPop);
			m_subPopIndex.resize(m_numSubPop+1);

			// check subpop size
			fill(m_subPopSize.begin(), m_subPopSize.end(), 0);
			for(IndIterator it = indBegin(), itEnd = indEnd(); it < itEnd;  ++it)
				m_subPopSize[ static_cast<UINT>(it->subPopID()) ] ++;
		}
		/// rebuild index
		size_t i = 1;
		for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
			m_subPopIndex[i] = m_subPopIndex[i-1] + m_subPopSize[i - 1];
	}

	void population::splitSubPop(UINT which, vectorlu sizes, vectoru subPopID)
	{
		DBG_ASSERT( accumulate(sizes.begin(), sizes.end(), 0UL) == subPopSize(which),
			ValueError,
			"Sum of subpopulation sizes does not equal to the size of subpopulation to be splitted.");

		DBG_FAILIF( !subPopID.empty() && subPopID.size() != sizes.size(), ValueError,
			"If subPopID is given, it should have the same length as subPOP");

		if( sizes.size() == 1)
			return;

		// set initial info
		setIndSubPopIDWithID();

		UINT spID;
		if(subPopID.empty())					  // starting sp number
			spID = which;
		else
		{
			spID = subPopID[0];
			DBG_WARNING( spID != which && spID < numSubPop(),
				"new subpop ID is already used. You are effectively merging two subpopulations")
		}
		ULONG sz=0;								  // idx within subpop
		size_t newSPIdx=0;
		for(IndIterator ind = indBegin(which); ind != indEnd(which); ++ind)
		{
			if( sz == sizes[newSPIdx])
			{
				sz = 0;
				newSPIdx++;
				if(subPopID.empty())
					spID = numSubPop()+newSPIdx-1;
				else
				{
					DBG_WARNING( subPopID[newSPIdx] != which && subPopID[newSPIdx] < numSubPop(),
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
		DBG_ASSERT( fcmp_eq(accumulate(proportions.begin(), proportions.end(), 0.), 1.), ValueError,
			"Proportions do not add up to one.");

		if( proportions.size() == 1)
			return;

		ULONG spSize = subPopSize(which);
		vectorlu subPop(proportions.size());
		for(size_t i=0; i< proportions.size()-1; ++i)
			subPop[i] = static_cast<ULONG>(floor(spSize*proportions[i]));
		// to avoid round off problem, calculate the last subpopulation
		subPop[ subPop.size()-1] = spSize - accumulate( subPop.begin(), subPop.end()-1, 0L);
		splitSubPop(which, subPop, subPopID);
	}

	void population::removeEmptySubPops()
	{
		// if remove empty subpops
		UINT newSPNum = m_numSubPop;
		vectorlu newSPSize;
		for(size_t sp=0; sp < m_numSubPop; ++sp)
		{
			if( m_subPopSize[sp] == 0 )
				newSPNum --;
			else
				newSPSize.push_back( m_subPopSize[sp]);
		}
		m_numSubPop = newSPNum;
		m_subPopSize.swap(newSPSize);
		m_subPopIndex.resize(m_numSubPop+1);
		/// rebuild index
		size_t i = 1;
		for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
			m_subPopIndex[i] = m_subPopIndex[i-1] + m_subPopSize[i - 1];
	}

	void population::removeSubPops(const vectoru& subPops, bool shiftSubPopID, bool removeEmptySubPops)
	{
#ifndef OPTIMIZED
		// check if subPops are valid
		for( vectoru::const_iterator sp = subPops.begin(); sp < subPops.end(); ++sp)
		{
			DBG_WARNING(*sp >= m_numSubPop, "Subpopulation " + toStr(*sp) + " does not exist.");
		}
#endif
		setIndSubPopIDWithID();
		int shift=0;
		for( size_t sp = 0; sp < m_numSubPop; ++sp)
		{
			if( find( subPops.begin(), subPops.end(), sp) != subPops.end())
			{
				shift++;
				for(IndIterator ind = indBegin(sp); ind != indEnd(sp); ++ind)
					ind->setSubPopID(-1);		  // remove
			}
			// other subpop shift left
			else if(shiftSubPopID)
			{
				for(IndIterator ind = indBegin(sp); ind != indEnd(sp); ++ind)
					ind->setSubPopID(sp-shift);	  // shift left
			}
		}

		UINT pendingEmptySubPops = 0;
		for(UINT i=m_numSubPop-1; i>=0 && (subPopSize(i) == 0
			|| find( subPops.begin(), subPops.end(), i) != subPops.end()); --i, ++pendingEmptySubPops);
		setSubPopByIndID();
		// what to do with pending empty subpops?
		if( pendingEmptySubPops != 0 && ! removeEmptySubPops )
		{
			vectorlu spSizes = subPopSizes();
			for(UINT i=0; i<pendingEmptySubPops; ++i)
				spSizes.push_back(0);
			setSubPopStru(spSizes, false);
		}
		if(removeEmptySubPops)
			this->removeEmptySubPops();
	}

	void population::removeIndividuals(const vectoru& inds, int subPop, bool removeEmptySubPops)
	{
		setIndSubPopIDWithID();
		if( subPop == -1 )
		{
			for(size_t i = 0; i < inds.size(); ++i)
				ind(inds[i]).setSubPopID(-1);	  // remove
		}
		else
		{
			for(size_t i = 0; i < inds.size(); ++i)
												  // remove
				ind(inds[i], subPop).setSubPopID(-1);
		}

		int oldNumSP = numSubPop();
		setSubPopByIndID();
		int pendingEmptySubPops = oldNumSP - numSubPop();
		// what to do with pending empty subpops?
		if( pendingEmptySubPops != 0 && ! removeEmptySubPops )
		{
			vectorlu spSizes = subPopSizes();
			for(int i=0; i<pendingEmptySubPops; ++i)
				spSizes.push_back(0);
			setSubPopStru(spSizes, false);
		}
		if(removeEmptySubPops)
			this->removeEmptySubPops();
	}

	void population::mergeSubPops(vectoru subPops, bool removeEmptySubPops)
	{
		// set initial info
		setIndSubPopIDWithID();

		// merge all subpopulations
		if(subPops.empty())
		{
			// [ popSize() ]
			vectorlu sz(1, popSize());
			setSubPopStru(sz, false);
			return;
		}

		UINT id = subPops[0];
		for(UINT sp=0; sp < numSubPop(); ++sp)
		{
			if( find( subPops.begin(), subPops.end(), sp) != subPops.end())
				for(IndIterator ind = indBegin(sp); ind != indEnd(sp); ++ind)
					ind->setSubPopID(id);
		}
		int oldNumSP = numSubPop();
		setSubPopByIndID();
		int pendingEmptySubPops = oldNumSP - numSubPop();
		// what to do with pending empty subpops?
		if( pendingEmptySubPops != 0 && ! removeEmptySubPops )
		{
			vectorlu spSizes = subPopSizes();
			for(int i=0; i<pendingEmptySubPops; ++i)
				spSizes.push_back(0);
			setSubPopStru(spSizes, false);
		}

		if( removeEmptySubPops)
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
		for(ULONG i=0; i< newPopSize; ++i, ptr+=step, infoPtr+=infoStep)
		{
			newInds[i].setGenoStruIdx(genoStruIdx());
			newInds[i].setGenoPtr( ptr );
			newInds[i].setInfoPtr( infoPtr );
		}
		// copy stuff over
		for(ULONG i = 0; i < popSize(); ++i)
			newInds[i].copyFrom(ind(i));
		ULONG start = popSize();
		for(ULONG i = 0; i < pop.popSize(); ++i)
			newInds[i+start].copyFrom(pop.ind(i));
		// now, switch!
		m_genotype.swap(newGenotype);
		m_info.swap(newInfo);
		m_inds.swap(newInds);
		m_popSize = newPopSize;
		setShallowCopied(false);
		setInfoOrdered(true);
		if (newSubPopSizes.empty())
			m_subPopSize = newSS;
		else
			m_subPopSize = newSubPopSizes;
		m_numSubPop = m_subPopSize.size();
		/// rebuild index
		m_subPopIndex.resize(m_numSubPop+1);
		size_t i = 1;
		for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
			m_subPopIndex[i] = m_subPopIndex[i-1] + m_subPopSize[i - 1];
	}

	void population::mergePopulation(const population & pop, const vectorlu & newSubPopSizes,
		int keepAncestralPops)
	{
		DBG_FAILIF(ancestralDepth() != pop.ancestralDepth(), ValueError,
			"Merged populations should have the same number of ancestral generations");
		UINT topGen;
		if(keepAncestralPops < 0 || static_cast<UINT>(keepAncestralPops) >= ancestralDepth())
			topGen = ancestralDepth();
		else
			topGen = keepAncestralPops;
		// go to the oldest generation
		useAncestralPop(topGen);
		const_cast<population &>(pop).useAncestralPop(topGen);
		mergePopulationPerGen(pop, newSubPopSizes);
		if(topGen > 0)
		{
			for(int depth = topGen - 1; depth >=0; --depth)
			{
				useAncestralPop(depth);
				const_cast<population &>(pop).useAncestralPop(depth);
				mergePopulationPerGen(pop, newSubPopSizes);
			}
		}
		useAncestralPop(0);
		const_cast<population &>(pop).useAncestralPop(0);
	}

	void population::mergePopulationByLoci(const population & pop,
		const vectoru & newNumLoci, const vectorf & newLociPos)
	{
		DBG_FAILIF(subPopSizes() != pop.subPopSizes(), ValueError,
			"Merged population should have the same number of individuals in each subpopulation");

		UINT gs1 = totNumLoci();
		UINT gs2 = pop.totNumLoci();

		// obtain new genotype structure and set it
		setGenoStructure(mergeGenoStru(pop.genoStruIdx()));

		DBG_FAILIF(ancestralDepth() != pop.ancestralDepth(), ValueError,
			"Merged populations should have the same number of ancestral generations");
		for(int depth = ancestralDepth(); depth >=0; --depth)
		{
			useAncestralPop(depth);
			const_cast<population &>(pop).useAncestralPop(depth);
			//
			ULONG newPopGenoSize = genoSize() * m_popSize;
			vectora newGenotype(newPopGenoSize);

			// copy data over
			GenoIterator ptr = newGenotype.begin();
			UINT pEnd = ploidy();
			for(ULONG i=0; i< m_popSize; ++i)
			{
				// set new geno structure
				m_inds[i].setGenoStruIdx(genoStruIdx());
				GenoIterator ptr1 = m_inds[i].genoPtr();
				GenoIterator ptr2 = pop.m_inds[i].genoPtr();
				// new genotype
				m_inds[i].setGenoPtr( ptr );
				// copy each allele
				for(UINT p=0; p < pEnd; ++p)
				{
					for(size_t j=0; j < gs1; ++j)
						*(ptr++) = *(ptr1++);
					for(size_t j=0; j < gs2; ++j)
						*(ptr++) = *(ptr2++);
				}
			}
			m_genotype.swap(newGenotype);
		}

		// reset genoStructure again
		if (!newNumLoci.empty() and ! newLociPos.empty())
			rearrangeLoci(newNumLoci, newLociPos);

		setShallowCopied(false);
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
		for(ULONG i=0; i< newPopSize; ++i, ptr+=step, infoPtr+=infoStep)
		{
			newInds[i].setGenoStruIdx(genoStruIdx());
			newInds[i].setGenoPtr( ptr );
			newInds[i].setInfoPtr( infoPtr );
		}
		// copy stuff over
		ULONG startSP = 0;
		for(UINT sp=0; sp < numSubPop(); ++sp)
		{
			ULONG spSize = subPopSize(sp);
			for(ULONG i =0, j = 0; i < newSubPopSizes[sp]; ++j, ++i)
			{
				// repeating?
				if ((j / spSize) > 0 && ! propagate)
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
		setShallowCopied(false);
		setInfoOrdered(true);
		m_subPopSize = newSubPopSizes;
		/// rebuild index
		size_t i = 1;
		for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
			m_subPopIndex[i] = m_subPopIndex[i-1] + m_subPopSize[i - 1];
	}

	void population::reorderSubPops(const vectoru& order, const vectoru& rank,
		bool removeEmptySubPops)
	{
		DBG_FAILIF( order.empty() && rank.empty(), ValueError,
			"Please specify one of order or rank.");

		DBG_FAILIF( !order.empty() && !rank.empty(), ValueError,
			"You can specify only one of order or rank.");

		if(removeEmptySubPops)
			this->removeEmptySubPops();

		if( ( !order.empty() && order.size() != m_numSubPop )
			|| ( ! rank.empty() && rank.size() != m_numSubPop))
			cout << "Warning: Given order or rank does not have the length of number of subpop." << endl;

		if( !order.empty())
		{
			// alow order[i] > numSubPop(). In a special case, I have last empty subpop...
			for(size_t i=0; i < order.size(); ++i)
			{
				if(order[i] >= numSubPop())
					continue;
				for(IndIterator ind = indBegin(order[i]); ind != indEnd(order[i]); ++ind)
					ind->setSubPopID(i);
			}
		}
		else
		{
			for(size_t i=0; i < rank.size(); ++i)
			{
				if(i >= numSubPop())
					continue;
				for(IndIterator ind = indBegin(i); ind != indEnd(i); ++ind)
					ind->setSubPopID(rank[i]);
			}
		}
		// reset ...
		setSubPopByIndID();
	}

	population& population::newPopByIndIDPerGen(const vectori& id, bool removeEmptySubPops)
	{
		// determine the size of needed individuals
		vectorlu sz;
		if(!id.empty())
		{
			DBG_ASSERT(id.size() == popSize(), ValueError, "Please assign id for each individual");
			for(ULONG i = 0; i != id.size(); ++i)
			{
				if(id[i] < 0)
					continue;
				if(static_cast<UINT>(id[i]) >= sz.size())
					sz.resize(id[i]+1);
				sz[id[i]]++;
			}
		}
		else
		{
			for(UINT sp = 0; sp < numSubPop(); ++sp)
			{
				for(IndIterator it = indBegin(sp); it != indEnd(sp); ++it)
				{
					int indID = it->subPopID();
					if(indID < 0)
						continue;
					if(static_cast<UINT>(indID) >= sz.size())
						sz.resize(indID+1);
					sz[indID]++;
				}
			}
		}
		DBG_DO(DBG_POPULATION, cout << "newPopByIndIDPerGen: New population size: " << sz << endl);

		// create a population with this size
		population * pop = new population(0, ploidy(), numLoci(), sexChrom(), lociPos(), sz, 0,
			alleleNames(), lociNames(), maxAllele(), infoFields(), chromMap());
		// copy individuals over
		IndIterator from = indBegin();
		vector<IndIterator> to;
		for(UINT sp=0; sp < sz.size(); ++sp)
			to.push_back(pop->indBegin(sp));
		if(!id.empty())
		{
			for(ULONG i = 0; i != id.size(); ++i)
			{
				if(id[i] >= 0)
				{
					to[id[i]]->copyFrom(ind(i));
					++to[id[i]];
				}
			}
		}
		else
		{
			for(; from != indEnd(); ++from)
			{
				int indID = from->subPopID();
				if(indID >= 0)
				{
					to[indID]->copyFrom(*from);
					++to[indID];
				}
			}
		}
		if(removeEmptySubPops)
			pop->removeEmptySubPops();
		return *pop;
	}

	/** form a new population according to info, info can be given directly */
	population& population::newPopByIndID(int keepAncestralPops,
		const vectori& id, bool removeEmptySubPops)
	{
		UINT topGen;
		if(keepAncestralPops < 0 || static_cast<UINT>(keepAncestralPops) >= ancestralDepth())
			topGen = ancestralDepth();
		else
			topGen = keepAncestralPops;
		// go to the oldest generation
		useAncestralPop(topGen);
		population & ret = newPopByIndIDPerGen(id, removeEmptySubPops);
		// prepare for push and discard
		ret.setAncestralDepth(topGen);
		if(topGen > 0)
		{
			for(int depth = topGen - 1; depth >=0; --depth)
			{
				useAncestralPop(depth);
				ret.pushAndDiscard(newPopByIndIDPerGen(id, removeEmptySubPops));
			}
		}
		return ret;
	}

	void population::removeLoci( const vectoru& remove, const vectoru& keep)
	{
		DBG_FAILIF( !keep.empty() && ! remove.empty(), ValueError,
			"Please specify one and only one of keep or remove.");

		if( keep.empty() && remove.empty() )
			return;

		vectoru loci;
		if( !keep.empty())
			loci = keep;
		else
		{
			for(size_t loc = 0; loc < this->totNumLoci(); ++loc)
				// if not removed
				if( find(remove.begin(), remove.end(), loc) == remove.end())
					loci.push_back(loc);
		}

#ifndef OPTIMIZED
		for(size_t i=0; i<loci.size(); ++i)
		{
			DBG_FAILIF( loci[i] >= this->totNumLoci(), ValueError,
				"Given loci " + toStr(loci[i]) + " exceed max number of loci." );
			DBG_FAILIF( i > 0 && loci[i] <= loci[i-1], ValueError,
				"Given loci should be in order.");
		}
#endif
		// adjust order before doing anything
		UINT newTotNumLoci = loci.size();
#ifdef SIMUMPI
		UINT oldTotNumLoci = localNumLoci();
#else
		UINT oldTotNumLoci = totNumLoci();
#endif

		// first, new genotype
		// get a GenoStructure with parameters. GenoStructure may be shared by some populations
		// a whole set of functions ploidy() etc in GenoStruTriat can be used after this step.
		vectoru newNumLoci;
		vectorf newLociDist;
		vectorstr newLociNames;
		UINT curCh = 9999;						  // not 0, will be set to 0 soon.
		for(vectoru::iterator loc = loci.begin();
			loc != loci.end(); ++loc)
		{
			UINT ch = this->chromLocusPair(*loc).first;
			if( newNumLoci.empty() || curCh != ch )
			{
				newNumLoci.push_back(1);
				curCh = ch;
			}
			else
				newNumLoci.back()++;
			newLociDist.push_back( this->locusPos(*loc));
			newLociNames.push_back( this->locusName(*loc));
		}

		// prepare data
		//
		// keep m_popSize;
		// keep m_numSubPop;
		// keep m_subPopSize;
		// keep m_subPopIndex;

		// genotype
		// allocate new genotype and inds
		// new geno structure is in effective now!
#ifdef SIMUMPI
		// for the MPI case, genotype may be moved from one node to another
		//
		// For this reason, geno structure is not set as in the non-MPI case.
		// A temp population is created to hold new geno structure
		population tmpPop = population(0, ploidy(), newNumLoci, sexChrom(), newLociDist,
			vectorlu(0), 0, vectorstr(), newLociNames, maxAllele(), infoFields(), chromMap());

		ULONG newPopGenoSize = tmpPop.localGenoSize() * m_popSize;
		vectora newGenotype(newPopGenoSize);
		// copy data over
		GenoIterator ptr = newGenotype.begin();
		UINT pEnd = ploidy();
		// rank map
		UINT rnk = mpiRank();
		UINT oldBeginLocus = beginLocus();
		vectoru oldRank(totNumLoci());
		vectoru newRank(totNumLoci());
		for(size_t i=0; i < totNumLoci(); ++i)
			oldRank[i] = rankOfLocus(i);
		for(size_t i=0; i < loci.size(); ++i)
			newRank[loci[i]] = tmpPop.rankOfLocus(i);
		//
		for(ULONG i=0; i < m_popSize; ++i)
		{
			GenoIterator oldPtr = m_inds[i].genoPtr();
			// new genotype
			m_inds[i].setGenoPtr( ptr );
			// copy each chromosome
			for(UINT p=0; p < pEnd; ++p)
			{
				for(vectoru::iterator loc = loci.begin();
					loc != loci.end(); ++loc)
				{
					// belong to the same rank
					if (oldRank[*loc] == newRank[*loc])
					{
						if(newRank[*loc] == rnk)
							*(ptr++) = oldPtr[*loc - oldBeginLocus];
					}
					else						  // cross rank copy
					{
						// bring allele from another node
						if(rnk == newRank[*loc])
						{
							Allele ale;
							// second parameter is tag
							mpiComm().recv(oldRank[*loc], *loc, ale);
							*(ptr++) = ale;
						}
						else if(rnk == oldRank[*loc])
							mpiComm().send(newRank[*loc], *loc,
									oldPtr[*loc - beginLocus()]);
					}
				}
				oldPtr += oldTotNumLoci;		  // next ploidy
			}
		}
		m_genotype.swap(newGenotype);

		// ancestral populations?
		for(size_t ap=0; ap < m_ancestralPops.size(); ++ap)
		{
			popData& p = m_ancestralPops[ap];
			// set pointers
			vector<individual>& inds = p.m_inds;
			ULONG ps = inds.size();
			vectora newGenotype(ps*pEnd*newTotNumLoci);
			ptr = newGenotype.begin();

			for(ULONG i=0; i< ps; ++i)
			{
				// set new geno structure
				inds[i].setGenoStruIdx(genoStruIdx());
				GenoIterator oldPtr = inds[i].genoPtr();
				// new genotype
				inds[i].setGenoPtr( ptr );
				// copy each chromosome
				for(UINT p=0; p < pEnd; ++p)
				{
					for(vectoru::iterator loc = loci.begin();
						loc != loci.end(); ++loc)
					{
						// belong to the same rank
						if (oldRank[*loc] = newRank[*loc])
						{
							if(newRank[*loc] == rnk)
							{
								*(ptr++) = oldPtr[*loc - oldBeginLocus];
							}
						}
						else					  // cross rank copy
						{
							// bring allele from another node
							if(rnk == newRank[*loc])
							{
								Allele ale;
								// second parameter is tag
								mpiComm().recv(oldRank[*loc], *loc, ale);
								*(ptr++) = ale;
							}
							else if(rnk == oldRank[*loc])
							{
								mpiComm().send(newRank[*loc], *loc,
									oldPtr[*loc - beginLocus()]);
							}
						}
					}
					oldPtr += oldTotNumLoci;	  // next ploidy
				}
			}
			p.m_genotype.swap(newGenotype);
		}
		this->setGenoStructure(this->ploidy(), newNumLoci, this->sexChrom(), newLociDist,
			this->alleleNames(), newLociNames, this->maxAllele(), this->infoFields(), this->chromMap() );
		// now set geno structure
		for(ULONG i=0; i< m_popSize; ++i)
		{
			// set new geno structure
			m_inds[i].setGenoStruIdx(genoStruIdx());
		}
#else
		this->setGenoStructure(this->ploidy(), newNumLoci, this->sexChrom(), newLociDist,
			this->alleleNames(), newLociNames, this->maxAllele(), this->infoFields(), this->chromMap() );
		ULONG newPopGenoSize = genoSize() * m_popSize;
		vectora newGenotype(newPopGenoSize);
		// keep newInds();

		// copy data over
		GenoIterator ptr = newGenotype.begin();
		UINT pEnd = ploidy();
		for(ULONG i=0; i< m_popSize; ++i)
		{
			// set new geno structure
			m_inds[i].setGenoStruIdx(genoStruIdx());
			GenoIterator oldPtr = m_inds[i].genoPtr();
			// new genotype
			m_inds[i].setGenoPtr( ptr );
			// copy each chromosome
			for(UINT p=0; p < pEnd; ++p)
			{
				for(vectoru::iterator loc = loci.begin();
					loc != loci.end(); ++loc)
				{
					*(ptr++) = oldPtr[*loc];
				}
				oldPtr += oldTotNumLoci;		  // next ploidy
			}
		}
		m_genotype.swap(newGenotype);

		// ancestral populations?
		for(size_t ap=0; ap < m_ancestralPops.size(); ++ap)
		{
			popData& p = m_ancestralPops[ap];
			// set pointers
			vector<individual>& inds = p.m_inds;
			ULONG ps = inds.size();
			vectora newGenotype(ps*pEnd*newTotNumLoci);
			ptr = newGenotype.begin();

			for(ULONG i=0; i< ps; ++i)
			{
				// set new geno structure
				inds[i].setGenoStruIdx(genoStruIdx());
				GenoIterator oldPtr = inds[i].genoPtr();
				// new genotype
				inds[i].setGenoPtr( ptr );
				// copy each chromosome
				for(UINT p=0; p < pEnd; ++p)
				{
					for(vectoru::iterator loc = loci.begin();
						loc != loci.end(); ++loc)
					{
						*(ptr++) = oldPtr[*loc];
					}
					oldPtr += oldTotNumLoci;	  // next ploidy
				}
			}
			p.m_genotype.swap(newGenotype);
		}										  // all ancestral
#endif
		setShallowCopied(false);
		setInfoOrdered(true);
	}

	/** get a new population with selected loci */
	population& population::newPopWithPartialLoci(
		const vectoru& remove,
		const vectoru& keep)
	{
		// copy the population over (info is also copied)
		population* pop = new population(*this);
		pop->removeLoci(remove, keep);
		return *pop;
	}

	void population::rearrangeLoci(const vectoru & newNumLoci, const vectorf & newLociPos)
	{
		/// total number of loci can not change
		DBG_FAILIF(std::accumulate(newNumLoci.begin(), newNumLoci.end(), 0U) != totNumLoci(), ValueError,
			"Re-arrange loci must keep the same total number of loci");
		setGenoStructure(ploidy(), newNumLoci, sexChrom(), newLociPos,
			alleleNames(), lociNames(), maxAllele(), infoFields(),
			chromMap());
		for(int depth = ancestralDepth(); depth >=0; --depth)
		{
			useAncestralPop(depth);

			// now set geno structure
			for(ULONG i=0; i < m_popSize; ++i)
				// set new geno structure
				m_inds[i].setGenoStruIdx(genoStruIdx());
		}
	}

	void population::pushAndDiscard(population& rhs, bool force)
	{
		// time consuming!
		DBG_ASSERT( rhs.genoStruIdx() == genoStruIdx(), ValueError,
			"Passed population has different genotypic structure");

		DBG_ASSERT( m_genotype.begin() != rhs.m_genotype.begin(), ValueError,
			"Passed population is a reference of current population, swapPop failed." );

		// front -1 pop, -2 pop, .... end
		//
		if( !force && m_ancestralDepth > 0
			&& ancestralDepth() == static_cast<size_t>(m_ancestralDepth) )
			m_ancestralPops.pop_back();

		// save current population
		if( force || m_ancestralDepth != 0 )
		{
			// add a empty popData
			m_ancestralPops.push_front(popData());
			// get its reference
			popData& pd = m_ancestralPops.front();
			// swap with real data
			// current population may *not* be in order
			pd.m_subPopSize.swap(m_subPopSize);
			// store starting geno ptr,
			// if m_genotype is re-allocated, reset pointers
			// in m_inds
#ifndef OPTIMIZED
			pd.m_startingGenoPtr = m_genotype.begin();
#endif
			pd.m_info.swap(m_info);
			pd.m_genotype.swap(m_genotype);
			pd.m_inds.swap(m_inds);
		}

		// then swap out data
#ifndef OPTIMIZED
		GenoIterator rhsStartingGenoPtr = rhs.m_genotype.begin();
		GenoIterator lhsStartingGenoPtr = m_genotype.begin();
#endif
		m_popSize = rhs.m_popSize;
		m_numSubPop = rhs.m_numSubPop;
		m_subPopSize.swap(rhs.m_subPopSize);
		m_subPopIndex.swap(rhs.m_subPopIndex);
		m_genotype.swap(rhs.m_genotype);
		m_info.swap(rhs.m_info);
		m_inds.swap(rhs.m_inds);
#ifndef OPTIMIZED
		DBG_FAILIF( rhsStartingGenoPtr != m_genotype.begin(),
			SystemError, "Starting genoptr has been changed.");
		DBG_FAILIF( lhsStartingGenoPtr != rhs.m_genotype.begin(),
			SystemError, "Starting genoptr has been changed.");
#endif
		// current population should be working well
		// (with all datamember copied form rhs
		// rhs may not be working well since m_genotype etc
		// may be from ancestral pops
		if( rhs.m_popSize != rhs.m_inds.size())
		{
			// keep size if pop size is OK.
			// remove all supopulation structure of rhs
			rhs.m_popSize = rhs.m_inds.size();
			rhs.m_numSubPop = 1;
			rhs.m_subPopSize.resize(1, rhs.m_popSize);
			rhs.m_subPopIndex.resize(2,0);
			rhs.m_subPopIndex[1] = rhs.m_popSize;
			// no need to set genoPtr or genoStru()
		}
	}

	/// add field
	void population::addInfoField(const string field, double init)
	{
#ifdef SIMUMPI
		if(mpiRank()==0)
		{
			DBG_ASSERT(m_info.size() == localInfoSize()*popSize(), SystemError,
				"Info size is wrong");
			UINT os = localInfoSize();
			UINT idx;
			// if this field exists, return directly
			try
			{
				// for node 0, local info index is the major index.
				idx = localInfoIdx(field);
				// only needs to initialize
				int oldAncPop = m_curAncestralPop;
				for(UINT anc=0; anc <= m_ancestralPops.size(); anc++)
				{
					useAncestralPop(anc);
					for(IndIterator ind=indBegin(); ind!=indEnd(); ++ind)
						ind->setInfo(init, idx);
				}
				useAncestralPop(oldAncPop);
				return;
			}
			catch(IndexError &)
			{
				// and we only add field to node 0
				struAddInfoField(field);
			}

			/// adjust information size.
			UINT is = localInfoSize();
			if(os != is)
			{
				int oldAncPop = m_curAncestralPop;
				for(UINT anc=0; anc <= m_ancestralPops.size(); anc++)
				{
					useAncestralPop(anc);
					vectorinfo newInfo(is*popSize());
					/// copy the old stuff in
					InfoIterator ptr = newInfo.begin();
					for(IndIterator ind=indBegin(); ind!=indEnd(); ++ind)
					{
						copy(ind->infoBegin(), ind->infoBegin() + is - 1, ptr);
						ind->setInfoPtr(ptr);
						fill(ind->infoBegin() + os, ind->infoEnd(), init);
						ptr += is;
					}
					m_info.swap(newInfo);
				}
				useAncestralPop(oldAncPop);
			}
			return;
		}
		// for other node, just return 0.
		else
			return;
#else

		DBG_ASSERT(m_info.size() == infoSize()*popSize(), SystemError,
			"Info size is wrong");

		UINT os = infoSize();
		UINT idx;
		// if this field exists, return directly
		try
		{
			idx = infoIdx(field);
			// only needs to initialize
			int oldAncPop = m_curAncestralPop;
			for(UINT anc=0; anc <= m_ancestralPops.size(); anc++)
			{
				useAncestralPop(anc);
				for(IndIterator ind=indBegin(); ind!=indEnd(); ++ind)
					ind->setInfo(init, idx);
			}
			useAncestralPop(oldAncPop);
			return;
		}
		catch(IndexError &)
		{
			// and we only add field to node 0
			struAddInfoField(field);
		}

		/// adjust information size.
		UINT is = infoSize();
		if(os != is)
		{
			int oldAncPop = m_curAncestralPop;
			for(UINT anc=0; anc <= m_ancestralPops.size(); anc++)
			{
				useAncestralPop(anc);
				vectorinfo newInfo(is*popSize());
				/// copy the old stuff in
				InfoIterator ptr = newInfo.begin();
				for(IndIterator ind=indBegin(); ind!=indEnd(); ++ind)
				{
					copy(ind->infoBegin(), ind->infoBegin() + is - 1, ptr);
					ind->setInfoPtr(ptr);
					fill(ind->infoBegin() + os, ind->infoEnd(), init);
					ptr += is;
				}
				m_info.swap(newInfo);
			}
			useAncestralPop(oldAncPop);
		}
		return;
#endif

	}

	void population::addInfoFields(const vectorstr& fields, double init)
	{
#ifdef SIMUMPI
		if(mpiRank()==0)
		{
			DBG_ASSERT(m_info.size() == localInfoSize()*popSize(), SystemError,
				"Info size is wrong");
			// oldsize, this is valid for rank 0
			UINT os = localInfoSize();
			for(vectorstr::const_iterator it=fields.begin(); it!=fields.end(); ++it)
			{
				try
				{
					// has field
					UINT idx = localInfoIdx(*it);
					// only needs to initialize
					int oldAncPop = m_curAncestralPop;
					for(UINT anc=0; anc <= m_ancestralPops.size(); anc++)
					{
						useAncestralPop(anc);

						for(IndIterator ind=indBegin(); ind!=indEnd(); ++ind)
							ind->setInfo(init, idx);
					}
					useAncestralPop(oldAncPop);
				}
				catch(IndexError &)
				{
					struAddInfoField(*it);
				}
			}

			// add these fields
			/// adjust information size.
			UINT is = localInfoSize();
			// need to extend.
			if(is != os)
			{
				int oldAncPop = m_curAncestralPop;
				for(UINT anc=0; anc <= m_ancestralPops.size(); anc++)
				{
					useAncestralPop(anc);
					vectorinfo newInfo(is*popSize(), 0.);
					/// copy the old stuff in
					InfoIterator ptr = newInfo.begin();
					for(IndIterator ind=indBegin(); ind!=indEnd(); ++ind)
					{
						copy(ind->infoBegin(), ind->infoBegin() + os, ptr);
						ind->setInfoPtr(ptr);
						fill(ind->infoBegin() + os, ind->infoEnd(), init);
						ptr += is;
					}
					m_info.swap(newInfo);
				}
				useAncestralPop(oldAncPop);
			}
		}
#else
		DBG_ASSERT(m_info.size() == infoSize()*popSize(), SystemError,
			"Info size is wrong");
		// oldsize, this is valid for rank 0
		UINT os = infoSize();
		for(vectorstr::const_iterator it=fields.begin(); it!=fields.end(); ++it)
		{
			try
			{
				// has field
				UINT idx = infoIdx(*it);
				// only needs to initialize
				int oldAncPop = m_curAncestralPop;
				for(UINT anc=0; anc <= m_ancestralPops.size(); anc++)
				{
					useAncestralPop(anc);

					for(IndIterator ind=indBegin(); ind!=indEnd(); ++ind)
						ind->setInfo(init, idx);
				}
				useAncestralPop(oldAncPop);
			}
			catch(IndexError &)
			{
				struAddInfoField(*it);
			}
		}

		// add these fields
		/// adjust information size.
		UINT is = infoSize();
		// need to extend.
		if(is != os)
		{
			int oldAncPop = m_curAncestralPop;
			for(UINT anc=0; anc <= m_ancestralPops.size(); anc++)
			{
				useAncestralPop(anc);
				vectorinfo newInfo(is*popSize(), 0.);
				/// copy the old stuff in
				InfoIterator ptr = newInfo.begin();
				for(IndIterator ind=indBegin(); ind!=indEnd(); ++ind)
				{
					copy(ind->infoBegin(), ind->infoBegin() + os, ptr);
					ind->setInfoPtr(ptr);
					fill(ind->infoBegin() + os, ind->infoEnd(), init);
					ptr += is;
				}
				m_info.swap(newInfo);
			}
			useAncestralPop(oldAncPop);
		}
#endif

	}

	/// set fields
	void population::setInfoFields(const vectorstr& fields, double init)
	{
#ifdef SIMUMPI
		if(mpiRank()==0)
		{
#endif
			struSetInfoFields(fields);

			/// reset info vector
			int oldAncPop = m_curAncestralPop;
#ifdef SIMUMPI
			UINT is = localInfoSize();
#else
			UINT is = infoSize();
#endif
			for(UINT anc=0; anc <= m_ancestralPops.size(); anc++)
			{
				useAncestralPop(anc);
				vectorinfo newInfo(is*popSize(), init);
				InfoIterator ptr = newInfo.begin();
				for(IndIterator ind=indBegin(); ind!=indEnd(); ++ind, ptr += is)
					ind->setInfoPtr(ptr);
				m_info.swap(newInfo);
			}
			useAncestralPop(oldAncPop);
#ifdef SIMUMPI
		}
#endif

	}

	/// set ancestral depth, can be -1
	void population::setAncestralDepth(int depth)
	{
		// just to make sure.
		useAncestralPop(0);
		//
		if(depth >=0 && m_ancestralPops.size() > static_cast<size_t>(depth))
		{
			int numRemove = m_ancestralPops.size() - depth;
			while( numRemove-- > 0 )
				m_ancestralPops.pop_back();
		}
		DBG_ASSERT( depth <0 || m_ancestralPops.size() <= static_cast<size_t>(depth), SystemError,
			"Failed to change ancestral Depth");

		m_ancestralDepth = depth;
	}

	void population::useAncestralPop(UINT idx)
	{
		if( m_curAncestralPop >= 0 && idx == static_cast<UINT>(m_curAncestralPop))
			return;

		DBG_DO(DBG_POPULATION, cout << "Use ancestralPop: " << idx <<
			"Curidx: " <<  m_curAncestralPop << endl);

		if( idx == 0 || m_curAncestralPop != 0)	  // recover pop.
		{
			popData& pd = m_ancestralPops[ m_curAncestralPop-1 ];
			pd.m_subPopSize.swap(m_subPopSize);
			pd.m_genotype.swap(m_genotype);
			pd.m_info.swap(m_info);
			pd.m_inds.swap(m_inds);
			m_curAncestralPop = 0;
#ifndef OPTIMIZED
			//DBG_FAILIF( pd.m_startingGenoPtr != m_genotype.begin(),
			//	SystemError, "Starting genoptr has been changed.");
			pd.m_startingGenoPtr = pd.m_genotype.begin();
#endif
			if( idx == 0)
			{									  // restore key paraemeters from data
				m_popSize = m_inds.size();
				m_numSubPop = m_subPopSize.size();
				m_subPopIndex.resize( m_numSubPop + 1);
				// build subPop index
				UINT i = 1;
				for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
					m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];

				return;
			}
		}

		// now m_curAncestralPop is zero.
		DBG_ASSERT( idx <= m_ancestralPops.size(),
			ValueError, "Ancestry population " + toStr(idx) + " does not exist.");

		// now idx should be at least 1
		m_curAncestralPop = idx;
		// swap  1 ==> 0, 2 ==> 1

		popData& pd = m_ancestralPops[ m_curAncestralPop -1];
		pd.m_subPopSize.swap(m_subPopSize);
		pd.m_genotype.swap(m_genotype);
		pd.m_info.swap(m_info);
		pd.m_inds.swap(m_inds);
#ifndef OPTIMIZED
		pd.m_startingGenoPtr = pd.m_genotype.begin();
#endif
		// use pd
		m_popSize = m_inds.size();
		m_numSubPop = m_subPopSize.size();
		m_subPopIndex.resize( m_numSubPop + 1);
		UINT i = 1;
		for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
			m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];
	}

	void population::savePopulation(const string& filename, const string& format, bool compress) const
	{
		io::filtering_ostream ofs;
		// get file extension
		string ext = fileExtension(filename);
#ifndef DISABLE_COMPRESSION
		if(compress || (ext.size() > 3 && ext.substr(ext.size() - 3, 3) == ".gz"))
			ofs.push(io::gzip_compressor());
#endif
		ofs.push(io::file_sink(filename));

		if(!ofs)
			throw ValueError("Can not open file " + filename );

		if( format == "text" || (format == "auto" && (ext == "txt" || ext == "txt.gz" )))
		{
			boost::archive::text_oarchive oa(ofs);
			oa << *this;
		}
		else if (format == "xml" || (format == "auto" && (ext == "xml" || ext == "xml.gz" )))
		{
			boost::archive::xml_oarchive oa(ofs);
			oa << boost::serialization::make_nvp("population",*this);
		}
		else if (format == "bin" ||  (format == "auto" && (ext == "bin" || ext == "bin.gz" )))
		{
			boost::archive::binary_oarchive oa(ofs);
			oa << *this;
		}
		else
			throw ValueError("Wrong format type. Use one of text, xml, bin or appropriate extension txt, xml or bin");
	}

	void population::loadPopulation(const string& filename, const string& format)
	{
		io::filtering_istream ifs;
		bool gzipped = isGzipped(filename);
		if(gzipped)
#ifdef DISABLE_COMPRESSION
			throw ValueError("This version of simuPOP can not handle compressed file");
#else
		ifs.push(io::gzip_decompressor());
#endif
		ifs.push(io::file_source(filename));
		// do not have to test again.
		if(!ifs)
			throw ValueError("Can not open file " + filename );

		// get file extension
		string ext = fileExtension(filename);

		// try to load the file, according to file extension.
		try
		{
			if( format == "text" || (format == "auto" && (ext == "txt" || ext == "txt.gz" ) ))
			{
				boost::archive::text_iarchive ia(ifs);
				ia >> *this;
			}
			else if (format == "xml" ||  (format == "auto" && (ext == "xml" || ext == "xml.gz" ) ))
			{
				boost::archive::xml_iarchive ia(ifs);
				ia >> boost::serialization::make_nvp("population",*this);
			}
			else if (format == "bin" || (format == "auto" && (ext == "bin" || ext == "bin.gz" ) ))
			{
				boost::archive::binary_iarchive ia(ifs);
				ia >> *this;
			}
			else								  // need special handling
				throw;
		}
		catch(...)								  // if any error happens, or can not determine format, try different methods
		{
			// first close the file handle.

			DBG_DO(DBG_POPULATION,
				cout << "Can not determine file type, or file type is wrong. Trying different ways." << endl);

			// open a fresh ifstream
			io::filtering_istream ifbin;
			if(gzipped)
				ifbin.push(io::gzip_decompressor());
			ifbin.push(io::file_source(filename));

			// try to load the file using different iarchives.
			try									  // binary?
			{
				boost::archive::binary_iarchive ia(ifbin);
				ia >> *this;
			}
			catch(...)							  // not binary, text?
			{
				io::filtering_istream iftxt;
				if(gzipped)
					iftxt.push(io::gzip_decompressor());
				iftxt.push(io::file_source(filename));
				try
				{
					boost::archive::text_iarchive ia(iftxt);
					ia >> *this;
				}
				catch(...)						  // then xml?
				{
					io::filtering_istream ifxml;
					if(gzipped)
						ifxml.push(io::gzip_decompressor());
					ifxml.push(io::file_source(filename));
					try
					{
						boost::archive::xml_iarchive ia(ifxml);
						ia >> boost::serialization::make_nvp("population",*this);
					}
					catch(...)
					{
						throw ValueError("Failed to load population " + filename + " in " + format + " format.\n");
					}
				}								  // try xml
			}									  // try text
		}										  // try bin
	}

	PyObject* population::vars(int subPop)
	{
		if(subPop < 0)
		{
			Py_INCREF(m_vars.dict());
			return m_vars.dict();
		}
		else
		{
			DBG_ASSERT( static_cast<UINT>(subPop) < numSubPop() ,
				IndexError, "Subpop index out of range of 0 ~ " + toStr(numSubPop()-1) );

			DBG_ASSERT( hasVar("subPop"), ValueError,
				"subPop statistics does not exist yet.");

			PyObject* spObj = getVar("subPop");
			spObj = PyList_GetItem(spObj, subPop);

			DBG_ASSERT(spObj != NULL, SystemError,
				"Something is wrong about the length of subPop list. ");

			Py_INCREF( spObj);
			return spObj;
		}
	}

	/// CPPONLY
	/// The same as vars(), but without increasing
	/// reference count.
	PyObject* population::dict(int subPop)
	{
		if(subPop < 0)
			return m_vars.dict();
		else
		{
			DBG_ASSERT( static_cast<UINT>(subPop) < numSubPop() ,
				IndexError, "Subpop index out of range of 0 ~ " + toStr(numSubPop()-1) );

			DBG_ASSERT( hasVar("subPop"), ValueError,
				"subPop statistics does not exist yet.");

			PyObject* spObj =  getVar("subPop");
			spObj = PyList_GetItem(spObj, subPop);

			DBG_ASSERT(spObj != NULL, SystemError,
				"Something is wrong about the length of subPop list. ");

			return spObj;
		}
	}

	/// CPPONLY
	void population::adjustGenoPosition(bool order)
	{
		DBG_DO(DBG_POPULATION, cout << "Adjust geno position " << endl);

		// everyone in strict order
		if(order)
		{
			DBG_DO(DBG_POPULATION, cout << "Refresh all order " << endl);
#ifdef SIMUMPI
			vectora tmpGenotype(m_popSize*localGenoSize());
			size_t sz = localGenoSize();
			vectorinfo tmpInfo(m_popSize*localInfoSize());
			UINT is = localInfoSize();
#else
			vectora tmpGenotype(m_popSize*genoSize());
			size_t sz = genoSize();
			vectorinfo tmpInfo(m_popSize*infoSize());
			UINT is = infoSize();
#endif
			vectora::iterator it = tmpGenotype.begin();
			vectorinfo::iterator infoPtr = tmpInfo.begin();

			for(IndIterator ind=indBegin(), indEd=indEnd(); ind!=indEd; ++ind)
			{
#ifdef BINARYALLELE
				copyGenotype(ind->genoBegin(), it, sz);
#else
				copy(ind->genoBegin(), ind->genoEnd(), it);
#endif
				ind->setGenoPtr(it);
				copy(ind->infoBegin(), ind->infoEnd(), infoPtr);
				ind->setInfoPtr(infoPtr);
				it += sz;
				infoPtr += is;
				ind->setShallowCopied(false);
			}
			// discard original genotype
			tmpGenotype.swap(m_genotype);
			tmpInfo.swap(m_info);
			// set geno pointer
			setShallowCopied(false);
			setInfoOrdered(true);
			return;
		}

		/// find out how many individuals are shallow copied.
		vectorl scIndex(0);
		ULONG j, k=0, jEnd;

		for(UINT sp=0, spEd = numSubPop(); sp < spEd;  sp++)
		{
#ifdef SIMUMPI
			GenoIterator spBegin = m_genotype.begin() + m_subPopIndex[sp]*localGenoSize();
			GenoIterator spEnd   = m_genotype.begin() + m_subPopIndex[sp+1]*localGenoSize();
#else
			GenoIterator spBegin = m_genotype.begin() + m_subPopIndex[sp]*genoSize();
			GenoIterator spEnd   = m_genotype.begin() + m_subPopIndex[sp+1]*genoSize();
#endif
			for(j=0, jEnd = subPopSize(sp); j < jEnd;  j++)
			{
				if(m_inds[k].shallowCopied() )
				{
					if(indGenoBegin(k) < spBegin || indGenoEnd(k) > spEnd)
						/// record individual index and genoPtr
						scIndex.push_back(k);
					else
						m_inds[k].setShallowCopied(false);
				}
				k++;
			}
		}

		if( scIndex.empty())
		{
			//setShallowCopied(false);
			return;
		}

		/// to further save time, deal with a special case that there are
		/// only two shallowCopied individuals
		if( scIndex.size() == 2 )
		{
			// swap!
			GenoIterator tmp = m_inds[ scIndex[0] ].genoPtr();
			m_inds[scIndex[0] ].setGenoPtr( m_inds[ scIndex[1] ].genoPtr() );
			m_inds[scIndex[1] ].setGenoPtr( tmp);

			Allele tmp1;
			for(UINT a =0; a < genoSize(); ++a)
			{
				tmp1 =  m_inds[ scIndex[0] ].allele(a);
				m_inds[scIndex[0] ].setAllele(m_inds[ scIndex[1] ].allele(a), a);
				m_inds[scIndex[1] ].setAllele(tmp1, a);
			}

			// copy info
			InfoType tmp2;
#ifdef SIMUMPI
			for(UINT a=0; a < localInfoSize(); ++a)
#else
				for(UINT a=0; a < infoSize(); ++a)
#endif
			{
				tmp2 =  m_inds[ scIndex[0] ].info(a);
				m_inds[scIndex[0] ].setInfo(m_inds[ scIndex[1] ].info(a), a);
				m_inds[scIndex[1] ].setInfo(tmp2, a);
			}

			m_inds[scIndex[0] ].setShallowCopied(false);
			m_inds[scIndex[1] ].setShallowCopied(false);
			setShallowCopied(false);
			setInfoOrdered(true);
			return;
		}

		/// save genotypic info
#ifdef SIMUMPI
		vectora scGeno(scIndex.size() * localNumLoci() * ploidy());
		vectorinfo scInfo(scIndex.size() * localInfoSize());
#else
		vectora scGeno(scIndex.size() * totNumLoci() * ploidy());
		vectorinfo scInfo(scIndex.size() * infoSize());
#endif
		vector<GenoIterator> scPtr( scIndex.size() );
		vector<InfoIterator> scInfoPtr( scIndex.size() );

		size_t i, iEnd;

		for(i=0, iEnd = scIndex.size(); i < iEnd;  i++)
		{
			scPtr[i] = m_inds[ scIndex[i]].genoPtr();
#ifdef SIMUMPI
#ifdef BINARYALLELE
			copyGenotype(indGenoBegin(scIndex[i]), scGeno.begin() + i* localGenoSize(), localGenoSize());
#else
			copy(indGenoBegin(scIndex[i]), indGenoEnd(scIndex[i]), scGeno.begin() + i* localGenoSize());
#endif
			scInfoPtr[i] = m_inds[ scIndex[i]].infoPtr();
			copy(ind(scIndex[i]).infoBegin(), ind(scIndex[i]).infoEnd(),
				scInfo.begin() + i* localInfoSize());
#else
#ifdef BINARYALLELE
			copyGenotype(indGenoBegin(scIndex[i]), scGeno.begin() + i* genoSize(), genoSize());
#else
			copy(indGenoBegin(scIndex[i]), indGenoEnd(scIndex[i]), scGeno.begin() + i* genoSize());
#endif
			scInfoPtr[i] = m_inds[ scIndex[i]].infoPtr();
			copy(ind(scIndex[i]).infoBegin(), ind(scIndex[i]).infoEnd(),
				scInfo.begin() + i* infoSize());
#endif
		}

		DBG_DO(DBG_POPULATION, cout << "Shallow copied" << scIndex << endl);

		/// sort the pointers!
		sort( scPtr.begin(), scPtr.end());
		sort( scInfoPtr.begin(), scInfoPtr.end());

		/// copy back.
		for(i=0, iEnd =scIndex.size(); i < iEnd;  i++)
		{
			m_inds[ scIndex[i] ].setGenoPtr( scPtr[i]);
#ifdef SIMUMPI
#ifdef BINARYALLELE
			copyGenotype(scGeno.begin() + i*localGenoSize(), indGenoBegin( scIndex[i] ), localGenoSize());
#else
			copy( scGeno.begin() + i*  localGenoSize(), scGeno.begin() + (i+1)*  localGenoSize(),
				indGenoBegin( scIndex[i] ));
#endif
			m_inds[ scIndex[i] ].setInfoPtr( scInfoPtr[i]);
			copy( scInfo.begin() + i*  localInfoSize(), scInfo.begin() + (i+1)*  localInfoSize(),
				ind(scIndex[i]).infoBegin());
#else
#ifdef BINARYALLELE
			copyGenotype(scGeno.begin() + i*genoSize(), indGenoBegin( scIndex[i] ), genoSize());
#else
			copy( scGeno.begin() + i*  genoSize(), scGeno.begin() + (i+1)*  genoSize(),
				indGenoBegin( scIndex[i] ));
#endif
			m_inds[ scIndex[i] ].setInfoPtr( scInfoPtr[i]);
			copy( scInfo.begin() + i*  infoSize(), scInfo.begin() + (i+1)*  infoSize(),
				ind(scIndex[i]).infoBegin());
#endif
			m_inds[ scIndex[i] ].setShallowCopied(false);
		}
		//setShallowCopied(false);
		return;
	}

	/// CPPONLY
	void population::adjustInfoPosition(bool order)
	{
		DBG_DO(DBG_POPULATION, cout << "Adjust info position " << endl);

		// everyone in strict order
		/*
		if(order)
		{
		*/
		DBG_DO(DBG_POPULATION, cout << "Refresh all order " << endl);
#ifdef SIMUMPI
		UINT is = localInfoSize();
#else
		UINT is = infoSize();
#endif
		size_t i;
		vectorinfo tmpInfo(m_popSize*is);
		vectorinfo::iterator infoPtr = tmpInfo.begin();
		vectorinfo::iterator tmp;

		for(IndIterator ind=indBegin(), indEd=indEnd(); ind!=indEd; ++ind)
		{
			tmp = ind->infoBegin();
			for(i=0; i<is; ++i)
				infoPtr[i] = tmp[i];
			ind->setInfoPtr(infoPtr);
			infoPtr += is;
		}
		// discard original genotype
		m_info.swap(tmpInfo);
		setInfoOrdered(true);
		return;
	}

	population& LoadPopulation(const string& file, const string& format)
	{
#ifndef _NO_SERIALIZATION_
		population *p = new population(1);
		p->loadPopulation(file, format);
		return *p;
#else
		cout << "This feature is not supported in this platform" << endl;
		return *new population(1);
#endif
	}

	vectorf testGetinfoFromInd(population& pop)
	{
		vectorf a(pop.popSize());
		size_t i=0;
		for(population::IndIterator ind=pop.indBegin(), indEnd = pop.indEnd();
			ind != indEnd; ++ind)
		a[i++] = ind->info(0);
		return a;
	}

	vectorf testGetinfoFromPop(population& pop, bool order)
	{
		vectorf a(pop.popSize());
		size_t i=0;

		if(order)
			pop.adjustInfoPosition(true);
		for(GappedInfoIterator it=pop.infoBegin(0, true),
			itEnd=pop.infoEnd(0, true); it != itEnd; ++it)
			a[i++] = *it;
		return a;
	}
}
