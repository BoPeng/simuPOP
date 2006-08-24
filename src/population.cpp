/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu                                                        *
 *                                                                         *
 *   $LastChangedDate: 2006-02-21 15:27:25 -0600 (Tue, 21 Feb 2006) $
 *   $Rev: 191 $
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

namespace io = boost::iostreams;

namespace simuPOP
{
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
		const vectorstr& infoFields)
		:
	GenoStruTrait(),
		m_popSize(size),
		m_numSubPop(subPop.size()),
		m_subPopSize(subPop),
		m_popGenoSize(0),
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
		m_fitness(0),
		m_shallowCopied(false)
	{
		DBG_FAILIF(maxAllele > MaxAllele, ValueError,
			"maxAllele is bigger than maximum allowed allele state of this library (" + toStr(MaxAllele) +
			")\nPlease use simuOpt.setOptions(longAllele=True) to use the long allele version of simuPOP.");

		DBG_FAILIF(maxAllele == 0, ValueError,
			"maxAllele should be at least 1 (0,1 two states). ");

		DBG_DO( DBG_POPULATION, cout << "Constructor of population is called\n");

		// if specify subPop but not m_popSize
		if( !subPop.empty() )
		{
			if( size == 0 )
				m_popSize = accumulate(subPop.begin(), subPop.end(), 0UL);
			else
				DBG_ASSERT( m_popSize == accumulate(subPop.begin(), subPop.end(), 0UL),
					ValueError, "If both size and subPop are specified, size should equal to sum(subPop)");
		}

		// get a GenoStructure with parameters. GenoStructure may be shared by some populations
		// a whole set of functions ploidy() etc in GenoStruTriat can be used after this step.
		this->setGenoStructure(ploidy, loci, sexChrom, lociPos, alleleNames, lociNames, maxAllele, infoFields);

		DBG_DO( DBG_DEVEL, cout << "individual size is " << sizeof(individual) << '+'
			<< sizeof(Allele) << '*' << genoSize() << endl
			<< "Info: " << sizeof(SubPop_ID) << ", Tag: " << sizeof(TagType)
			<< ", GenoPtr: " << sizeof(Allele*) << ", Flag: " << sizeof(unsigned char)
			<< ", plus genoStru" << endl );

		// size of genotypic data for the whole population
		m_popGenoSize = totNumLoci() * ploidy * m_popSize;

		try
		{
			// allocate memory here (not in function definition)
			m_inds.resize(m_popSize);

			// create genotype vector holding alleles for all individuals.
			m_genotype.resize( m_popGenoSize);

			/// allocate info
			size_t is = infoSize();
			DBG_ASSERT(is == infoFields.size(), SystemError, "Wrong geno structure");
			m_info.resize(m_popSize*is);

			// set subpopulation indexes, do not allow popsize change
			setSubPopStru(subPop, false);

			// set individual pointers
			// reset individual pointers
			GenoIterator ptr = m_genotype.begin();
			InfoIterator infoPtr = m_info.begin();
			UINT step = genoSize();
			for(ULONG i=0; i< m_popSize; ++i, ptr+=step, infoPtr+=is)
			{
				m_inds[i].setGenoPtr(ptr);
				m_inds[i].setGenoStruIdx(this->genoStruIdx());
				m_inds[i].setShallowCopied(false);
				m_inds[i].setInfoPtr(infoPtr);
			}
		}
		catch(...)
		{
			cout << "You are creating a population of size " << m_popSize << endl
				<< " which requires approximately " << static_cast<double>(m_popGenoSize)/1024/1024
				+ static_cast<double>(m_popSize) * sizeof( individual ) /1024/1024 << "M RAM." << endl;
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
		m_popGenoSize(rhs.m_popGenoSize),
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
		m_fitness(0),							  // do not copy fitness
		m_shallowCopied(false)
	{
		DBG_DO(DBG_POPULATION,
			cout << "Copy constructor of population is called\n" << endl);

		try
		{
			m_inds.resize(rhs.m_popSize);
			m_genotype.resize(rhs.m_popGenoSize);
			m_info.resize(rhs.m_popSize*infoSize());
		}
		catch(...)
		{
			cout << "You are creating a population of size " << m_popSize << endl
				<< " which requires approximately " << static_cast<double>(m_popGenoSize)/1024/1024
				+ static_cast<double>(m_popSize) * sizeof( individual ) /1024/1024 << "M RAM." << endl;
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
		UINT step = this->genoSize();
		UINT infoStep = this->infoSize();
		for(ULONG i=0; i< m_popSize; ++i, ptr+=step, infoPtr+=infoStep)
		{
			m_inds[i].setGenoPtr( ptr );
			m_inds[i].setInfoPtr( infoPtr );
			m_inds[i].copyFrom( rhs.m_inds[i]);
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

	int population::__cmp__(const population& rhs) const
	{
		if( genoStruIdx() != rhs.genoStruIdx() )
			return 1;

		if(popSize() != rhs.popSize() )
			return 1;

		for( ULONG i=0, iEnd = popSize(); i < iEnd; ++i)
			if( m_inds[i] != rhs.m_inds[i])
				return 1;

		// FIXME: also compare ancestral populations
		return 0;
	}

	void population::setSubPopStru(const vectorlu& newSubPopSizes, bool allowPopSizeChange)
	{
		if( allowPopSizeChange && !m_fitness.empty() )
			throw SystemError("individual order can not be changed with non-empty fitness vector\n"
				"Please put selector after migrator or other such operators.");

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
				m_popGenoSize = genoSize() * m_popSize;

				try
				{
					m_genotype.resize( m_popGenoSize );
					m_inds.resize(m_popSize);
					m_info.resize(m_popSize*infoSize());
				}
				catch(...)
				{
					cout << "You are creating a population of size " << m_popSize << endl
						<< " which requires approximately " << static_cast<double>(m_popGenoSize)/1024/1024
						+ static_cast<double>(m_popSize) * sizeof( individual ) /1024/1024 << "M RAM." << endl;
					cout << "Memory allocation fail. " << endl;
					throw OutOfMemory("Memory allocation fail");
				}
				// reset individual pointers
				GenoIterator ptr = m_genotype.begin();
				InfoIterator infoPtr = m_info.begin();
				UINT step = genoSize();
				UINT is = infoSize();
				for(ULONG i=0; i< m_popSize; ++i, ptr+=step, infoPtr+=is)
				{
					m_inds[i].setGenoPtr( ptr );
					m_inds[i].setInfoPtr( infoPtr );
					m_inds[i].setGenoStruIdx(genoStruIdx());
					m_inds[i].setShallowCopied(false);
				}
				m_shallowCopied=false;
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
			ULONG newPopGenoSize = genoSize() * newPopSize;
			vectora newGenotype(newPopGenoSize);
			vector<InfoType> newInfo(newPopSize*infoSize());
			vector<individual> newInds(newPopSize);

			DBG_ASSERT( indEnd()== newPopSize+it, SystemError,
				"Pointer misplaced. ");

			// assign genotype location and set structure information for individuals
			GenoIterator ptr = newGenotype.begin();
			InfoIterator infoPtr = newInfo.begin();
			UINT step = genoSize();
			UINT infoStep = infoSize();
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
			m_popGenoSize = newPopGenoSize;
			setShallowCopied(false);
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
		if( ! m_fitness.empty() )
			throw SystemError("individual order can not be changed with non-empty fitness vector\n"
				"Please put selector after migrator or other such operators.");

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
					ind->setSubPopID(-1);			  // remove
			}
			// other subpop shift left
			else if(shiftSubPopID)
			{
				for(IndIterator ind = indBegin(sp); ind != indEnd(sp); ++ind)
					ind->setSubPopID(sp-shift);		  // shift left
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
		if( ! m_fitness.empty() )
			throw SystemError("individual order can not be changed with non-empty fitness vector\n"
				"Please put selector after migrator or other such operators.");

		setIndSubPopIDWithID();
		if( subPop == -1 )
		{
			for(size_t i = 0; i < inds.size(); ++i)
				ind(inds[i]).setSubPopID(-1);		  // remove
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

	void population::reorderSubPops(const vectoru& order, const vectoru& rank,
		bool removeEmptySubPops)
	{
		if( ! m_fitness.empty() )
			throw SystemError("individual order can not be changed with non-empty fitness vector\n"
				"Please put selector after migrator or other such operators.");

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

	/** form a new population according to info, info can be given directly */
	population& population::newPopByIndID(bool keepAncestralPops,
		const vectori& id, bool removeEmptySubPops)
	{
		// copy the population over (info is also copied)
		population& pop = this->clone(keepAncestralPops);
		// and shrink them
		pop.setSubPopByIndID(id);
		if( removeEmptySubPops)
			pop.removeEmptySubPops();
		return pop;
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
		UINT oldTotNumLoci = this->totNumLoci();

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

		// new geno structure is in effective now!
		this->setGenoStructure(this->ploidy(), newNumLoci, this->sexChrom(), newLociDist,
			this->alleleNames(), newLociNames, this->maxAllele(), this->infoFields() );
		// prepare data
		//
		// keep m_popSize;
		// keep m_numSubPop;
		// keep m_subPopSize;
		// keep m_subPopIndex;

		// genotype
		// allocate new genotype and inds
		ULONG newPopGenoSize = genoSize() * m_popSize;
		vectora newGenotype(newPopGenoSize);
		// keep newInds();

		// copy data over
		GenoIterator ptr = newGenotype.begin();
		UINT pEnd = this->ploidy();
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
		m_popGenoSize = newPopGenoSize;
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
		setShallowCopied(false);
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
		m_subPopSize.swap( rhs.m_subPopSize);
		m_popGenoSize = rhs.m_popGenoSize;
		m_subPopIndex.swap( rhs.m_subPopIndex);
		m_genotype.swap( rhs.m_genotype);
		m_info.swap( rhs.m_info);
		m_inds.swap(rhs.m_inds);
		m_fitness.swap(rhs.m_fitness);
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
			rhs.m_popGenoSize = genoSize()*rhs.m_popSize;
			rhs.m_subPopIndex.resize(2,0);
			rhs.m_subPopIndex[1] = rhs.m_popSize;
			// no need to set genoPtr or genoStru()
		}
	}

	/// request info field
	/// code: In this version, I do not allow dynamic addition
	/// of information field
	/*
	int population::requestInfoField(const string name)
	{
		// if this field does not yet exist, create a new one
		int index = infoIdx(name);
		if (index >= 0)
			return index;
			
		/// the code does not exist yet,
		/// crete a new one
		index = addInfoField(name);
		///
		vector<InfoType> newInfo((index+1)*popSize());
		/// copy the old stuff in
		InfoIterator ptr = newInfo.begin();
		for(IndIterator ind=indBegin(); ind!=indEnd(); ++ind)
		{
			if(index>0)
				copy(ind->infoBegin(), ind->infoBegin() + index, ptr);
			ind->setInfoPtr(ptr);
			ptr += index + 1;
		}
		m_info.swap(newInfo);
		return index;
	}
	*/

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

	void population::useAncestralPop(int idx)
	{
		if( idx == m_curAncestralPop)
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
			DBG_FAILIF( pd.m_startingGenoPtr != m_genotype.begin(),
				SystemError, "Starting genoptr has been changed.");
			pd.m_startingGenoPtr = pd.m_genotype.begin();
#endif
			if( idx == 0)
			{									  // restore key paraemeters from data
				m_popSize = m_inds.size();
				m_popGenoSize = genoSize() * m_popSize;
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
		DBG_ASSERT( idx >0 && static_cast<size_t>(idx-1) < m_ancestralPops.size(),
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
		m_popGenoSize = genoSize() * m_popSize;
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
		if(compress || (ext.size() > 3 && ext.substr(ext.size() - 3, 3) == ".gz"))
			ofs.push(io::gzip_compressor());
		ofs.push(io::file_sink(filename));

		if(!ofs)
			throw ValueError("Can not open file " + filename );

		if( format == "text" || (format == "auto" && (ext == "txt" || ext == "txt.gz" )))
		{
			boost::archive::text_oarchive oa(ofs);
			oa << *this;
		}
#ifndef __NO_XML_SUPPORT__
		else if (format == "xml" || (format == "auto" && (ext == "xml" || ext == "xml.gz" )))
		{
			boost::archive::xml_oarchive oa(ofs);
			oa << boost::serialization::make_nvp("population",*this);
		}
#endif
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
			ifs.push(io::gzip_decompressor());
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
#ifndef __NO_XML_SUPPORT__
			else if (format == "xml" ||  (format == "auto" && (ext == "xml" || ext == "xml.gz" ) ))
			{
				boost::archive::xml_iarchive ia(ifs);
				ia >> boost::serialization::make_nvp("population",*this);
			}
#endif
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
#ifndef __NO_XML_SUPPORT__
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
						throw ValueError("Failed to load population. Your file may be corrupted, "
							"or being a copy of non-transferrable file (.bin)");
					}
#else
					throw ValueError("Failed to load population. Your file may be corrupted, "
						"or being a copy of non-transferrable file (.bin)");
#endif
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
	void population::adjustGenoPosition(bool deep)
	{
		DBG_DO(DBG_POPULATION, cout << "Adjust geno position " << endl);

		// everyone in strict order
		if(deep)
		{
			vectora tmpGenotype(m_popGenoSize);
			vector<InfoType> tmpInfo(m_popSize*infoSize());
			vectora::iterator it = tmpGenotype.begin();
			vector<InfoType>::iterator infoPtr = tmpInfo.begin();
			UINT is = infoSize();

			for(IndIterator ind=indBegin(), indEd=indEnd(); ind!=indEd; ++ind)
			{
				copy(ind->genoBegin(), ind->genoEnd(), it);
				ind->setGenoPtr(it);
				copy(ind->infoBegin(), ind->infoEnd(), infoPtr);
				ind->setInfoPtr(infoPtr);
				it += this->totNumLoci()*this->ploidy();
				infoPtr += is;
				ind->setShallowCopied(false);
			}
			// discard original genotype
			tmpGenotype.swap(m_genotype);
			tmpInfo.swap(m_info);
			// set geno pointer
			setShallowCopied(false);
			return;
		}

		/// find out how many individuals are shallow copied.
		vectorl scIndex(0);
		ULONG j, k=0, jEnd;

		for(UINT sp=0, spEd = numSubPop(); sp < spEd;  sp++)
		{
			GenoIterator spBegin = m_genotype.begin() + m_subPopIndex[sp]*genoSize();
			GenoIterator spEnd   = m_genotype.begin() + m_subPopIndex[sp+1]*genoSize();
			for( j=0, jEnd = subPopSize(sp); j < jEnd;  j++)
			{
				if( m_inds[k].shallowCopied() )
				{
					if( indGenoBegin(k) < spBegin
						|| indGenoEnd(k) > spEnd )
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
			setShallowCopied(false);
			return;
		}

		/// to further save time, deal with a special case that there are
		/// only two shallowCopied individuals
		if( scIndex.size() == 2 )
		{
			// swap!
			GenoIterator tmp = m_inds[ scIndex[0] ].genoPtr();
			m_inds[ scIndex[0] ].setGenoPtr( m_inds[ scIndex[1] ].genoPtr() );
			m_inds[ scIndex[1] ].setGenoPtr( tmp);

			Allele tmp1;
			for(UINT a =0; a < genoSize(); ++a)
			{
				tmp1 =  m_inds[ scIndex[0] ].allele(a);
				m_inds[ scIndex[0] ].setAllele(m_inds[ scIndex[1] ].allele(a), a);
				m_inds[ scIndex[1] ].setAllele(tmp1, a);
			}

			// copy info
			InfoType tmp2;
			for(UINT a=0; a < infoSize(); ++a)
			{
				tmp2 =  m_inds[ scIndex[0] ].info(a);
				m_inds[ scIndex[0] ].setInfo(m_inds[ scIndex[1] ].info(a), a);
				m_inds[ scIndex[1] ].setInfo(tmp2, a);
			}
			
			m_inds[ scIndex[0] ].setShallowCopied(false);
			m_inds[ scIndex[1] ].setShallowCopied(false);
			setShallowCopied(false);
			return;
		}

		/// save genotypic info
		vectora scGeno(scIndex.size() * totNumLoci() * ploidy());
		vector<GenoIterator> scPtr( scIndex.size() );
		vector<InfoType> scInfo(scIndex.size() * infoSize());
		vector<InfoIterator> scInfoPtr( scIndex.size() );

		size_t i, iEnd;

		for(i=0, iEnd = scIndex.size(); i < iEnd;  i++)
		{
			scPtr[i] = m_inds[ scIndex[i]].genoPtr();
			copy( indGenoBegin(scIndex[i]), indGenoEnd(scIndex[i]), scGeno.begin() + i* genoSize());
			scInfoPtr[i] = m_inds[ scIndex[i]].infoPtr();
			copy( ind(scIndex[i]).infoBegin(), ind(scIndex[i]).infoEnd(), 
				scInfo.begin() + i* infoSize());
		}

		DBG_DO(DBG_POPULATION, cout << "Shallow copied" << scIndex << endl);

		/// sort the pointers!
		sort( scPtr.begin(), scPtr.end());
		sort( scInfoPtr.begin(), scInfoPtr.end());

		/// copy back.
		for(i=0, iEnd =scIndex.size(); i < iEnd;  i++)
		{
			m_inds[ scIndex[i] ].setGenoPtr( scPtr[i]);
			copy( scGeno.begin() + i*  genoSize(), scGeno.begin() + (i+1)*  genoSize(),
				indGenoBegin( scIndex[i] ));
			m_inds[ scIndex[i] ].setInfoPtr( scInfoPtr[i]);
			copy( scInfo.begin() + i*  infoSize(), scInfo.begin() + (i+1)*  infoSize(),
				ind(scIndex[i]).infoBegin());
			m_inds[ scIndex[i] ].setShallowCopied(false);
		}
		setShallowCopied(false);
		return;
	}

	population& LoadPopulation(const string& file,
		const string& format)
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
}
