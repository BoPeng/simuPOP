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

#include "individual.h"

namespace simuPOP
{

#ifdef SIMUMPI
	Allele individual::allele(UINT index) const
	{
		CHECKRANGEGENOSIZE(index);
		// find out which node has the allele and broadcast it.
		UINT p = index / totNumLoci();
		UINT locus = index - p * totNumLoci();
		std::pair<UINT, UINT> chIdx = chromLocusPair(locus);
		UINT rank = rankOfChrom(chIdx.first);
		Allele val = 0;
		if (mpiRank() == rank)
			val = *(m_genoPtr + chIdx.second +	  // infex in chrom
				localChromBegin(chIdx.first) +	  // local chrom begin
				p*localNumLoci());
		broadcast(mpiComm(), val, rank);
		return val;
	}

	Allele individual::allele(UINT index, UINT p) const
	{
		CHECKRANGEABSLOCUS(index);
		CHECKRANGEPLOIDY(p);

		std::pair<UINT, UINT> chIdx = chromLocusPair(index);
		UINT rank = rankOfChrom(chIdx.first);
		Allele val = 0;
		if (mpiRank() == rank)
			val = *(m_genoPtr
				+ index - beginLocus()			  // index is within a  block
				+ p*localNumLoci());
		broadcast(mpiComm(), val, rank);
		return val;
	}

	Allele individual::allele(UINT index, UINT p, UINT ch) const
	{
		CHECKRANGELOCUS(ch, index);
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(ch);

		UINT rank = rankOfChrom(ch);
		Allele val = 0;
		if(mpiRank() == rank)
			val = *(m_genoPtr + index + localChromBegin(ch)
				+ p*localNumLoci());
		broadcast(mpiComm(), val, rank);
		return val;
	}

	void individual::setAllele(Allele allele, UINT index)
	{
		CHECKRANGEGENOSIZE(index);

		// find out which node has the allele and broadcast it.
		UINT p = index / totNumLoci();
		UINT locus = index - p * totNumLoci();
		std::pair<UINT, UINT> chIdx = chromLocusPair(locus);
		UINT rank = rankOfChrom(chIdx.first);
		if (mpiRank() == rank)
			*(m_genoPtr							  // infex in chrom
			+ localChromBegin(chIdx.first)
			+ chIdx.second
			+ p*localNumLoci()) = allele;
	}

	/// set allele from an index.
	/** \param allele allele to set
	\param index index from the begining of genotype
	\param p on p'th set of chromosome, p=0 by default
	 */
	void individual::setAllele(Allele allele, UINT index, UINT p)
	{
		CHECKRANGEABSLOCUS(index);
		CHECKRANGEPLOIDY(p);

		std::pair<UINT, UINT> chIdx = chromLocusPair(index);
		UINT rank = rankOfChrom(chIdx.first);
		if (mpiRank() == rank)
			*(m_genoPtr + chIdx.second + localChromBegin(chIdx.first)
			+ p*localNumLoci()) = allele;
	}

	void individual::setAllele(Allele allele, UINT index, UINT p, UINT ch)
	{
		CHECKRANGELOCUS(ch, index);
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(ch);

		UINT rank = rankOfChrom(ch);
		if(mpiRank() == rank)
			*(m_genoPtr + index + localChromBegin(ch) +
			p*localNumLoci()) = allele;
	}

	Sex individual::sex() const
	{
		bool isFemale = ISSETFLAG(m_flags, m_flagFemale) ;
		broadcast(mpiComm(), isFemale, 0);
		if(isFemale)
			return Female;
		else
			return Male;
	};

	void individual::setSex(Sex sex)
	{
		if(mpiRank() == 0)
		{
			CHECKRANGESEX(sex);
			if( sex == Male )
				RESETFLAG(m_flags, m_flagFemale);
			else
				SETFLAG(m_flags, m_flagFemale);
		}
	}

	bool individual::affected() const
	{
		bool isAffected = ISSETFLAG(m_flags, m_flagAffected);
		broadcast(mpiComm(), isAffected, 0);
		return isAffected;
	}

	void individual::setAffected(bool affected)
	{
		if(mpiRank() == 0)
		{
			if(affected)
				SETFLAG(m_flags, m_flagAffected);
			else
				RESETFLAG(m_flags, m_flagAffected);
		}
	}

	/// get info
	InfoType individual::info(UINT idx) const
	{
		InfoType info = 0;
		CHECKRANGEINFO(idx);
		if (mpiRank() == 0)
		{
			// broad cast the value to all nodes
			info = m_infoPtr[idx];
		}
		broadcast(mpiComm(), info, 0);
		return info;
	}

	/// set info
	void individual::setInfo(InfoType value, UINT idx)
	{
		if (mpiRank()==0)
		{
			CHECKRANGELOCALINFO(idx);
			m_infoPtr[idx] = value;
		}
	}

	/// get info
	InfoType individual::info(const string& name) const
	{
		InfoType info = 0;
		// use glocal infoIdx (with broadcast to trigger
		// IndexError properly
		int idx = infoIdx(name);
		if (mpiRank()==0)
			info = m_infoPtr[idx];
		broadcast(mpiComm(), info, 0);
		return info;
	}

	/// set info
	void individual::setInfo(InfoType value, const string& name)
	{
		if (mpiRank()==0)
		{
			int idx = localInfoIdx(name);
			DBG_ASSERT(idx>=0, IndexError,
				"Info name " + name + " is not a valid info field name");
			m_infoPtr[idx] = value;
		}
	}
#endif

	PyObject* individual::arrGenotype()
	{
		// this &* is to avoid any possible type mismatch thing.
		// there is some magic here, for MPI module
		// PyObject are different from node to node, but the
		// referred value is the same. The right one is used and
		// the value is boardcasted.
#ifdef SIMUMPI
		// which portion is this piece of array in?
		return Allele_Vec_As_NumArray(m_genoPtr, m_genoPtr + localGenoSize(),
			genoSize(), totNumLoci(), beginLocus(), endLocus());
#else
		return Allele_Vec_As_NumArray(m_genoPtr, m_genoPtr + genoSize());
#endif
	}

	/// return genotype as python Numeric.array object
	/// This is the p'th copy of chromosomes
	PyObject* individual::arrGenotype(UINT p)
	{
		CHECKRANGEPLOIDY(p);
#ifdef SIMUMPI
		return Allele_Vec_As_NumArray( m_genoPtr + p*localNumLoci(),
			m_genoPtr + (p+1)*localNumLoci(),
			totNumLoci(), totNumLoci(), beginLocus(), endLocus());
#else
		return Allele_Vec_As_NumArray( m_genoPtr + p*totNumLoci(),
			m_genoPtr + (p+1)*totNumLoci() );
#endif
	}

	/// return genotype as python Numeric.array object
	/// This is the ch chromosome of the pth copy of chromosome
	PyObject* individual::arrGenotype(UINT p, UINT ch)
	{
		CHECKRANGEPLOIDY(p);
#ifdef SIMUMPI
		return Allele_Vec_As_NumArray( m_genoPtr + p*localNumLoci() + localChromBegin(ch),
			m_genoPtr + p*localNumLoci() + localChromEnd(ch),
			numLoci(ch), totNumLoci(), 0, localNumLoci());
#else
		return Allele_Vec_As_NumArray( m_genoPtr + p*totNumLoci() + chromBegin(ch),
			m_genoPtr + p*totNumLoci() +chromEnd(ch));
#endif
	}

	PyObject* individual::arrInfo()
	{
#ifdef SIMUMPI
		if (mpiRank() == 0)
			return Info_Vec_As_NumArray(m_infoPtr, m_infoPtr + localInfoSize() );
		else
			throw ValueError("arrInfo() is not available for non-head nodes");
#else
		return Info_Vec_As_NumArray(m_infoPtr, m_infoPtr + infoSize() );
#endif
	}

	individual& individual::operator= (const individual& rhs)
	{
		setShallowCopied(true);

		/// when the source is being moved. Its relative position may change
		/// so it also becomes shallowCopied.
		const_cast<individual&>(rhs).setShallowCopied(true);

		m_flags = rhs.m_flags;
		setSubPopID(rhs.subPopID());
		setGenoPtr(rhs.genoPtr());
		setInfoPtr(rhs.infoPtr());
		// also copy genoStru pointer...
		this->setGenoStruIdx(rhs.genoStruIdx());
		return *this;
	}

	individual& individual::copyFrom( const individual& rhs)
	{
		m_flags = rhs.m_flags;
		setSubPopID(rhs.subPopID());
		copy(rhs.genoBegin(), rhs.genoEnd(), genoBegin());
		copy(rhs.infoBegin(), rhs.infoEnd(), infoBegin());
		// also copy genoStru pointer...
		this->setGenoStruIdx(rhs.genoStruIdx());
		setShallowCopied(false);
		return *this;
	}

	bool individual::operator== (const individual& rhs) const
	{
#ifdef SIMUMPI
		bool equal = true;
		if( genoStruIdx() != rhs.genoStruIdx() )
		{
			DBG_DO(DBG_POPULATION, cout << "Geno stru different" << endl);
			equal = false;
		}

		if(mpiRank() == 0)
		{
			if(ISSETFLAG(m_flags, m_flagFemale) != ISSETFLAG(rhs.m_flags, m_flagFemale)
				|| ISSETFLAG(m_flags, m_flagAffected) != ISSETFLAG(rhs.m_flags, m_flagAffected) )
			{
				DBG_DO(DBG_POPULATION, cout << "Flags different: sex "
					<< ISSETFLAG(m_flags, m_flagFemale) << " vs " << ISSETFLAG(rhs.m_flags, m_flagFemale) << ", aff "
					<< ISSETFLAG(m_flags, m_flagAffected) << " vs " << ISSETFLAG(rhs.m_flags, m_flagAffected)
					<< endl);
				equal = false;
			}

			if(equal)
			{
				for( UINT i=0, iEnd = localInfoSize(); i < iEnd;  ++i)
				{
					if( info(i) != rhs.info(i) )
					{
						DBG_DO(DBG_POPULATION, cout << "Info different" << endl);
						equal = false;
						break;
					}
				}
			}
		}

		for( UINT i=0, iEnd = localGenoSize(); i < iEnd;  ++i)
		{
			// no shift
			if( *(m_genoPtr+i) != *(rhs.m_genoPtr+i) )
			{
				DBG_DO(DBG_POPULATION, cout << "Genotype different" << endl);
				equal = false;
				break;
			}
		}
		bool res;
		reduce(mpiComm(), equal, res, std::logical_and<bool>(), 0);
		broadcast(mpiComm(), res, 0);
		return res;
#else
		if( genoStruIdx() != rhs.genoStruIdx() )
			return false;

		if(ISSETFLAG(m_flags, m_flagFemale) != ISSETFLAG(rhs.m_flags, m_flagFemale)
			|| ISSETFLAG(m_flags, m_flagAffected) != ISSETFLAG(rhs.m_flags, m_flagAffected) )
			return false;

		for( UINT i=0, iEnd = infoSize(); i < iEnd;  ++i)
			if( info(i) != rhs.info(i) )
				return false;

		for( UINT i=0, iEnd = genoSize(); i < iEnd;  ++i)
			if( *(m_genoPtr+i) != *(rhs.m_genoPtr+i) )
				return false;
#endif
		return true;
	}

	int individual::__cmp__(const individual& rhs) const
	{
#ifdef SIMUMPI
		bool equal = true;
		if( genoStruIdx() != rhs.genoStruIdx() )
		{
			DBG_DO(DBG_POPULATION, cout << "Geno stru different" << endl);
			equal = false;
		}

		if(mpiRank() == 0)
		{
			if( equal && m_flags != rhs.m_flags )
			{
				DBG_DO(DBG_POPULATION, cout << "Flags different" << endl);
				equal = false;
			}

			if(equal)
			{
				for( UINT i=0, iEnd = localInfoSize(); i < iEnd;  ++i)
				{
					if( info(i) != rhs.info(i) )
					{
						DBG_DO(DBG_POPULATION, cout << "Info different" << endl);
						equal = false;
						break;
					}
				}
			}
		}

		if(equal)
		{
			for( UINT i=0, iEnd = localGenoSize(); i < iEnd;  ++i)
				if( allele(i) != rhs.allele(i) )
			{
				DBG_DO(DBG_POPULATION, cout << "Genotype different" << endl);
				equal = false;
				break;
			}
		}
		bool res;
		reduce(mpiComm(), equal, res, std::logical_and<bool>(), 0);
		broadcast(mpiComm(), res, 0);
		return res;

#else

		if( genoStruIdx() != rhs.genoStruIdx() )
			return 1;

		if( m_flags != rhs.m_flags )
			return 1;

		for( UINT i=0, iEnd = infoSize(); i < iEnd;  ++i)
			if( info(i) != rhs.info(i) )
				return 1;

		for( UINT i=0, iEnd = genoSize(); i < iEnd;  ++i)
			if( allele(i) != rhs.allele(i) )
				return 1;

		return 0;
#endif
	}

	void individual::swap(individual& ind, bool swapContent)
	{
		if( genoStruIdx() != ind.genoStruIdx() )
			throw SystemError("Can only swap individuals with different geno structure.");

		std::swap(m_subPopID, ind.m_subPopID);
		std::swap(m_infoPtr, ind.m_infoPtr);

		if(swapContent)
		{
			Allele tmp;
			for(UINT i=0, iEnd = genoSize(); i < iEnd;  i++)
			{
				tmp = m_genoPtr[i];
				m_genoPtr[i] = ind.m_genoPtr[i];
				ind.m_genoPtr[i] = tmp;
			}
		}
		else
		{
			setShallowCopied(true);
			ind.setShallowCopied(true);
			std::swap(m_genoPtr, ind.m_genoPtr);
		}
	}

	void individual::display( ostream& out, int width, const vectori& chrom, const vectori& loci)
	{
		out << sexChar() << affectedChar() << " ";
		DBG_DO(DBG_POPULATION, out <<  subPopID() << " ");
		for(UINT p=0, pEnd = ploidy(); p < pEnd;  ++p)
		{
			//      copy( genoBegin()+i, genoBegin()+i+totNumLoci(),
			//        std::ostream_iterator<string>(out, outputSeparator()) );
			if(chrom.empty() && loci.empty())
			{
				for(UINT ch=0, chEnd=numChrom(); ch<chEnd; ++ch)
				{
					for(UINT j = 0, jEnd = numLoci(ch); j < jEnd;  ++j)
						out << setw(width) << alleleChar(j, p, ch);
					out << " ";
				}
			}
			else if(! chrom.empty() && loci.empty())
			{
				for(vectori::const_iterator ch=chrom.begin(); ch != chrom.end(); ++ch)
				{
					for(UINT j = 0, jEnd = numLoci(*ch); j < jEnd;  ++j)
						out << setw(width) << alleleChar(j, p, *ch);
					out << " ";
				}
			}
			else if( chrom.empty() && ! loci.empty())
			{
				for(vectori::const_iterator loc=loci.begin(); loc != loci.end(); ++loc)
					out << setw(width) << alleleChar(*loc, p);
				out << " ";
			}
			else								  // both specified
				throw ValueError("Please specify only one of chrom and loci.");

			if( p != pEnd-1)
				out << "| ";
		}
	}
}
