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
		if( genoStruIdx() != rhs.genoStruIdx() )
			return false;

		if(ISSETFLAG(m_flags, m_flagFemale) != ISSETFLAG(rhs.m_flags, m_flagFemale)
			|| ISSETFLAG(m_flags, m_flagAffected) != ISSETFLAG(rhs.m_flags, m_flagAffected) )
			return false;

		for( UINT i=0, iEnd = infoSize(); i < iEnd;  ++i)
			if( info(i) != rhs.info(i) )
				return false;

		for( UINT i=0, iEnd = genoSize(); i < iEnd;  ++i)
			if( allele(i) != rhs.allele(i) )
				return false;

		return true;
	}

	int individual::__cmp__(const individual& rhs) const
	{
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
