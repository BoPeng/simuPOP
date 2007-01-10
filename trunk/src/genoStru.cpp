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

#include "genoStru.h"

namespace simuPOP
{
	GenoStructure::GenoStructure(UINT ploidy, const vectoru& loci, bool sexChrom,
		const vectorf& lociPos, const vectorstr& alleleNames,
		const vectorstr& lociNames, UINT maxAllele, const vectorstr& infoFields,
		const vectori& chromMap)
		:m_ploidy(ploidy),  m_numChrom(loci.size()), m_numLoci(loci), m_sexChrom(sexChrom),
		m_lociPos(lociPos), m_chromIndex(loci.size()+1),
		m_alleleNames(alleleNames), m_lociNames(lociNames),
		m_maxAllele(maxAllele),
		m_infoFields(infoFields),
		m_chromMap(chromMap)
	#ifdef SIMUMPI
		, m_beginChrom(0), m_endChrom(0), m_beginLocus(0), m_endLocus(0),
		m_localNumLoci(0), m_localGenoSize(0)
	#endif
	{
#ifdef BINARYALLELE
		DBG_ASSERT(maxAllele == 1, ValueError,
			"max allele must be 1 for binary modules");
#endif
		DBG_ASSERT( ploidy >= 1, ValueError,
			"Ploidy must be >= 1. Given " + toStr(ploidy) );

		// default: one chromosome, one locus
		// otherwise, Loci copies from loci
		if (loci.empty())
		{
			m_numChrom = 1;
			m_numLoci.resize(1);
			m_numLoci[0] = 1;
			m_chromIndex.resize(2);
		}

		// build chromosome index
		ULONG i, j;
		for (m_chromIndex[0] = 0, i = 1; i <= m_numChrom; ++i)
			m_chromIndex[i] = m_chromIndex[i - 1] + m_numLoci[i - 1];

		m_totNumLoci = m_chromIndex[m_numChrom];
		m_genoSize = m_totNumLoci*m_ploidy;

		// if lociPos not specified, use 1,2,3.. 1,2,3. .. on each chromosome.
		if( m_lociPos.empty() )
		{
			m_lociPos.resize(m_totNumLoci);
			for (i = 0; i < m_numChrom; ++i)
				for (j = 0; j < m_numLoci[i]; j++)
					m_lociPos[m_chromIndex[i]+j] = j + 1;
		}
#ifndef OPTIMIZED
		else									  // check loci distance
		{
			// loci distance, if specified, chould have length of chromosome.
			DBG_FAILIF( m_lociPos.size() != m_totNumLoci, ValueError,
				"You should specify loci distance for every locus (" + toStr(m_totNumLoci) + ")");

			for( i=0; i < m_numChrom; ++i)
				for (j = 0; j < m_numLoci[i]; ++j)
					DBG_FAILIF( j > 0 && fcmp_lt( m_lociPos[m_chromIndex[i]+j], m_lociPos[m_chromIndex[i]+j-1]),
						ValueError, "Loci distance must be in order.");
		}
#endif

		if( m_lociNames.empty())
		{
			m_lociNames.resize(m_totNumLoci);
			for(i = 0; i < m_numChrom; ++i)
				for (j = 0; j < m_numLoci[i]; j++)
					m_lociNames[m_chromIndex[i]+j] = "loc"+ toStr(i+1) + "-" + toStr(j+1);
		}
		DBG_ASSERT( m_lociNames.size() == m_totNumLoci, ValueError,
			"Loci names, if specified, should be given to every loci");
		DBG_WARNING( (!m_alleleNames.empty()) && m_alleleNames.size() != m_maxAllele+1,
			"Not all allele names are given. ");

#ifdef SIMUMPI
		// no information fields for non-head nodes
		if (mpiRank() != 0)
			m_infoFields.clear();

		if( m_chromMap.empty())
			m_chromMap = vectori(m_numChrom, 1);
		// begining and end chromosome?
		UINT rank = mpiRank();
		if(rank == 0)
		{
			m_beginChrom = 0;
			m_endChrom = 0;
		}
		else if(rank == 1)
		{
			m_beginChrom = 0;
			m_endChrom = m_chromMap[0];
		}
		else if(rank > 1)
		{
			size_t sum = 0;
			for(i=0; i < rank-1; ++i)
				sum += m_chromMap[i];
			m_beginChrom = sum;
			m_endChrom = sum + m_chromMap[i];
		}
		m_beginLocus = m_chromIndex[m_beginChrom];
		m_endLocus = m_chromIndex[m_endChrom];
		m_localNumLoci = m_endLocus - m_beginLocus;
		m_localGenoSize = m_localNumLoci*m_ploidy;
		// local chromosome index
		m_localChromIndex.resize(m_endChrom-m_beginChrom+1);
		for (m_localChromIndex[0] = 0, i = 1; i <= m_endChrom - m_beginChrom; ++i)
			m_localChromIndex[i] = m_localChromIndex[i - 1] + m_numLoci[i - 1 + m_beginChrom];
		DBG_DO(DBG_POPULATION, cout << "rank " << rank
			<< " begin Locus " << m_beginLocus
			<< " end Locus " << m_endLocus
			<< " begin Chrom " << m_beginChrom
			<< " ebd Chrom " << m_endChrom
			<< " local numLoci " << m_localNumLoci
			<< " loca genosize " << m_localGenoSize << endl);
#endif
	}

	bool GenoStructure::operator== (const GenoStructure& rhs)
	{
		// compare pointer directly will be fastest
		if(this == &rhs || (
			( m_ploidy == rhs.m_ploidy) &&
			( m_numLoci == rhs.m_numLoci) &&
			( m_sexChrom == rhs.m_sexChrom) &&
			( m_lociPos == rhs.m_lociPos) &&
			( m_alleleNames == rhs.m_alleleNames) &&
			( m_lociNames == rhs.m_lociNames) &&
			( m_maxAllele == rhs.m_maxAllele) &&
			( m_infoFields == rhs.m_infoFields)
			// chromosome map etc are not compared since
			// they can be different from run to run.
			))
			return true;
		else
			return false;
	}

	// initialize static variable s)genoStruRepository.
	vector<GenoStructure> GenoStruTrait::s_genoStruRepository = vector<GenoStructure>();

	void GenoStruTrait::setGenoStructure(UINT ploidy, const vectoru& loci, bool sexChrom,
		const vectorf& lociPos, const vectorstr& alleleNames,
		const vectorstr& lociNames, UINT maxAllele, const vectorstr& infoFields,
		const vectori& chromMap)
	{
		/// only allow for TraitMaxIndex-1 different genotype structures
		/// As a matter of fact, most simuPOP scripts have only one
		/// population type.
		if( s_genoStruRepository.size() == TraitMaxIndex-1 )
		{
			throw SystemError("This simuPOP library only allows " + toStr(TraitMaxIndex-1)
				+ " different genotype structures. \n" +
				+ "If you do need more structures, modify individual.h/TraitMaxType and " +
				+ "recompile simuPOP.");
		}

		GenoStructure tmp = GenoStructure( ploidy, loci, sexChrom,
			lociPos, alleleNames, lociNames, maxAllele, infoFields, chromMap);

		for(TraitIndexType it = 0; it < s_genoStruRepository.size();
			++it)
		{
			// object comparison
			if( s_genoStruRepository[it] == tmp )
			{
				m_genoStruIdx = it;
				return;
			}
		}
		// if not found
		s_genoStruRepository.push_back(tmp);
		// the last one.
		m_genoStruIdx = s_genoStruRepository.size()-1;
	}

	void GenoStruTrait::setGenoStructure(GenoStructure& rhs)
	{
		for(TraitIndexType it = 0; it < s_genoStruRepository.size();
			++it)
		{
												  // object comparison
			if( s_genoStruRepository[it] == rhs )
			{
				m_genoStruIdx = it;
				return;
			}
		}

		// if not found, make a copy and store it.
		s_genoStruRepository.push_back( rhs );
		m_genoStruIdx = s_genoStruRepository.size() - 1;
	}

	string GenoStruTrait::ploidyName() const
	{
		DBG_FAILIF( m_genoStruIdx == TraitMaxIndex, SystemError,
			"PloidyName: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 1)
			return "haploid";
		else if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 2)
			return "diploid";
		else if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 3)
			return "triploid";
		else if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 4)
			return "tetraploid";
		else
			return toStr(s_genoStruRepository[m_genoStruIdx].m_ploidy) + "-polid";
	}

	std::pair<UINT, UINT> GenoStruTrait::chromLocusPair(UINT locus) const
	{
		CHECKRANGEABSLOCUS(locus);

		pair<UINT, UINT> loc;

		for(UINT i=1, iEnd =numChrom(); i <= iEnd;  ++i)
		{
			if( s_genoStruRepository[m_genoStruIdx].m_chromIndex[i] > locus)
			{
				loc.first = i-1;
				loc.second = locus - s_genoStruRepository[m_genoStruIdx].m_chromIndex[i-1];
				break;
			}
		}
		return loc;
	}

	string GenoStruTrait::alleleName(const Allele allele) const
	{
#ifndef BINARYALLELE
		DBG_FAILIF(allele > s_genoStruRepository[m_genoStruIdx].m_maxAllele,
			IndexError, "Allele out of range of 0 ~ " +
			toStr(s_genoStruRepository[m_genoStruIdx].m_maxAllele));
		if( allele < s_genoStruRepository[m_genoStruIdx].m_alleleNames.size() )
		{
			DBG_FAILIF( allele >= s_genoStruRepository[m_genoStruIdx].m_alleleNames.size() ,
				IndexError, "No name for allele " + toStr(static_cast<UINT>(allele)));

			return s_genoStruRepository[m_genoStruIdx].m_alleleNames[allele];
		}
		else
			return toStr(static_cast<int>(allele));
#else
		if( static_cast<unsigned>(allele) < s_genoStruRepository[m_genoStruIdx].m_alleleNames.size() )
			return s_genoStruRepository[m_genoStruIdx].m_alleleNames[allele];
		else if(allele)
			return "1";
		else
			return "0";
#endif

	}

	UINT GenoStruTrait::infoIdx(const string& name) const
	{
		vectorstr& names = s_genoStruRepository[m_genoStruIdx].m_infoFields;

		for(UINT i=0; i< names.size(); ++i)
		{
			if(names[i] == name)
				return i;
		}
		throw IndexError("Info field '" + name + "' is not found. "
			"Plese use infoFields=['" + name + "'] option of population() during construction\n"
			"or use addInfoField('" + name + "') to add to an existing population.");
		// this should never be reached.
		return 0;
	}

#ifdef SIMUMPI
	/// return node rank by chromosome number, according to map on setChromMap
	UINT GenoStruTrait::rankOfChrom(UINT chrom) const
	{
		vectori & map = s_genoStruRepository[m_genoStruIdx].m_chromMap;

		for(size_t i=0, sum = 0; i<map.size(); ++i)
		{
			sum += map[i];
			if(chrom < sum)
				return i+1;
		}
		DBG_FAILIF(true, IndexError, "Chromosome " + toStr(chrom) + " is not on chromosome map");
	}

	/// begin chromosome for a given rank
	UINT GenoStruTrait::beginChromOfRank(UINT rank) const
	{
		if (rank <= 1)
			return 0;

		vectori & map = s_genoStruRepository[m_genoStruIdx].m_chromMap;

		DBG_ASSERT(rank <= map.size(), IndexError, "Given rank " + toStr(rank) + " is invalid.");

		size_t sum = 0;
		for(size_t i=0; i<rank-1; ++i)
			sum += map[i];
		return sum;
	}

	/// end chromosome for a given rank (actually begin chromosome for the next rank)
	UINT GenoStruTrait::endChromOfRank(UINT rank) const
	{
		vectori & map = s_genoStruRepository[m_genoStruIdx].m_chromMap;

		DBG_ASSERT(rank <= map.size(), IndexError, "Given rank " + toStr(rank) + " is invalid.");
		size_t sum = 0;
		for(size_t i=0; i < rank; ++i)
			sum += map[i];
		return sum;
	}
#endif

}
