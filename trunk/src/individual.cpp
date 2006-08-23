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
	GenoStructure::GenoStructure(UINT ploidy, const vectoru& loci, bool sexChrom,
		const vectorf& lociPos, const vectorstr& alleleNames,
		const vectorstr& lociNames, UINT maxAllele, const vectorstr& infoName)
		:m_ploidy(ploidy),  m_numChrom(loci.size()), m_numLoci(loci), m_sexChrom(sexChrom),
		m_lociPos(lociPos), m_chromIndex(loci.size()+1),
		m_alleleNames(alleleNames), m_lociNames(lociNames), 
		m_maxAllele(maxAllele), m_infoName(infoName)
	{
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
			"Not all allele names are given. ")
	}

	// initialize static variable s)genoStruRepository.
	vector<GenoStructure> GenoStruTrait::s_genoStruRepository = vector<GenoStructure>();

	void GenoStruTrait::setGenoStructure(UINT ploidy, const vectoru& loci, bool sexChrom,
		const vectorf& lociPos, const vectorstr& alleleNames,
		const vectorstr& lociNames, UINT maxAllele, const vectorstr& infoName)
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
			lociPos, alleleNames, lociNames, maxAllele, infoName);

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

}
