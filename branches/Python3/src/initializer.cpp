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

#include "initializer.h"

namespace simuPOP {

void initializer::setRanges(population & pop)
{
	m_ranges = m_indRange;

	if (!m_ranges.empty())
		return;

	if (m_subPop.empty()) {
		m_ranges.resize(pop.numSubPop());
		for (size_t i = 0; i < pop.numSubPop(); ++i) {
			m_ranges[i].resize(2, pop.subPopBegin(i));
			m_ranges[i][1] = pop.subPopEnd(i);
		}
	} else {
		m_ranges.resize(m_subPop.size());
		for (size_t i = 0; i < m_subPop.size(); ++i) {
			m_ranges[i].resize(2, pop.subPopBegin(m_subPop[i]));
			m_ranges[i][1] = pop.subPopEnd(m_subPop[i]);
		}
	}
}


bool initSex::apply(population & pop)
{
	setRanges(pop);

	for (size_t rg = 0; rg < m_ranges.size(); ++rg) {
		IndIterator ind = pop.indBegin() + m_ranges[rg][0];
		IndIterator right = pop.indBegin() + m_ranges[rg][1];

		if (m_sex.empty())
			for (; ind != right; ++ind)
				ind->setSex(rng().randUniform01() < m_maleFreq ? Male : Female);
		else
			for (size_t idx = 0, sexSz = m_sex.size(); ind != right; ++ind, ++idx)
				ind->setSex(m_sex[idx % sexSz] == 1 ? Male : Female);
	}
	return true;
}


bool initByFreq::apply(population & pop)
{
	// initialize m_ranges
	setRanges(pop);

	DBG_FAILIF(m_alleleFreq.size() > 1 && m_alleleFreq.size() != m_ranges.size(),
		ValueError, "Ranges and values should have the same length");

	for (size_t rg = 0; rg < m_ranges.size(); ++rg) {
		vectorf & alleleFreq = m_alleleFreq.size() == 1 ? m_alleleFreq[0] : m_alleleFreq[rg];

		ULONG left = m_ranges[rg][0];
		ULONG right = m_ranges[rg][1];

		DBG_FAILIF(left > pop.popSize() || right > pop.popSize() || left > right,
			ValueError, "Invaid range boundary: " + toStr(left) + " - " + toStr(right - 1));

		// Weightedsampler ws(rng(), incFreq);
		Weightedsampler ws(rng(), alleleFreq);

		DBG_ASSERT(fcmp_eq(std::accumulate(alleleFreq.begin(), alleleFreq.end(), 0.), 1),
			SystemError, "Allele frequecies should add up to one.");

		vectoru atLoci = m_atLoci;
		if (m_atLoci.empty())
			for (size_t i = 0 ; i < pop.totNumLoci(); ++i)
				atLoci.push_back(i);

		pop.sortIndividuals();
		if (m_identicalInds) {
			for (vectoru::iterator locus = atLoci.begin(); locus != atLoci.end(); ++locus) {
				if (m_atPloidy == -1) {
					ws.get(pop.alleleBegin(*locus) + left * pop.ploidy(),
						pop.alleleBegin(*locus) + (left + 1) * pop.ploidy() );

					for (ULONG ind = left + 1; ind < right; ++ind)
						copy(pop.alleleBegin(*locus) + left * pop.ploidy(),
							pop.alleleBegin(*locus) + (left + 1) * pop.ploidy(),
							pop.alleleBegin(*locus) + ind * pop.ploidy());
				} else {                                          // only one of the copies (do it one by one?)
					UINT a = ws.get();
					for (ULONG ind = left; ind != right; ++ind)
						pop.ind(ind).setAllele(a, *locus, m_atPloidy);
				}
			}
		} else {                                                                        // not idential individuals
			if (m_atPloidy == -1) {
				for (vectoru::iterator locus = atLoci.begin(); locus != atLoci.end(); ++locus)
					ws.get(pop.alleleBegin(*locus) + left * pop.ploidy(),
						pop.alleleBegin(*locus) + right * pop.ploidy() );
			} else {                                                  // for only one ploidy
				for (vectoru::iterator locus = atLoci.begin(); locus != atLoci.end(); ++locus) {
					for (ULONG ind = left; ind != right; ++ind)
						// ws.get() return ULONG from 0 to len(allele frequency).
						pop.ind(ind).setAllele(static_cast<Allele>(ws.get()), *locus, m_atPloidy);
				}
			}
		}
	}                                                                                               // range
	return initSex::apply(pop);
}


bool initByValue::apply(population & pop)
{
	for (size_t i = 0; i < m_value.size(); ++i) {
		// extend m_value
		if (!m_proportion.empty() || m_atPloidy != -1)
			continue;

		if (m_atLoci.empty() && m_value[i].size() == pop.totNumLoci()) {
			m_value[i].resize(pop.totNumLoci() * pop.ploidy());
			for (UINT p = 1; p < pop.ploidy(); ++p)
				copy(m_value[i].begin(), m_value[i].begin() + pop.totNumLoci(),
					m_value[i].begin() + pop.totNumLoci() * p);
		}
		if (!m_atLoci.empty() && (m_value[i].size() == m_atLoci.size() ||
		                          m_value[i].size() == m_atLoci.size() * pop.ploidy())) {
			m_value[i].resize(m_atLoci.size() * pop.ploidy());
			for (UINT p = 1; p < pop.ploidy(); ++p)
				copy(m_value[i].begin(), m_value[i].begin() + m_atLoci.size(),
					m_value[i].begin() + m_atLoci.size() * p);
		}                                                                                           // no proportion
	}                                                                                               // check m_value[i]

	DBG_DO(DBG_INITIALIZER, cout << "Size of src is " << m_value[0].size() << endl);

#ifndef OPTIMIZED
	UINT gSz = m_value[0].size();
	for (size_t v = 1; v < m_value.size(); ++v)
		if (m_value[v].size() != gSz)
			throw ValueError("Given values should have the same length (either one copy of chromosomes or the whole genotype.");
#endif

	setRanges(pop);

	vectoru atLoci = m_atLoci;
	if (m_atLoci.empty())
		for (size_t i = 0 ; i < pop.totNumLoci(); ++i)
			atLoci.push_back(i);

	DBG_FAILIF(m_proportion.empty() && m_value.size() > 1
		&& m_value.size() != m_ranges.size(),
		ValueError, "Ranges and values should have the same length");

	// for each range
	for (size_t rg = 0; rg < m_ranges.size(); ++rg) {
		// we have left and right
		ULONG left = m_ranges[rg][0];
		ULONG right = m_ranges[rg][1];

		DBG_FAILIF(left > pop.popSize() || right > pop.popSize() || left > right,
			ValueError, "Invaid m_ranges boundary: " + toStr(left) + " - " + toStr(right));

		size_t lociSz = atLoci.size();
		size_t totNumLoci = pop.totNumLoci();
		if (m_proportion.empty()) {
			vectori & src = m_value.size() > 1 ? m_value[rg] : m_value[0];
			size_t srcSz = src.size();

			for (ULONG ind = left; ind < right; ++ind) {
				if (m_atPloidy == -1) {                                     // all copied of chromosome
					DBG_ASSERT(src.size() == atLoci.size() ||
						(src.size() == atLoci.size() * pop.ploidy() ),
						ValueError, "Length of value does not match atLoci size");
					for (size_t loc = 0; loc != srcSz ; ++loc)
						*(pop.indGenoBegin(ind) + atLoci[loc % lociSz] + loc / lociSz * totNumLoci) = src[loc];
				} else {                                          // one of the copies.
					DBG_ASSERT(src.size() == atLoci.size(), ValueError,
						"Ploidy is specified but the length of alleles do not match length of given allele array.");
					for (size_t loc = 0; loc != srcSz ; ++loc)
						*(pop.ind(ind).genoBegin(m_atPloidy) +
						  atLoci[loc % lociSz] + loc / lociSz * totNumLoci) = src[loc];
				}
			}
		} else {                                                                          // use proportion.
			size_t srcSz = m_value[0].size();
			Weightedsampler ws(rng(), m_proportion);

			for (ULONG ind = left; ind < right; ++ind) {
				if (srcSz == lociSz) {                                  // one by one
					if (m_atPloidy == -1) {
						for (UINT p = 0; p < pop.ploidy(); ++p) {
							UINT idx = ws.get();
							for (size_t loc = 0; loc != srcSz ; ++loc)
								*(pop.indGenoBegin(ind) + p * totNumLoci + atLoci[loc]) = m_value[idx][loc];
						}
					} else {
						UINT idx = ws.get();
						for (size_t loc = 0; loc != srcSz ; ++loc)
							*(pop.indGenoBegin(ind) + m_atPloidy * totNumLoci +
							  atLoci[loc]) = m_value[idx][loc];
					}
				} else {                                          // who geno (at loci .. though)
					if (m_atPloidy == -1) {
						UINT idx = ws.get();
						for (size_t loc = 0; loc != srcSz ; ++loc)
							*(pop.indGenoBegin(ind) + atLoci[loc % lociSz] +
							  loc / lociSz * totNumLoci) = m_value[idx][loc];
					} else {
						UINT idx = ws.get();
						for (size_t loc = 0; loc != srcSz ; ++loc)
							*(pop.ind(ind).genoBegin(m_atPloidy)
							  + atLoci[loc % lociSz] +
							  loc / lociSz * totNumLoci) = m_value[idx][loc];
					}
				}
			}
		}
	}
	return initSex::apply(pop);
}


bool pyInit::apply(population & pop)
{
	setRanges(pop);

	for (size_t rg = 0; rg < m_ranges.size(); ++rg) {
		ULONG it = m_ranges[rg][0];
		ULONG right = m_ranges[rg][1];
		for (; it != right; ++it) {
			UINT sp = pop.subPopIndPair(it).first;
			for (UINT al = 0, alEnd = pop.totNumLoci(); al < alEnd; ++al) {
				for (UINT p = 0, pEnd = pop.ploidy(); p < pEnd; ++p) {
					int resInt;
					PyCallFunc3(m_func, "(iii)", al, p, sp, resInt, PyObj_As_Int);
					pop.ind(it).setAllele(static_cast<Allele>(resInt), al, p);
				}
			}
		}
	}
	return initSex::apply(pop);
}


}
