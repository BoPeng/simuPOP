/**
 *  $File: qtrait.cpp $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
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

#include "qtrait.h"

namespace simuPOP {
bool quanTrait::apply(population & pop)
{
	UINT idx = pop.infoIdx(infoField(0));

	UINT ansGen = 0;

	if (m_ancGen == -1)
		ansGen = pop.ancestralGens();
	else if (m_ancGen > 0) {
		if (static_cast<UINT>(m_ancGen) > pop.ancestralGens())
			ansGen = pop.ancestralGens();
		else
			ansGen = m_ancGen;
	}

	for (UINT i = 0; i <= ansGen; ++i) {
		pop.useAncestralGen(i);

		// we need info to be in order
		IndInfoIterator traitIt = pop.infoBegin(idx);
		for (IndIterator it = pop.indIterator(); it.valid(); ++it)
			*traitIt++ = qtrait(& * it) ;
	}
	pop.useAncestralGen(0);

	return true;
}


double mapQuanTrait::qtrait(individual * ind)
{
	string key;

	for (vectorlu::iterator loc = m_loci.begin(); loc != m_loci.end(); ++loc) {
		// get genotype of ind
		Allele a = ToAllele(ind->allele(*loc, 0));
		Allele b = ToAllele(ind->allele(*loc, 1));

		if (loc != m_loci.begin() )
			key += '|';

		if (!m_phase && a > b)  // ab=ba
			key += toStr(static_cast<int>(b)) + "-" + toStr(static_cast<int>(a));
		else
			key += toStr(static_cast<int>(a)) + "-" + toStr(static_cast<int>(b));
	}
	strDict::iterator pos = m_dict.find(key);

	DBG_ASSERT(pos != m_dict.end(), ValueError,
		"No qtrait value for genotype " + key);

	return rng().randNormal(pos->second, m_sigma) ;
}


maQuanTrait::maQuanTrait(const uintList & loci, const vectorf & qtrait, const uintList & wildtype,
	const floatList & sigma, int ancGen, int stage, int begin, int end, int step,
	const intList & at, const repList & rep, const subPopList & subPops,
	const vectorstr & infoFields) :
	quanTrait(ancGen, stage, begin, end, step, at, rep, subPops, infoFields),
	m_loci(loci.elems()), m_qtrait(qtrait), m_sigma(sigma.elems()), m_wildtype(wildtype.elems())
{
	if (m_sigma.empty())
		m_sigma.resize(m_qtrait.size(), 0.);
	else if (m_sigma.size() == 1)
		m_sigma.resize(m_qtrait.size(), m_sigma[0]);

	DBG_ASSERT(m_qtrait.size() == static_cast<UINT>(pow(static_cast<double>(3),
														static_cast<double>(m_loci.size()))),
		ValueError, "Please specify qtrait for every combination of genotype.");
	DBG_ASSERT(m_sigma.size() == m_qtrait.size(), ValueError,
		"Size of sigma does not match that of qtrait");
};


double maQuanTrait::qtrait(individual * ind)
{
	UINT index = 0;

	for (vectorlu::iterator loc = m_loci.begin(); loc != m_loci.end(); ++loc) {
		// get genotype of ind
		Allele a = ToAllele(ind->allele(*loc, 0));
		Allele b = ToAllele(ind->allele(*loc, 1));

		int numWildtype = 0;

		// count number of wildtype
		if (find(m_wildtype.begin(), m_wildtype.end(), AlleleUnsigned(a)) != m_wildtype.end() )
			numWildtype++;

		if (find(m_wildtype.begin(), m_wildtype.end(), AlleleUnsigned(b)) != m_wildtype.end() )
			numWildtype++;
		index = index * 3 + 2 - numWildtype;
	}

	return rng().randNormal(m_qtrait[index], m_sigma[index]);
}


double mlQuanTrait::qtrait(individual * ind)
{
	if (m_mode == Multiplicative) {
		double fit = 1;
		for (vectorop::iterator s = m_qtraits.begin(), sEnd = m_qtraits.end();
		     s != sEnd; ++s)
			fit *= static_cast<quanTrait *>(*s)->qtrait(ind);
		return rng().randNormal(fit, m_sigma);
	} else if (m_mode == Additive) {
		double fit = 0;
		for (vectorop::iterator s = m_qtraits.begin(), sEnd = m_qtraits.end();
		     s != sEnd; ++s)
			fit += static_cast<quanTrait *>(*s)->qtrait(ind);
		return rng().randNormal(fit, m_sigma);
	}
	return 0.;
}


double pyQuanTrait::qtrait(individual * ind)
{
	if (m_len == 0) {
		m_len = m_loci.size() * ind->ploidy();
		m_alleles.resize(m_len);
		m_numArray = Allele_Vec_As_NumArray(m_alleles.begin(), m_alleles.end() );
	}

	DBG_FAILIF(static_cast<size_t>(m_len) != ind->ploidy() * m_loci.size(),
		SystemError,
		"Length of m_len is wrong. Have you changed pop type?");

	UINT pEnd = ind->ploidy();
	for (size_t i = 0, iEnd = m_loci.size(), j = 0; i < iEnd; ++i)
		for (UINT p = 0; p < pEnd; ++p)
			m_alleles[j++] = ToAllele(ind->allele(m_loci[i], p));

	return m_func(PyObj_As_Double, "(O)", m_numArray);
}


}
