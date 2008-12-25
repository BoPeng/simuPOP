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

#include "qtrait.h"

namespace simuPOP {
bool quanTrait::apply(population & pop)
{
	UINT idx = pop.infoIdx(infoField(0));

	UINT ansGen = 0;

	if (m_ancestralGen == -1)
		ansGen = pop.ancestralGens();
	else if (m_ancestralGen > 0) {
		if (static_cast<UINT>(m_ancestralGen) > pop.ancestralGens())
			ansGen = pop.ancestralGens();
		else
			ansGen = m_ancestralGen;
	}

	for (UINT i = 0; i <= ansGen; ++i) {
		pop.useAncestralGen(i);

		// we need info to be in order
		IndInfoIterator traitIt = pop.infoBegin(idx);
		for (IndIterator it = pop.indBegin(); it.valid(); ++it)
			*traitIt++ = qtrait(& * it) ;
	}
	pop.useAncestralGen(0);

	return true;
}


double mapQuanTrait::qtrait(individual * ind)
{
	string key;

	for (vectoru::iterator loc = m_loci.begin(); loc != m_loci.end(); ++loc) {
		// get genotype of ind
		Allele a = ind->allele(*loc, 0);
		Allele b = ind->allele(*loc, 1);

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


double maQuanTrait::qtrait(individual * ind)
{
	UINT index = 0;

	for (vectoru::iterator loc = m_loci.begin(); loc != m_loci.end(); ++loc) {
		// get genotype of ind
		Allele a = ind->allele(*loc, 0);
		Allele b = ind->allele(*loc, 1);

		int numWildtype = 0;

		// count number of wildtype
		if (find(m_wildtype.begin(), m_wildtype.end(), a) != m_wildtype.end() )
			numWildtype++;

		if (find(m_wildtype.begin(), m_wildtype.end(), b) != m_wildtype.end() )
			numWildtype++;
		index = index * 3 + 2 - numWildtype;
	}

	return rng().randNormal(m_qtrait[index], m_sigma[index]);
}


double mlQuanTrait::qtrait(individual * ind)
{
	if (m_mode == QT_Multiplicative) {
		double fit = 1;
		for (vectorop::iterator s = m_qtraits.begin(), sEnd = m_qtraits.end();
		     s != sEnd; ++s)
			fit *= static_cast<quanTrait *>(*s)->qtrait(ind);
		return rng().randNormal(fit, m_sigma);
	} else if (m_mode == QT_Additive) {
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
			m_alleles[j++] = ind->allele(m_loci[i], p);

	return m_func.call(PyObj_As_Double, "(O)", m_numArray);
}


}
