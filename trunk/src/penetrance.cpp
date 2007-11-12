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

#include "penetrance.h"

namespace simuPOP {
/// set pentrance to all individuals and record penetrance if requested.
bool penetrance::apply(population & pop)
{
	double p;

	bool savePene = infoSize() > 0;

	UINT ansGen = 0;

	if (m_ancestralGen == -1)
		ansGen = pop.ancestralDepth();
	else if (m_ancestralGen > 0) {
		if (static_cast<UINT>(m_ancestralGen) > pop.ancestralDepth())
			ansGen = pop.ancestralDepth();
		else
			ansGen = m_ancestralGen;
	}

	for (UINT i = 0; i <= ansGen; ++i) {
		pop.useAncestralPop(i);
		IndInfoIterator penIt;
		if (savePene) {
			UINT idx = pop.infoIdx(infoField(0));
			penIt = pop.infoBegin(idx);
		}
		for (IndIterator it = pop.indBegin(); it.valid(); ++it) {
			p = penet(& * it);

			if (rng().randUniform01() < p)
				it->setAffected(true);
			else
				it->setAffected(false);
			if (savePene)
				*penIt++ = p;
		}
	}
	pop.useAncestralPop(0);

	return true;
}


bool penetrance::applyDuringMating(population & pop, IndIterator offspring,
                                   individual * dad, individual * mom)
{
	double p = penet(& * offspring);

	if (infoSize() > 0)
		(*offspring).setInfo(p, infoField(0));
	if (rng().randUniform01() < p)
		offspring->setAffected(true);
	else
		offspring->setAffected(false);
	return true;
}


double mapPenetrance::penet(individual * ind)
{
	string key;

	for (vectoru::iterator loc = m_loci.begin(); loc != m_loci.end(); ++loc) {

		/// get genotype of ind
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
	    "No penetrance value for genotype " + key);

	return pos->second;
}


double maPenetrance::penet(individual * ind)
{
	UINT index = 0;

	for (vectoru::iterator loc = m_loci.begin(); loc != m_loci.end(); ++loc) {

		/// get genotype of ind
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

	return m_penetrance[index];
}


double mlPenetrance::penet(individual * ind)
{
	if (m_mode == PEN_Multiplicative) {
		// x1 x2 x3 ...
		double pen = 1;
		for (vectorop::iterator s = m_peneOps.begin(), sEnd = m_peneOps.end();
		     s != sEnd; ++s)
			pen *= static_cast<penetrance *>(*s)->penet(ind);
		return pen;
	} else if (m_mode == PEN_Additive) {
		// x1 + x2 + x3
		double pen = 0;
		for (vectorop::iterator s = m_peneOps.begin(), sEnd = m_peneOps.end();
		     s != sEnd; ++s)
			pen += static_cast<penetrance *>(*s)->penet(ind);
		return pen > 1 ? 1 : pen;
	} else if (m_mode == PEN_Heterogeneity) {
		// 1-(1-x1)(1-x2)
		double pen = 1;
		for (vectorop::iterator s = m_peneOps.begin(), sEnd = m_peneOps.end();
		     s != sEnd; ++s)
			pen *= 1 - static_cast<penetrance *>(*s)->penet(ind);
		return 1 - pen;
	}

	return 0.0;
}


double pyPenetrance::penet(individual * ind)
{
	int len = m_loci.size() * ind->ploidy();

	if (m_len != len) {
		m_len = len;
		m_alleles.resize(m_len);
		if (m_numArray != NULL)
			Py_DECREF(m_numArray);
#ifdef SIMUMPI
		m_numArray = Allele_Vec_As_NumArray(m_alleles.begin(), m_alleles.end(),
		                 m_alleles.size(), m_alleles.size(), 0, m_alleles.size());
#else
		m_numArray = Allele_Vec_As_NumArray(m_alleles.begin(), m_alleles.end() );
#endif
	}

	if (infoSize() > 1) {
		if (m_info.size() + 1 != infoSize() ) {
			m_info.resize(infoSize() - 1);
			m_infoArray = Double_Vec_As_NumArray(m_info.begin(), m_info.end());
		}
		// assign information fields from individusl
		for (size_t i = 1; i < infoSize(); ++i)
			m_info[i - 1] = ind->info(infoField(i));
	}

	UINT pEnd = ind->ploidy();
	for (size_t i = 0, iEnd = m_loci.size(), j = 0; i < iEnd; ++i)
		for (UINT p = 0; p < pEnd; ++p)
			m_alleles[j++] = ind->allele(m_loci[i], p);

	double resDouble;
	if (infoSize() <= 1) {
		PyCallFunc(m_func, "(O)", m_numArray, resDouble, PyObj_As_Double);
	} else {
		PyCallFunc2(m_func, "(OO)", m_numArray, m_infoArray, resDouble, PyObj_As_Double);
	}

	// make sure the returned value is legitimate.
	DBG_ASSERT(fcmp_ge(resDouble, 0.) && fcmp_le(resDouble, 1.),
	    ValueError, "Returned fitness " + toStr(resDouble) + " is out of range [0,1]");

	return resDouble;
}


}
