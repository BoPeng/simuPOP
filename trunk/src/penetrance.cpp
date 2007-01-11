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

#include "penetrance.h"

namespace simuPOP
{
	double pyPenetrance::penet(individual * ind)
	{
		int len = m_loci.size() * ind->ploidy();
		if( m_len != len )
		{
			m_len = len;
			m_alleles.resize(m_len);
			if(m_numArray != NULL)
				Py_DECREF(m_numArray);
#ifdef SIMUMPI
			m_numArray = Allele_Vec_As_NumArray( m_alleles.begin(), m_alleles.end(),
				m_alleles.size(), m_alleles.size(), 0, m_alleles.size());
#else
			m_numArray = Allele_Vec_As_NumArray( m_alleles.begin(), m_alleles.end() );
#endif
		}

		UINT pEnd = ind->ploidy();
		for(size_t i=0, iEnd=m_loci.size(), j=0; i < iEnd; ++i)
			for(UINT p=0; p < pEnd; ++p)
				m_alleles[j++] = ind->allele(m_loci[i], p);

		double resDouble;
		PyCallFunc(m_func, "(O)", m_numArray, resDouble, PyObj_As_Double);

		// make sure the returned value is legitimate.
		DBG_ASSERT( fcmp_ge( resDouble, 0.) && fcmp_le( resDouble, 1.),
			ValueError, "Returned fitness " + toStr(resDouble) + " is out of range [0,1]" );

		return resDouble;
	}

}
