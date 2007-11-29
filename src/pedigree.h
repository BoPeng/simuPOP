/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu
*                                                                         *
*   $LastChangedDate$
*   $Rev$                                                       *
*                                                                         *
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

#ifndef _PEDIGREE_H
#define _PEDIGREE_H
/**
 \file
 \brief head file of class mating and its subclasses
 */
#include "utility.h"

namespace simuPOP {

/// A pedigree manipulation class
class pedigree
{
	typedef vector<vectorlu > Pedigree;
	typedef vector<vectorlu > PedSize;
	typedef vector<vector<double> > PedInfo;

public:
	// numParents can be automatically determined
	pedigree() : m_numParents(0)
	{
	}


	ULONG size(ULONG gen)
	{
		DBG_FAILIF(gen >= m_paternal.size(), IndexError,
		    "Generation number " + toStr(gen) + " out of bound (<"
		    + toStr(m_paternal.size()) + ")");
		return m_paternal[gen].size();
	}


	UINT numParents()
	{
		return m_numParents;
	}


	ULONG gen()
	{
		return m_paternal.size();
	}


	ULONG father(ULONG gen, ULONG idx)
	{
		DBG_FAILIF(gen >= m_paternal.size(), IndexError,
		    "Generation number " + toStr(gen) + " out of bound (<"
		    + toStr(m_paternal.size()) + ")");
		DBG_FAILIF(idx >= m_paternal[gen].size(), IndexError,
		    "Father index out of bound");
		return m_paternal[gen][idx];
	}


	ULONG mother(ULONG gen, ULONG idx)
	{
		DBG_FAILIF(gen >= m_paternal.size(), IndexError,
		    "Generation number " + toStr(gen) + " out of bound (<"
		    + toStr(m_paternal.size()) + ")");
		DBG_FAILIF(idx >= m_maternal[gen].size(), IndexError,
		    "Mother index out of bound");
		return m_maternal[gen][idx];
	}


	ULONG father(ULONG gen, SubPopID subPop, ULONG idx)
	{
		DBG_FAILIF(gen >= m_paternal.size(), IndexError,
		    "Generation number " + toStr(gen) + " out of bound (<"
		    + toStr(m_paternal.size()) + ")");
		DBG_FAILIF(idx >= m_pedSize[gen][subPop], IndexError,
		    "Father index out of bound");
		size_t shift = 0;
		for (SubPopID i = 0; i < subPop; ++i)
			shift += m_pedSize[gen][i];
		return m_paternal[gen][shift + idx];
	}


	ULONG mother(ULONG gen, SubPopID subPop, ULONG idx)
	{
		DBG_FAILIF(gen >= m_paternal.size(), IndexError,
		    "Generation number " + toStr(gen) + " out of bound (<"
		    + toStr(m_paternal.size()) + ")");
		DBG_FAILIF(idx >= m_pedSize[gen][subPop], IndexError,
		    "Mother index out of bound");
		size_t shift = 0;
		for (SubPopID i = 0; i < subPop; ++i)
			shift += m_pedSize[gen][i];
		return m_maternal[gen][shift + idx];
	}


	vectorlu subPopSizes(ULONG gen)
	{
		DBG_FAILIF(gen >= m_paternal.size(), IndexError,
		    "Generation number " + toStr(gen) + " out of bound (<"
		    + toStr(m_paternal.size()) + ")");
		return m_pedSize[gen];
	}


	ULONG subPopSize(ULONG gen, SubPopID subPop)
	{
		DBG_FAILIF(gen >= m_paternal.size(), IndexError,
		    "Generation number " + toStr(gen) + " out of bound (<"
		    + toStr(m_paternal.size()) + ")");
		DBG_FAILIF(static_cast<size_t>(subPop) >= m_pedSize[gen].size(), IndexError,
		    "SubPop index out of bound");
		return m_pedSize[gen][subPop];
	}


	void read(const string & filename, const string & aux_filename = string());

	pedigree * clone()
	{
		return new pedigree();
	}


	/// write the pedigree (and its auxillary information) to files
	void write(const string & filename, const string & aux_filename = string());

	/// choose given individuals from the last generation
	/// The last generation will be shrinked to only have these individuals
	void selectIndividuals(const vectorlu & inds);

	/// mark individuals that are unrelated to the last generation from the pedigree
	void markUnrelated();

	/// remove individuals that are unrelated to the last generation from the pedigree
	/// indexes will be adjusted.
	void removeUnrelated();

private:
	int m_numParents;

	Pedigree m_paternal;
	Pedigree m_maternal;
	PedSize m_pedSize;
	PedInfo m_info;
};

}
#endif
