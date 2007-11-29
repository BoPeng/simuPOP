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
private:
	/// The reason why I do not want to use a
	/// vector<vector< struct<dad, mom, vector<info>>>
	/// structure is because mom and info are optional.
	/// Because pedigree can be large, it makes sense to treat
	/// them separately.

	typedef vector<vectorlu > Pedigree;
	// pedigree subpopulation size information ...
	typedef Pedigree PedSize;
	typedef vector<vector<vectorf> > PedInfo;

public:
	/// create a pedigree. If a filename \c pedfile is given, the pedgree
	/// will be loaded from this file.
	pedigree(const string & pedfile = string());

	/// population size at generation \c gen
	ULONG popSize(ULONG gen)
	{
		DBG_FAILIF(gen >= m_paternal.size(), IndexError,
		    "Generation number " + toStr(gen) + " out of bound (<"
		    + toStr(m_paternal.size()) + ")");
		return m_paternal[gen].size();
	}

	/// return the number of parents for each individual
	UINT numParents()
	{
		return m_numParents;
	}

	/// Set the number of parents for each individual
	void setNumParents(int numParents)
	{
		DBG_ASSERT(numParents == 1 || numParents == 2, ValueError,
			"Number of parents has to be 1 or 2");
		m_numParents = numParents;	
	}
	
	/// Return the number of generations of this pedigree
	ULONG gen()
	{
		return m_paternal.size();
	}

	/// Return the index of the father of individual \c idx at generation \c gen
	/// The returned index is the absolute index of father in the parental generation
	ULONG father(ULONG gen, ULONG idx);

	/// Return the index of the mother of individual \c idx at generation \c gen
	/// The returned index is the absolute index of mother in the parental generation
	ULONG mother(ULONG gen, ULONG idx);

	/// Return the index of the father of individual \c idx of subpopulation
	/// \c subPop at generation \c gen
	/// The returned index is the absolute index of father in the parental generation
	ULONG father(ULONG gen, SubPopID subPop, ULONG idx);

	/// Return the index of the mother of individual \c idx of subpopulation
	/// \c subPop at generation \c gen
	/// The returned index is the absolute index of mother in the parental generation
	ULONG mother(ULONG gen, SubPopID subPop, ULONG idx);

	/// Set the index of the father of individual  \c idx at generation \c gen
	void setFather(ULONG parent, ULONG gen, ULONG idx);

	/// Set the index of the mother of individual  \c idx at generation \c gen
	void setMother(ULONG parent, ULONG gen, ULONG idx);

	/// Set the index of the father of individual \c idx of subpopulation \c subPop at generation \c gen
	void setFather(ULONG parent, ULONG gen, SubPopID subPop, ULONG idx);

	/// Set the index of the mother of individual \c idx of subpopulation \c subPop at generation \c gen
	void setMother(ULONG parent, ULONG gen, SubPopID subPop, ULONG idx);
	
	/// Return information \c name of individual \c idx at generation \c gen
	double info(ULONG gen, ULONG idx, const string & name);

	/// Return information \c name of individual \c idx of subpopulation \c subPop at generation \c gen
	double info(ULONG gen, SubPopID subPop, ULONG idx, const string & name);
	
	/// Set information \c name  of individual \c idx at generation \c gen
	void setInfo(double info, ULONG gen, ULONG idx, const string & name);

	/// Set information \c name of individual \c idx of subpopulation \c subPop at generation \c gen
	void setInfo(double info, ULONG gen, SubPopID subPop, ULONG idx, const string & name);

	/// Add a generation to the existing pedigree, with given subpopulation sizes
	/// \c subPopSize . All parental indexes and information will be set to zero 
	/// for the new generation.
	void addGen(const vectorlu & subPopSize);
	
	/// Return the subpopulation sizes of generation \c gen
	vectorlu subPopSizes(ULONG gen);

	/// Return the subpopulation size of subpopulation \c subPop of generation \c gen
	ULONG subPopSize(ULONG gen, SubPopID subPop);

	/// Make a copy of this pedigree
	pedigree * clone()
	{
		return new pedigree();
	}

	/// load pedigree from a file, the file is usually saved by \c parentTagger or
	/// \c parentsTagger. The format is described in the simuPOP reference manual
	void load(const string & filename);

	/// load information \c name from a information pedigree file and add to this pedigree
	/// Information \c name should not have existed in the pedigree. The information
	/// pedigree file \c filename is usually produced by taggers such as \c sexTagger
	/// \c affectionTagger, \c pyTagger and \c infoTagger.
	void loadInfo(const string & filename, const string & name);

	/// load information \c names from a information pedigree file and add to this pedigree
	/// Information names in \c names should not have existed in the pedigree.
	/// pedigree file \c filename is usually produced by taggers such as \c sexTagger
	/// \c affectionTagger, \c pyTagger and \c infoTagger.
	void loadInfo(const string & filename, const vectorstr & names);

	/// add an information field to the pedigree, with given initial value
	void addInfo(const string & name, double init = 0);

	/// Write the pedigree to a file
	void save(const string & filename);

	/// Save auxiliary information \c name to an information pedigree file
	void saveInfo(const string & filename, const string & name);

	/// save auxiliary information \c names to an information pedigree file
	void saveInfo(const string & filename, const vectorstr & names);

	/// Mark individuals \c inds as 'unrelated' in the last generation
	/// \c removeUnrelated function will remove these individuals from the pdeigree
	void selectIndividuals(const vectorlu & inds);

	/// mark individuals that are unrelated to the visible individuals at the
	/// last generation (not marked by selectIndividuals) from the pedigree.
	void markUnrelated();

	/// remove individuals that are unrelated to the last generation from the pedigree
	/// indexes will be adjusted.
	void removeUnrelated();

private:
	int m_numParents;

	Pedigree m_paternal;
	///
	Pedigree m_maternal;
	///
	PedSize m_pedSize;
	///
	PedInfo m_info;
	///
	vectorstr m_infoNames;
};

}
#endif
