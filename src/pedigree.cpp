/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu                                                        *
*                                                                         *
*   $LastChangedDate: 2006-02-21 15:27:25 -0600 (Tue, 21 Feb 2006)        *
*   $Rev: 191$                                                            *
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

#include "pedigree.h"
#include <fstream>
using std::ifstream;
using std::ofstream;

#include <sstream>
using std::istringstream;
using std::ws;

#include <numeric>
using std::accumulate;

namespace simuPOP {

const unsigned long UnusedIndividual = std::numeric_limits<unsigned long>::max();

pedigree::pedigree(const string & pedfile)
	: m_numParents(0)
{
	if (!pedfile.empty())
		load(pedfile);
}


ULONG pedigree::father(ULONG gen, ULONG idx)
{
	DBG_FAILIF(gen >= m_paternal.size(), IndexError,
	    "Generation number " + toStr(gen) + " out of bound (<"
	    + toStr(m_paternal.size()) + ")");
	DBG_FAILIF(idx >= m_paternal[gen].size(), IndexError,
	    "Father index out of bound");
	return m_paternal[gen][idx];
}


ULONG pedigree::mother(ULONG gen, ULONG idx)
{
	DBG_FAILIF(gen >= m_paternal.size(), IndexError,
	    "Generation number " + toStr(gen) + " out of bound (<"
	    + toStr(m_paternal.size()) + ")");
	DBG_FAILIF(idx >= m_maternal[gen].size(), IndexError,
	    "Mother index out of bound");
	return m_maternal[gen][idx];
}


ULONG pedigree::father(ULONG gen, SubPopID subPop, ULONG idx)
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


ULONG pedigree::mother(ULONG gen, SubPopID subPop, ULONG idx)
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

void pedigree::setFather(ULONG parent, ULONG gen, ULONG idx)
{
	DBG_FAILIF(m_numParents == 0, ValueError,
		"Please set number of parents before change the pedigree");
	DBG_FAILIF(gen >= m_paternal.size(), IndexError,
	    "Generation number " + toStr(gen) + " out of bound (<"
	    + toStr(m_paternal.size()) + ")");
	DBG_FAILIF(idx >= m_paternal[gen].size(), IndexError,
	    "Father index out of bound");
	DBG_FAILIF(gen >= 1 && parent >= m_paternal[gen-1].size(),
		IndexError, "Given parent index is out of bound of the parental population");
	m_paternal[gen][idx] = parent;
}


void pedigree::setMother(ULONG parent, ULONG gen, ULONG idx)
{
	DBG_FAILIF(m_numParents == 0, ValueError,
		"Please set number of parents before change the pedigree");
		
	DBG_FAILIF(gen >= m_paternal.size(), IndexError,
	    "Generation number " + toStr(gen) + " out of bound (<"
	    + toStr(m_paternal.size()) + ")");
	DBG_FAILIF(idx >= m_maternal[gen].size(), IndexError,
	    "Mother index out of bound");
	DBG_FAILIF(gen >= 1 && parent >= m_maternal[gen-1].size(),
		IndexError, "Given parent index is out of bound of the parental population");
	m_maternal[gen][idx] = parent;
}


void pedigree::setFather(ULONG parent, ULONG gen, SubPopID subPop, ULONG idx)
{
	DBG_FAILIF(m_numParents == 0, ValueError,
		"Please set number of parents before change the pedigree");
	DBG_FAILIF(gen >= m_paternal.size(), IndexError,
	    "Generation number " + toStr(gen) + " out of bound (<"
	    + toStr(m_paternal.size()) + ")");
	DBG_FAILIF(idx >= m_pedSize[gen][subPop], IndexError,
	    "Father index out of bound");
	DBG_FAILIF(gen >= 1 && parent >= m_paternal[gen-1].size(),
		IndexError, "Given parent index is out of bound of the parental population");
	size_t shift = 0;
	for (SubPopID i = 0; i < subPop; ++i)
		shift += m_pedSize[gen][i];
	m_paternal[gen][shift + idx] = parent;
}


void pedigree::setMother(ULONG parent, ULONG gen, SubPopID subPop, ULONG idx)
{
	DBG_FAILIF(m_numParents == 0, ValueError,
		"Please set number of parents before change the pedigree");
	DBG_FAILIF(gen >= m_paternal.size(), IndexError,
	    "Generation number " + toStr(gen) + " out of bound (<"
	    + toStr(m_paternal.size()) + ")");
	DBG_FAILIF(idx >= m_pedSize[gen][subPop], IndexError,
	    "Mother index out of bound");
	DBG_FAILIF(gen >= 1 && parent >= m_maternal[gen-1].size(),
		IndexError, "Given parent index is out of bound of the parental population");
	size_t shift = 0;
	for (SubPopID i = 0; i < subPop; ++i)
		shift += m_pedSize[gen][i];
	m_maternal[gen][shift + idx] = parent;
}


double pedigree::info(ULONG gen, ULONG idx, const string & name)
{
	DBG_FAILIF(gen >= m_paternal.size(), IndexError,
	    "Generation number " + toStr(gen) + " out of bound (<"
	    + toStr(m_paternal.size()) + ")");
	DBG_FAILIF(idx >= m_paternal[gen].size(), IndexError,
	    "Father index out of bound");
	for (size_t i = 0; i < m_infoNames.size(); ++i)
		if (m_infoNames[i] == name)
			return m_info[gen][idx][i];
	DBG_FAILIF(true, ValueError, "Information " + name + " is not found");
	return 0;
}


double pedigree::info(ULONG gen, SubPopID subPop, ULONG idx, const string & name)
{
	DBG_FAILIF(gen >= m_paternal.size(), IndexError,
	    "Generation number " + toStr(gen) + " out of bound (<"
	    + toStr(m_paternal.size()) + ")");
	DBG_FAILIF(idx >= m_pedSize[gen][subPop], IndexError,
	    "Father index out of bound");
	size_t shift = 0;
	for (SubPopID i = 0; i < subPop; ++i)
		shift += m_pedSize[gen][i];
	for (size_t i = 0; i < m_infoNames.size(); ++i)
		if (m_infoNames[i] == name)
			return m_info[gen][shift + idx][i];
	DBG_FAILIF(true, ValueError, "Information " + name + " is not found");
	return 0;
}


void pedigree::setInfo(double info, ULONG gen, ULONG idx, const string & name)
{
	DBG_FAILIF(m_numParents == 0, ValueError,
		"Please set number of parents before change the pedigree");
	DBG_FAILIF(gen >= m_paternal.size(), IndexError,
	    "Generation number " + toStr(gen) + " out of bound (<"
	    + toStr(m_paternal.size()) + ")");
	DBG_FAILIF(idx >= m_paternal[gen].size(), IndexError,
	    "Father index out of bound");
	for (size_t i = 0; i < m_infoNames.size(); ++i)
		if (m_infoNames[i] == name) {
			m_info[gen][idx][i] = info;
			return;
		}
	DBG_FAILIF(true, ValueError, "Information " + name + " is not found");
}


void pedigree::setInfo(double info, ULONG gen, SubPopID subPop, ULONG idx, const string & name)
{
	DBG_FAILIF(m_numParents == 0, ValueError,
		"Please set number of parents before change the pedigree");
	DBG_FAILIF(gen >= m_paternal.size(), IndexError,
	    "Generation number " + toStr(gen) + " out of bound (<"
	    + toStr(m_paternal.size()) + ")");
	DBG_FAILIF(idx >= m_pedSize[gen][subPop], IndexError,
	    "Father index out of bound");
	size_t shift = 0;
	for (SubPopID i = 0; i < subPop; ++i)
		shift += m_pedSize[gen][i];
	for (size_t i = 0; i < m_infoNames.size(); ++i)
		if (m_infoNames[i] == name) {
			m_info[gen][shift + idx][i] = info;
			return;
		}
	DBG_FAILIF(true, ValueError, "Information " + name + " is not found");
}

vectorlu pedigree::subPopSizes(ULONG gen)
{
	DBG_FAILIF(gen >= m_paternal.size(), IndexError,
	    "Generation number " + toStr(gen) + " out of bound (<"
	    + toStr(m_paternal.size()) + ")");
	return m_pedSize[gen];
}


ULONG pedigree::subPopSize(ULONG gen, SubPopID subPop)
{
	DBG_FAILIF(gen >= m_paternal.size(), IndexError,
	    "Generation number " + toStr(gen) + " out of bound (<"
	    + toStr(m_paternal.size()) + ")");
	DBG_FAILIF(static_cast<size_t>(subPop) >= m_pedSize[gen].size(), IndexError,
	    "SubPop index out of bound");
	return m_pedSize[gen][subPop];
}

void pedigree::addGen(const vectorlu & sizes)
{
	DBG_FAILIF(m_numParents == 0, ValueError,
		"Please set number of parents before change the pedigree");
	
	ULONG popSize = accumulate(sizes.begin(), sizes.end(), 0UL);
	m_paternal.push_back(vectorlu(popSize));
	if (m_numParents == 2)
		m_maternal.push_back(vectorlu(popSize));
	if (m_info.size() > 0) {
		UINT infoSize = m_infoNames.size();
		m_info.push_back(vector<vectorf>(popSize));
		for (size_t i = 0; i < popSize; ++i)
			m_info.back()[i].resize(infoSize);
	}
	m_pedSize.push_back(sizes);
}

void pedigree::load(const string & filename)
{
	m_paternal.clear();
	m_maternal.clear();
	m_pedSize.clear();

	ifstream ifs(filename.c_str());

	DBG_FAILIF(!ifs, SystemError, "Can not open pedigree file " + filename + " to read");

	string line;
	while (getline(ifs, line)) {
		istringstream input(line);
		long int idx;
		vectorlu values;
		vectorlu sizes;
		while (!input.eof()) {
			input >> idx;
			values.push_back(idx == -1 ? UnusedIndividual : static_cast<ULONG>(idx));
			input >> ws;
			if (input.peek() == '#') {
				// ignore '#'
				input.ignore();
				// start reading population size.
				while (!input.eof()) {
					input >> idx;
					sizes.push_back(idx);
					input >> ws;
				}
				// exit the outer loop
				break;
			}
		}
		// check if things are all right
		// determine number of parents...
		ULONG popSize = accumulate(sizes.begin(), sizes.end(), 0UL);
		DBG_FAILIF(popSize != 0 && values.size() == 0,
		    ValueError, "No parent is read");
		if (popSize == 0 && values.size() == 0)
			continue;
		if (m_numParents == 0)
			m_numParents = values.size() / popSize;
		DBG_ASSERT(m_numParents * popSize == values.size(), ValueError,
		    "Number of parents does not match subpopulation sizes.\n"
		    "Line: " + toStr(m_paternal.size() + 1) + ", Individuals read: "
		    + toStr(values.size()) +
		    ", Pop size: " + toStr(popSize));
		m_pedSize.push_back(vectorlu());
		m_pedSize.back().swap(sizes);
		if (m_numParents == 1) {
			m_paternal.push_back(vectorlu());
			m_paternal.back().swap(values);
		} else if (m_numParents == 2) {
			m_paternal.push_back(vectorlu(popSize));
			m_maternal.push_back(vectorlu(popSize));
			for (size_t i = 0; i < popSize; ++i) {
				m_paternal.back()[i] = values[2 * i];
				m_maternal.back()[i] = values[2 * i + 1];
			}
		} else {
			DBG_ASSERT(false, SystemError,
			    "Sorry, pedigree does not support more than two parents");
		}
	}
	ifs.close();
}


void pedigree::loadInfo(const string & filename, const string & name)
{
	vectorstr names(1, name);
	loadInfo(filename, names);
}


void pedigree::loadInfo(const string & filename, const vectorstr & names)
{
	ifstream afs(filename.c_str());

	DBG_FAILIF(!afs, SystemError, "Can not open auxiliary information pedigree" + filename + " to read");
	size_t gen = 0;
	size_t numInfo = names.size();

	for (size_t i = 0; i < numInfo; ++i) {
		DBG_ASSERT(find(m_infoNames.begin(), m_infoNames.end(), names[i]) == m_infoNames.end(),
		    ValueError, "Information " + names[i] + " has already been loaded");
		m_infoNames.push_back(names[i]);
	}
	string line;
	while (getline(afs, line)) {
		istringstream input(line);
		vectorf values;
		double value;
		while (!input.eof()) {
			input >> value;
			values.push_back(value);
			input >> ws;
		}
		DBG_FAILIF(gen >= m_paternal.size(), ValueError,
		    "Information pedigree is larger than parental pedigree");
		ULONG size = m_paternal[gen].size();

		DBG_FAILIF(numInfo * size != values.size(), ValueError,
		    "At generation " + toStr(gen) + ", number of information read is "
		    + toStr(values.size()) + ", which is not a multiple of number of individuals "
		    + toStr(m_paternal[gen].size()));
		//
		if (m_info.size() <= gen)
			m_info.push_back(vector<vectorf>(size));
		DBG_FAILIF(m_info.size() < gen + 1, ValueError,
		    "Error loading information pedigree");
		size_t idx = 0;
		for (size_t i = 0; i < size; ++i)
			for (size_t j = 0; j < numInfo; ++j, ++idx)
				m_info[gen][i].push_back(values[idx]);
		//
		gen++;
	}
	afs.close();
}


void pedigree::addInfo(const string & name, double init)
{
	DBG_ASSERT(find(m_infoNames.begin(), m_infoNames.end(), name) == m_infoNames.end(),
	    ValueError, "Information " + name + " has already been loaded");
	m_infoNames.push_back(name);

	for (size_t gen = 0; gen < m_paternal.size(); ++gen) {
		ULONG size = m_paternal[gen].size();
		if (m_info.size() <= gen)
			m_info.push_back(vector<vectorf>(size));
		DBG_FAILIF(m_info.size() < gen + 1, ValueError,
		    "Error loading information pedigree");
		for (size_t i = 0; i < size; ++i)
			m_info[gen][i].push_back(init);
	}
}


void pedigree::save(const string & filename)
{
	ofstream ofs(filename.c_str());

	DBG_FAILIF(!ofs, SystemError, "Can not open pedigree file " + filename + " to write.");

	for (size_t gen = 0; gen < m_paternal.size(); ++gen) {
		size_t sz = m_paternal[gen].size();
		for (size_t idx = 0; idx < sz; ++idx) {
			if (m_paternal[gen][idx] == UnusedIndividual)
				ofs << -1;
			else
				ofs << m_paternal[gen][idx];
			if (m_numParents == 2) {
				if (m_maternal[gen][idx] == UnusedIndividual)
					ofs << "\t-1";
				else
					ofs << '\t' << m_maternal[gen][idx];
			}
			ofs << '\t';
		}
		ofs << "#\t";
		for (size_t idx = 0; idx < m_pedSize[gen].size(); ++idx)
			ofs << m_pedSize[gen][idx] << '\t';
		ofs << '\n';
	}
	ofs.close();
}


void pedigree::saveInfo(const string & filename, const string & name)
{
	vectorstr names(1, name);

	saveInfo(filename, names);
}


/// save auxiliary information \c name to a file
void pedigree::saveInfo(const string & filename, const vectorstr & names)
{
	vectoru idx;

	for (size_t i = 0; i < names.size(); ++i) {
		size_t j = 0;
		for (; j < m_infoNames.size(); ++j) {
			if (m_infoNames[j] == names[i]) {
				idx.push_back(j);
				break;
			}
		}
		DBG_FAILIF(j == m_infoNames.size(), ValueError, "Invalid information name: " + names[i]);
	}

	ofstream afs(filename.c_str());
	DBG_FAILIF(!afs, SystemError, "Can not open information pedigree file " + filename + " to write.");

	for (size_t gen = 0; gen < m_paternal.size(); ++gen) {
		size_t sz = m_paternal[gen].size();
		for (size_t i = 0; i < sz; ++i)
			for (size_t j = 0; j < idx.size(); ++j)
				afs << m_info[gen][i][j] << '\t';
		afs << '\n';
	}
	afs.close();
}


void pedigree::selectIndividuals(const vectorlu & inds)
{
	DBG_FAILIF(m_paternal.empty(), ValueError,
	    "Can not select individuals from an empty pedigree");

	vector<bool> used(m_paternal.back().size(), false);
	vectorlu::const_iterator it = inds.begin();
	vectorlu::const_iterator it_end = inds.end();
	for (; it != it_end; ++it) {
		DBG_FAILIF(*it >= m_paternal.back().size(), IndexError,
		    "Index exceeded the size of the last generation");
		used[*it] = true;
	}
	for (size_t idx = 0; idx < m_paternal.back().size(); ++idx)
		if (!used[idx]) {
			m_paternal.back()[idx] = UnusedIndividual;
			m_maternal.back()[idx] = UnusedIndividual;
		}
}


void pedigree::markUnrelated()
{
	if (m_paternal.size() <= 1)
		return;

	vector<bool> used;
	// starting from the last generation, gen=0 etc will be replaced.
	for (size_t gen = m_paternal.size() - 1; gen > 0; --gen) {
		used.clear();
		used.resize(m_paternal[gen - 1].size(), false);
		for (size_t idx = 0; idx < m_paternal[gen].size(); ++idx)
			if (m_paternal[gen][idx] != UnusedIndividual)
				used[m_paternal[gen][idx]] = true;
		if (m_numParents == 2)
			for (size_t idx = 0; idx < m_maternal[gen].size(); ++idx)
				if (m_maternal[gen][idx] != UnusedIndividual)
					used[m_maternal[gen][idx]] = true;
		for (size_t idx = 0; idx < m_paternal[gen - 1].size(); ++idx)
			if (!used[idx]) {
				m_paternal[gen - 1][idx] = UnusedIndividual;
				if (m_numParents == 2)
					m_maternal[gen - 1][idx] = UnusedIndividual;
			}
	}
}


void pedigree::removeUnrelated()
{
	if (m_paternal.size() <= 1)
		return;

	// starting from the last generation, gen=0 etc will be replaced.
	for (size_t gen = 0; gen < m_paternal.size() - 1; ++gen) {
		vectorlu & curGen = m_paternal[gen];
		vectorlu & nextGen = m_paternal[gen + 1];
		size_t shift = 0;
		for (size_t idx = 0; idx < curGen.size(); ++idx)
			if (curGen[idx] == UnusedIndividual) {
				// the next generation, with value > idx will be shifted by 1.
				for (size_t idx1 = 0; idx1 < nextGen.size(); ++idx1)
					if (nextGen[idx1] != UnusedIndividual && nextGen[idx1] + shift > idx)
						nextGen[idx1]--;
				shift++;
			}
	}
	// maternal
	if (m_numParents == 2) {
		for (size_t gen = 0; gen < m_maternal.size() - 1; ++gen) {
			vectorlu & curGen = m_maternal[gen];
			vectorlu & nextGen = m_maternal[gen + 1];
			size_t shift = 0;
			for (size_t idx = 0; idx < curGen.size(); ++idx)
				if (curGen[idx] == UnusedIndividual) {
					// the next generation, with value > idx will be shifted by 1.
					for (size_t idx1 = 0; idx1 < nextGen.size(); ++idx1)
						if (nextGen[idx1] != UnusedIndividual && nextGen[idx1] + shift > idx)
							nextGen[idx1]--;
					shift++;
				}
		}
	}
	// adjust m_pedSize
	for (size_t gen = 0; gen < m_paternal.size(); ++gen) {
		size_t idx = 0;
		for (UINT sp = 0; sp < m_pedSize[gen].size(); ++sp) {
			UINT spSize = m_pedSize[gen][sp];
			for (size_t i = 0; i < spSize; ++i, ++idx)
				if (m_paternal[gen][idx] == UnusedIndividual)
					m_pedSize[gen][sp]--;
		}
	}
	// remove individuals
	// new pedigree generation entries
	vectorlu l_pat;
	vectorlu l_mat;
	vector<vectorf> l_info;
	for (size_t gen = 0; gen < m_paternal.size(); ++gen) {
		l_pat.clear();
		l_mat.clear();
		l_info.clear();
		for (size_t idx = 0; idx < m_paternal[gen].size(); ++idx)
			if (m_paternal[gen][idx] != UnusedIndividual) {
				l_pat.push_back(m_paternal[gen][idx]);
				if (m_numParents == 2) {
					DBG_ASSERT(m_maternal[gen][idx] != UnusedIndividual,
					    ValueError, "Inconsistent maternal and matermal pedigree");
					l_mat.push_back(m_maternal[gen][idx]);
				}
				if (!m_info.empty())
					l_info.push_back(m_info[gen][idx]);
			}
		m_paternal[gen].swap(l_pat);
		if (m_numParents == 2)
			m_maternal[gen].swap(l_mat);
		if (!m_info.empty())
			m_info[gen].swap(l_info);
	}
}


}
