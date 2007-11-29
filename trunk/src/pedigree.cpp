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

void pedigree::read(const string & filename, const string & aux_filename)
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

	if (aux_filename.empty())
		return;

	ifstream afs(aux_filename.c_str());
	DBG_FAILIF(!afs, SystemError, "Can not open auxiliary information pedigree" + aux_filename + " to read");
	size_t gen = 0;
	while (getline(afs, line)) {
		m_info.push_back(vector<double>());

		istringstream values(line);
		double value;
		while (!values.eof()) {
			values >> value;
			m_info.back().push_back(value);
		}
		if (!m_info.empty() && m_info.back().empty())
			m_info.pop_back();
		DBG_FAILIF(m_info.size() > m_maternal.size(), ValueError,
		    "Auxiliary information pedigree has different size from the main pedigree");
		DBG_FAILIF(m_info.back().size() != m_maternal[gen].size(), ValueError,
		    "Auxiliary information pedigree has different size from the main pedigree");
		gen++;
	}
	afs.close();
}


void pedigree::write(const string & filename, const string & aux_filename)
{
	ofstream ofs(filename.c_str());

	DBG_FAILIF(!ofs, SystemError, "Can not open pedigree file " + filename + " to write.");

	bool hasInfo = !aux_filename.empty();
	ofstream afs;
	if (hasInfo) {
		afs.open(aux_filename.c_str());
		DBG_FAILIF(!ofs, SystemError, "Can not open auxiliary information file " + filename + " to write.");
	}

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
			if (hasInfo) {
				afs << m_info[gen][idx];
				if (idx != sz - 1)
					afs << '\t';
			}
		}
		ofs << "#\t";
		for (size_t idx = 0; idx < m_pedSize[gen].size(); ++idx)
			ofs << m_pedSize[gen][idx] << '\t';
		ofs << '\n';
		if (hasInfo)
			afs << '\n';
	}
	ofs.close();
	if (hasInfo)
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
	vector<double> l_info;
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
