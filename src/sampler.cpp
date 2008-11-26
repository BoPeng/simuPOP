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

#include "sampler.h"

namespace simuPOP {


// 
// void sample::resetParentalIndex(population & pop, const string & fatherField,
//                                 const string & motherField, const string & indexField)
// {
// 	UINT fatherIdx = pop.infoIdx(fatherField);
// 	UINT motherIdx = pop.infoIdx(motherField);
// 	UINT indexIdx = pop.infoIdx(indexField);
// 
// 	// for the top generation, no parents
// 	pop.useAncestralGen(pop.ancestralGens());
// 	for (IndIterator it = pop.indBegin(); it.valid(); ++it) {
// 		// -1 means no parents.
// 		// we do not use 0 since 0 is valid
// 		it->setInfo(-1, fatherIdx);
// 		it->setInfo(-1, motherIdx);
// 	}
// 	// for other generations
// 	for (size_t ans = 0; ans < pop.ancestralGens(); ++ans) {
// 		// parents... get their old index
// 		pop.useAncestralGen(ans + 1);
// 		vectorf oldindex = pop.indInfo(indexIdx);
// 		// children
// 		pop.useAncestralGen(ans);
// 		for (IndIterator it = pop.indBegin(); it.valid(); ++it) {
// 			// what is the idx of my parents now?
// 			vectorf::iterator tmp = find(oldindex.begin(), oldindex.end(), it->info(fatherIdx));
// 			// no parents
// 			if (tmp == oldindex.end())
// 				it->setInfo(-1, fatherIdx);
// 			else
// 				it->setInfo(tmp - oldindex.begin(), fatherIdx);
// 			tmp = find(oldindex.begin(), oldindex.end(), it->info(motherIdx));
// 			if (tmp == oldindex.end())
// 				it->setInfo(-1, motherIdx);
// 			else
// 				it->setInfo(tmp - oldindex.begin(), motherIdx);
// 		}
// 	}
// 	pop.useAncestralGen(0);
// }


// bool caseControlSample::prepareSample(population & pop)
// {
// 	// record the index in the old population, for samples.
// 	pop.addInfoField("oldindex", -1);
// 	pop.locateRelatives(REL_Self, vectorstr(1, "oldindex"));
// 	if (!m_spSample) {                                                        // sample from the whole population.
// 		DBG_FAILIF(m_numCases.size() > 1 || m_numControls.size() > 1,
// 			ValueError, "Cases, controls need to be a number if sample from the whole population.");
// 
// 		// first get number of affected.
// 		int numAffected = count_if(pop.indBegin(), pop.indEnd(),
// 			isAffected<individual>());
// 		int numUnaffected = pop.popSize() - numAffected;
// 
// 		if (m_numCases.size() == 1 && m_numCases[0] > numAffected)
// 			cout << "Warning: Not enough affected individuals to be sampled: " <<
// 			"expected " << m_numCases[0] << " available: " << numAffected;
// 
// 		if (m_numControls.size() == 1 && m_numControls[0] > numUnaffected)
// 			cout << "Warning: Not enough unaffected individuals to be sampled." <<
// 			"expected " << m_numControls[0] << " available: " << numUnaffected;
// 
// 		// save individual indexes
// 		m_caseIdx.resize(1);
// 		m_controlIdx.resize(1);
// 		m_caseIdx[0].resize(numAffected);
// 		m_controlIdx[0].resize(numUnaffected);
// 		for (size_t i = 0, j = 0, k = 0, iEnd = pop.popSize(); i < iEnd; ++i) {
// 			if (pop.ind(i).affected() )
// 				m_caseIdx[0][j++] = i;
// 			else
// 				m_controlIdx[0][k++] = i;
// 		}
// 	} else {
// 		UINT numSP = pop.numSubPop();
// 		vectori numAffected(numSP);
// 		vectori numUnaffected(numSP);
// 		m_caseIdx.resize(numSP);
// 		m_controlIdx.resize(numSP);
// 
// 		if (m_numCases.size() != numSP || m_numControls.size() != numSP)
// 			throw ValueError("Size of cases/controls does not match number of subpopulations.");
// 
// 		for (UINT sp = 0; sp < numSP; ++sp) {
// 			// first get number of affected.
// 			numAffected[sp] = count_if(pop.indBegin(sp), pop.indEnd(sp),
// 				isAffected<individual>());
// 			numUnaffected[sp] = pop.subPopSize(sp) - numAffected[sp];
// 
// 			if (m_numCases[sp] > numAffected[sp])
// 				cout << "Warning: Not enough affected individuals to be sampled: " <<
// 				"expected " << m_numCases[sp] << " available: " << numAffected[sp];
// 
// 			if (m_numControls[sp] > numUnaffected[sp])
// 				cout << "Warning: Not enough unaffected individuals to be sampled." <<
// 				"expected " << m_numControls[sp] << " available: " << numUnaffected[sp];
// 
// 			// save indexes
// 			m_caseIdx[sp].resize(numAffected[sp]);
// 			m_controlIdx[sp].resize(numUnaffected[sp]);
// 			for (size_t i = 0, j = 0, k = 0, iEnd = pop.subPopSize(sp); i < iEnd; ++i) {
// 				if (pop.ind(i, sp).affected() )
// 					m_caseIdx[sp][j++] = i;
// 				else
// 					m_controlIdx[sp][k++] = i;
// 			}
// 		}
// 	}
// 	return true;
// }
// 
// 
// population & caseControlSample::drawsample(population & pop)
// {
// 	if (!m_spSample) {                                                        // draw sample from the whole population
// 		DBG_DO(DBG_SELECTOR, cout << "Selecting from the whole population" << endl);
// 		// now choose m_caseIdx and m_controlIdx
// 		// random shuffle the index array.
// 		random_shuffle(m_caseIdx[0].begin(), m_caseIdx[0].end());
// 		random_shuffle(m_controlIdx[0].begin(), m_controlIdx[0].end());
// 
// 		int numAffected = m_caseIdx[0].size();
// 		int numUnaffected = m_controlIdx[0].size();
// 
// 		// keep first m_size individuals of shuffled indexes
// 		int nCase, nControl;
// 		if (m_numCases.empty() || m_numCases[0] == 0 || m_numCases[0] > numAffected)
// 			nCase = numAffected;
// 		else
// 			nCase = m_numCases[0];
// 
// 		if (m_numControls.empty() || m_numControls[0] == 0 || m_numControls[0] > numUnaffected)
// 			nControl = numUnaffected;
// 		else
// 			nControl = m_numControls[0];
// 
// 		DBG_DO(DBG_SELECTOR, cout << "nCase: " << nCase << " nControl: " << nControl << endl);
// 
// 		// keep track of which how many case/control from each subpop
// 		int i;
// 		vectori nCaseInSP(pop.numSubPop()), nControlInSP(pop.numSubPop());
// 		for (i = 0; i < nCase; ++i) {
// 			nCaseInSP[pop.subPopIndPair(m_caseIdx[0][i]).first]++;
// 			pop.ind(m_caseIdx[0][i]).setSubPopID(0);
// 		}
// 		// remove others
// 		for (i = nCase; i < numAffected; ++i)
// 			pop.ind(m_caseIdx[0][i]).setSubPopID(-1);
// 
// 		// keep first m_size individuals of shuffled indexes
// 		for (i = 0; i < nControl; ++i) {
// 			nControlInSP[pop.subPopIndPair(m_controlIdx[0][i]).first]++;
// 			pop.ind(m_controlIdx[0][i]).setSubPopID(1);
// 		}
// 		// remove others
// 		for (i = nControl; i < numUnaffected; ++i)
// 			pop.ind(m_controlIdx[0][i]).setSubPopID(-1);
// 
// 		DBG_DO(DBG_SELECTOR, cout << "Getting sample population" << endl);
// 		population & sample = pop.newPopByIndID(0);
// 
// 		sample.setIntVectorVar("nCases", nCaseInSP);
// 		sample.setIntVectorVar("nControls", nControlInSP);
// 		// determine exactly how many cases and controls in the final sample
// 		return sample;
// 	} else {                                                                          // sample from each subpop
// 		UINT numSP = pop.numSubPop();
// 		int nCase, nControl;
// 
// 		for (UINT sp = 0; sp < numSP; ++sp) {
// 			// now choose m_caseIdx and m_controlIdx
// 			// random shuffle the index array.
// 			random_shuffle(m_caseIdx[sp].begin(), m_caseIdx[sp].end());
// 			random_shuffle(m_controlIdx[sp].begin(), m_controlIdx[sp].end());
// 
// 			// keep first m_size individuals of shuffled indexes
// 			nCase = std::min(m_numCases[sp], static_cast<int>(m_caseIdx[sp].size()));
// 			for (int i = 0; i < nCase; ++i)
// 				pop.ind(m_caseIdx[sp][i], sp).setSubPopID(0);
// 			// remove others
// 			for (int i = nCase, iEnd = m_caseIdx[sp].size(); i < iEnd; ++i)
// 				pop.ind(m_caseIdx[sp][i], sp).setSubPopID(-1);
// 
// 			// keep first m_size individuals of shuffled indexes
// 			nControl = std::min(m_numControls[sp], static_cast<int>(m_controlIdx[sp].size()));
// 			for (int i = 0; i < nControl; ++i)
// 				pop.ind(m_controlIdx[sp][i], sp).setSubPopID(1);
// 			// remove others
// 			for (int i = nControl, iEnd = m_controlIdx[sp].size(); i < iEnd; ++i)
// 				pop.ind(m_controlIdx[sp][i], sp).setSubPopID(-1);
// 		}
// 		// newPop .... but ignore ancestral populations
// 		population & sample = pop.newPopByIndID(0);
// 		sample.setIntVectorVar("nCases", m_numCases);
// 		sample.setIntVectorVar("nControls", m_numControls);
// 		return sample;
// 	}
// }
// 
// 
// bool affectedSibpairSample::prepareSample(population & pop)
// {
// 	// get parental index for each subpop
// 	DBG_FAILIF(m_size.size() > 1 && m_size.size() != pop.numSubPop(),
// 		ValueError,
// 		"Length of array size and number of subpopulations do not match.");
// 
// 	m_validSibs.clear();
// 
// 	m_father_id = pop.infoIdx(infoField(0));
// 	m_mother_id = pop.infoIdx(infoField(1));
// 	vectorstr fields(5);
// 	fields[0] = "oldindex";
// 	fields[1] = "pedindex";
// 	fields[2] = "offspring0";
// 	fields[3] = "offspring1";
// 	fields[4] = "spouse";
// 	pop.addInfoFields(fields, -1);
// 	pop.locateRelatives(REL_Self, vectorstr(1, "oldindex"));
// 	UINT pedindexIdx = pop.infoIdx("pedindex");
// 	UINT off0Idx = pop.infoIdx("offspring0");
// 	UINT off1Idx = pop.infoIdx("offspring1");
// 	UINT spouseIdx = pop.infoIdx("spouse");
// 
// 	pop.locateRelatives(REL_Offspring, vectorstr(fields.begin() + 2, fields.begin() + 4));
// 	pop.locateRelatives(REL_Spouse, vectorstr(1, "spouse"));
// 	//
// 	// find sibpairs from the parental generation.
// 	pop.useAncestralGen(1);
// 	vectorlu off;
// 	int pedIdx = 0;
// 	// valid sibpairs for each subpopulation
// 	for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
// 		for (int it = 0; it < static_cast<int>(pop.subPopSize(sp)); ++it) {
// 			individual & ind = pop.ind(it, sp);
// 			// individual already belongs to another family
// 			if (ind.info(pedindexIdx) != -1.)
// 				continue;
// 			int spouse = ind.intInfo(spouseIdx);
// 			// if no spouse, or spouse belongs to anther pedigree, or there are less than two kids
// 			if (spouse == -1 || pop.ind(spouse).info(pedindexIdx) >= 0
// 			    || pop.ind(spouse).info(spouseIdx) != it
// 			    || ind.info(off1Idx) == -1.)
// 				continue;
// 			int child0 = ind.intInfo(off0Idx);
// 			int child1 = ind.intInfo(off1Idx);
// 			int child0_f = pop.ancestor(child0, 0).intInfo(m_father_id);
// 			int child0_m = pop.ancestor(child0, 0).intInfo(m_mother_id);
// 			int child1_f = pop.ancestor(child1, 0).intInfo(m_father_id);
// 			int child1_m = pop.ancestor(child1, 0).intInfo(m_mother_id);
// 			// invalid child
// 			if ((child0_f != it && child0_f != spouse) ||
// 			    (child0_m != it && child0_m != spouse) ||
// 			    (child1_f != it && child1_f != spouse) ||
// 			    (child1_m != it && child1_m != spouse))
// 				continue;
// 			ind.setInfo(pedIdx, pedindexIdx);
// 			pop.ind(spouse).setInfo(pedIdx, pedindexIdx);
// 			off.push_back(ind.intInfo(off0Idx));
// 			off.push_back(ind.intInfo(off1Idx));
// 			pedIdx++;
// 		}
// 		DBG_DO(DBG_SELECTOR, cout << "Number of sibpairs in subpop " << sp << " is "
// 			                      << m_validSibs[sp].size() << endl);
// 	}                                                                                         // each subpop
// 	pop.useAncestralGen(0);
// 	m_validSibs.resize(pop.numSubPop());
// 	for (UINT sp = 0; sp < pop.numSubPop(); ++sp)
// 		m_validSibs[sp].clear();
// 	UINT nSibs = 0;
// 	for (size_t i = 0; i < off.size() / 2; ++i) {
// 		pop.ind(off[2 * i]).setInfo(i, pedindexIdx);
// 		pop.ind(off[2 * i + 1]).setInfo(i, pedindexIdx);
// 		if (m_affectedness == pop.ind(off[2 * i]).affected()
// 		    && m_affectedness == pop.ind(off[2 * i + 1]).affected()) {
// 			m_validSibs[pop.subPopIndPair(off[2 * i]).first].push_back(i);
// 			nSibs++;
// 		}
// 	}
// 	pop.setIntVar("numAffectedSibpairs", nSibs);
// 	// do not do sampling if countOnly
// 	if (m_countOnly)
// 		return false;
// 	else
// 		return true;
// }
// 
// 
// population & affectedSibpairSample::drawsample(population & pop)
// {
// 	// mark to remove everyone
// 	pop.setIndSubPopID(vectori(1, -1));
// 	vectori acceptedSibs;
// 
// 	if (m_size.size() <= 1) {                                         // draw from the whole population
// 		// collect all families
// 		vector<UINT> allSibs;
// 
// 		for (UINT sp = 0; sp < pop.numSubPop(); ++sp)
// 			allSibs.insert(allSibs.end(), m_validSibs[sp].begin(), m_validSibs[sp].end());
// 
// 		if (!m_size.empty() && m_size[0] > allSibs.size())
// 			cout << "Warning: Not enough sibpairs to be sampled. Requested "
// 			     << m_size[0] << ", existing " << allSibs.size() << endl;
// 
// 		UINT N = 0;
// 		if (m_size.empty() || m_size[0] > allSibs.size())
// 			N = allSibs.size();
// 		else {
// 			N = m_size[0];
// 			// random shuffle.
// 			random_shuffle(allSibs.begin(), allSibs.end());
// 		}
// 		acceptedSibs.insert(acceptedSibs.end(), allSibs.begin(), allSibs.begin() + N);
// 	} else {                                                                          // for each subpop
// 		for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
// 			vectorlu & sibpairs = m_validSibs[sp];
// 
// 			UINT N = sibpairs.size();
// 			if (N > m_size[sp]) {
// 				N = m_size[sp];
// 				// sample sibpairs
// 				random_shuffle(sibpairs.begin(), sibpairs.end());
// 			}
// 			acceptedSibs.insert(acceptedSibs.end(), sibpairs.begin(), sibpairs.begin() + N);
// 		}                                                                                 // sp
// 	}
// 	// now, we have acepted sibs, set subpopid in preparation for a new
// 	// population
// 	pop.useAncestralGen(1);
// 
// 	vectorlu offspring;
// 	int pedIdx = 0;
// 	vectorlu off;
// 	UINT spouseIdx = pop.infoIdx("spouse");
// 	UINT pedindexIdx = pop.infoIdx("pedindex");
// 	UINT off0Idx = pop.infoIdx("offspring0");
// 	UINT off1Idx = pop.infoIdx("offspring1");
// 
// 	for (size_t i = 0; i < pop.popSize(); ++i) {
// 		individual & ind = pop.ind(i);
// 		int infoPedIdx = ind.intInfo(pedindexIdx);
// 		if (infoPedIdx == -1)
// 			continue;
// 		int spouse = ind.intInfo(spouseIdx);
// 		// only look forward
// 		if (spouse < static_cast<int>(i))
// 			continue;
// 		// if this family is selected
// 		if (find(acceptedSibs.begin(), acceptedSibs.end(), infoPedIdx) != acceptedSibs.end()) {
// 			ind.setSubPopID(pedIdx);
// 			pop.ind(spouse).setSubPopID(pedIdx);
// 			off.push_back(ind.intInfo(off0Idx));
// 			off.push_back(ind.intInfo(off1Idx));
// 			pedIdx++;
// 		}
// 	}
// 	pop.useAncestralGen(0);
// 	for (size_t i = 0; i < off.size() / 2; ++i) {
// 		pop.ind(off[2 * i]).setSubPopID(i);
// 		pop.ind(off[2 * i + 1]).setSubPopID(i);
// 	}
// 
// 	// this is offspring population with copy of ancestral pop
// 	// (1 means keep only one ancestral population
// 	population & newPop = pop.newPopByIndID(1);
// 	//
// 	resetParentalIndex(newPop, "father_idx", "mother_idx", "oldindex");
// 
// 	DBG_DO(DBG_SELECTOR, cout << "Offspring selection done." << endl);
// 	return newPop;
// }
// 
// 
// bool largePedigreeSample::prepareSample(population & pop)
// {
// 	DBG_FAILIF(pop.ancestralGens() < 2, ValueError,
// 		"At least two ancestral populations are needed to draw large pedigrees");
// 	//
// 	m_validPedigrees.clear();
// 	UINT nPed = 0;
// 	//
// 	vectorstr fields(3 + m_maxOffspring);
// 	vectori offspringIdx(m_maxOffspring);
// 	for (size_t i = 0; i < m_maxOffspring; ++i)
// 		fields[i] = "offspring" + toStr(i);
// 	fields[m_maxOffspring] = "pedindex";
// 	fields[m_maxOffspring + 1] = "spouse";
// 	fields[m_maxOffspring + 2] = "oldindex";
// 	// add info fields
// 	pop.addInfoFields(fields, -1);
// 	UINT pedindexIdx = pop.infoIdx("pedindex");
// 	UINT spouseIdx = pop.infoIdx("spouse");
// 	UINT fatherIdx = pop.infoIdx("father_idx");
// 	UINT motherIdx = pop.infoIdx("mother_idx");
// 	for (size_t i = 0; i < m_maxOffspring; ++i)
// 		offspringIdx[i] = pop.infoIdx(fields[i]);
// 	// save old index
// 	pop.locateRelatives(REL_Self, vectorstr(1, "oldindex"));
// 	pop.locateRelatives(REL_Offspring, vectorstr(fields.begin(), fields.begin() + m_maxOffspring));
// 	pop.locateRelatives(REL_Spouse, vectorstr(1, "spouse"));
// 
// 	//
// 	pop.useAncestralGen(2);
// 	m_validPedigrees.resize(pop.numSubPop());
// 	size_t pedindex = 0;
// 	DBG_DO(DBG_SELECTOR, cout << "Finding all three-generation pedigrees" << endl);
// 	for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
// 		m_validPedigrees[sp].clear();
// 		int g3start = pop.subPopBegin(sp);
// 		int g3end = pop.subPopEnd(sp);
// 		//
// 		for (int idx = g3start; idx < g3end; ++idx) {
// 			pop.useAncestralGen(2);
// 			unsigned pedSize = 2;
// 			unsigned numAffected = 0;
// 			//
// 			// already belong to other pedigree
// 			int grandspouse = pop.ind(idx).intInfo(spouseIdx);
// 			// no spuse? one of grandparents belong to another pedigree?
// 			if (grandspouse < 0 || pop.ind(idx).info(pedindexIdx) >= 1 ||
// 			    pop.ind(grandspouse).intInfo(spouseIdx) != idx || // spouse not paired
// 			    pop.ind(grandspouse).info(pedindexIdx) >= 1)
// 				continue;
// 			if (pop.ind(idx).affected())
// 				numAffected++;
// 			if (pop.ind(grandspouse).affected())
// 				numAffected++;
// 			//
// 			vectorlu parentsTmp;
// 			for (size_t x = 0; x < m_maxOffspring; ++x) {
// 				int off = pop.ind(idx).intInfo(offspringIdx[x]);
// 				if (off != -1)
// 					parentsTmp.push_back(off);
// 			}
// 			// parents generation
// 			//
// 			if (parentsTmp.empty())
// 				continue;
// 			// go to parent generation
// 			pop.useAncestralGen(1);
// 			vectorlu spouseofparents;
// 			vectorlu childrenTmp;
// 			// verify parents
// 			vectorlu parents;
// 			for (vectorlu::iterator par = parentsTmp.begin(); par != parentsTmp.end(); ++par) {
// 				if (pop.ind(*par).info(pedindexIdx) == -1.) {
// 					int ind_p = pop.ind(*par).intInfo(fatherIdx);
// 					int ind_m = pop.ind(*par).intInfo(motherIdx);
// 					// child with another spouse
// 					if (grandspouse != ind_p && grandspouse != ind_m)
// 						continue;
// 					if (idx != ind_p && idx != ind_m)
// 						continue;
// 
// 					parents.push_back(*par);
// 					pedSize++;
// 					if (pop.ind(*par).affected())
// 						numAffected++;
// 				}
// 			}
// 			//
// 			for (vectorlu::iterator par = parents.begin(); par != parents.end(); ++par) {
// 				int spouse = pop.ind(*par).intInfo(spouseIdx);
// 				// if there is spouse, add it in
// 				if (spouse >= 0 && pop.ind(*par).info(pedindexIdx) == -1.) {
// 					// if spouse relationship is mutual
// 					if (pop.ind(spouse).info(spouseIdx) != *par)
// 						continue;
// 					spouseofparents.push_back(spouse);
// 					pedSize++;
// 					if (pop.ind(spouse).affected())
// 						numAffected++;
// 					// there are children only when there is spouse
// 					for (size_t x = 0; x < m_maxOffspring; ++x) {
// 						int off = pop.ind(*par).intInfo(offspringIdx[x]);
// 						if (off != -1)
// 							childrenTmp.push_back(off);
// 					}
// 				}
// 			}
// 			if (childrenTmp.empty())
// 				continue;
// 			pop.useAncestralGen(0);
// 			// verify children
// 			vectorlu children;
// 			for (vectorlu::iterator child = childrenTmp.begin(); child != childrenTmp.end(); ++child) {
// 				// unoccupied.
// 				if (pop.ind(*child).info(pedindexIdx) == -1.) {
// 					ULONG ind_p = pop.ind(*child).intInfo(fatherIdx);
// 					ULONG ind_m = pop.ind(*child).intInfo(motherIdx);
// 					//
// 					if (find(parents.begin(), parents.end(), ind_m) == parents.end() &&
// 					    find(spouseofparents.begin(), spouseofparents.end(), ind_m) == spouseofparents.end())
// 						continue;
// 					if (find(parents.begin(), parents.end(), ind_p) == parents.end() &&
// 					    find(spouseofparents.begin(), spouseofparents.end(), ind_p) == spouseofparents.end())
// 						continue;
// 					children.push_back(*child);
// 					pedSize++;
// 					if (pop.ind(*child).affected())
// 						numAffected++;
// 				}
// 			}
// 			if (children.empty())
// 				continue;
// 			if (pedSize < m_minPedSize || (m_minAffected > 0 && numAffected < m_minAffected))
// 				continue;
// 			// everything seems to be fine, add them to the family
// 			pop.useAncestralGen(2);
// 			// grandparents
// 			pop.ind(idx).setInfo(pedindex, pedindexIdx);
// 			pop.ind(grandspouse).setInfo(pedindex, pedindexIdx);
// 			// parents and their spouse
// 			pop.useAncestralGen(1);
// 			for (vectorlu::iterator it = parents.begin(); it != parents.end(); ++it)
// 				pop.ind(*it).setInfo(pedindex, pedindexIdx);
// 			// spouse of parents
// 			for (vectorlu::iterator it = spouseofparents.begin(); it != spouseofparents.end(); ++it)
// 				pop.ind(*it).setInfo(pedindex, pedindexIdx);
// 			// now children
// 			pop.useAncestralGen(0);
// 			for (vectorlu::iterator it = children.begin(); it != children.end(); ++it)
// 				pop.ind(*it).setInfo(pedindex, pedindexIdx);
// 			// is this family qualified?
// 			m_validPedigrees[sp].push_back(boost::tie(pedindex, pedSize));
// 			pedindex += 1;
// 		}
// 		nPed += m_validPedigrees[sp].size();
// 	}
// 	pop.useAncestralGen(0);
// 	pop.setIntVar("numPedigrees", nPed);
// 
// 	// do not do sampling if countOnly
// 	if (m_countOnly)
// 		return false;
// 	else
// 		return true;
// }
// 
// 
// population & largePedigreeSample::drawsample(population & pop)
// {
// 	// sample sibpairs
// 	DBG_DO(DBG_SELECTOR, cout << "Generating sample" << endl);
// 	pop.setIndSubPopID(vectori(1, -1));
// 	pedArray acceptedPeds;
// 
// 	if (m_size.size() <= 1) {                                         // draw from the whole population
// 		// collect all families
// 		pedArray allPeds;
// 
// 		for (UINT sp = 0; sp < pop.numSubPop(); ++sp)
// 			allPeds.insert(allPeds.end(), m_validPedigrees[sp].begin(), m_validPedigrees[sp].end());
// 
// 		if (!m_size.empty() && m_size[0] > allPeds.size())
// 			cout << "Warning: Not enough sibpairs to be sampled. Requested "
// 			     << m_size[0] << ", existing " << allPeds.size() << endl;
// 
// 		random_shuffle(allPeds.begin(), allPeds.end());
// 		UINT N = 0;
// 		// consider total size.
// 		if (m_size.empty() || m_minTotalSize > 0) {
// 			size_t totalSize = 0;
// 			for (pedArray::iterator it = allPeds.begin(); it != allPeds.end(); ++it) {
// 				totalSize += boost::get<1>(*it);
// 				N++;
// 				if (totalSize > m_minTotalSize)
// 					break;
// 			}
// 			if (totalSize < m_minTotalSize)
// 				cout << "Warning: can not reach min total size" << endl;
// 		} else
// 			N = m_size[0];
// 		acceptedPeds.insert(acceptedPeds.end(), allPeds.begin(), allPeds.begin() + N);
// 	} else {                                                                          // for each subpop
// 		for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
// 			pedArray & peds = m_validPedigrees[sp];
// 
// 			UINT N = peds.size();
// 			if (N > m_size[sp]) {
// 				N = m_size[sp];
// 				// sample sibpairs
// 				random_shuffle(peds.begin(), peds.end());
// 			}
// 			acceptedPeds.insert(acceptedPeds.end(), peds.begin(), peds.begin() + N);
// 		}                                                                                 // sp
// 	}
// 
// 	DBG_DO(DBG_SELECTOR, cout << "Sampling " << acceptedPeds.size() << " pedigrees " << endl);
// 
// 	//
// 	UINT spouseIdx = pop.infoIdx("spouse");
// 	UINT pedindexIdx = pop.infoIdx("pedindex");
// 	vectori offspringIdx(m_maxOffspring);
// 	for (size_t i = 0; i < m_maxOffspring; ++i)
// 		offspringIdx[i] = pop.infoIdx("offspring" + toStr(i));
// 	pop.useAncestralGen(2);
// 	vectorf grandIdx = pop.indInfo(pedindexIdx);
// 	int newPedID = 0;
// 	for (pedArray::iterator ped = acceptedPeds.begin(); ped != acceptedPeds.end(); ++ped, ++newPedID) {
// 		double pedID = boost::get<0>(*ped);
// 		DBG_DO(DBG_SELECTOR, cout << "Getting pedigree " << pedID << " of size " << boost::get<1>(*ped) << endl);
// 		// real pedsize, for verification purpose
// 		int ps = 0;
// 		// find these offspring....
// 		pop.useAncestralGen(2);
// 		// find grandparent
// 		vectorf::iterator tmp = find(grandIdx.begin(), grandIdx.end(), pedID);
// 		DBG_FAILIF(tmp == grandIdx.end(), ValueError, "Can not find pedigree");
// 		//
// 		int grandpar1 = tmp - grandIdx.begin();
// 		DBG_FAILIF(pop.ind(grandpar1).info(spouseIdx) == -1., SystemError,
// 			"Grand parent's spouse is invalid");
// 		int grandpar2 = pop.ind(grandpar1).intInfo(spouseIdx);
// 		DBG_FAILIF(pop.ind(grandpar2).info(pedindexIdx) != pedID, SystemError,
// 			"Grand parent's spouse is invalid");
// 		pop.ind(grandpar1).setSubPopID(newPedID);
// 		pop.ind(grandpar2).setSubPopID(newPedID);
// 		ps += 2;
// 		// find parents
// 		vectori parents;
// 		for (size_t x = 0; x < m_maxOffspring; ++x) {
// 			int off = pop.ind(grandpar1).intInfo(offspringIdx[x]);
// 			if (off != -1)
// 				parents.push_back(off);
// 		}
// 		//
// 		pop.useAncestralGen(1);
// 		vectori children;
// 		for (vectori::iterator it = parents.begin(); it != parents.end(); ++it) {
// 			if (pop.ind(*it).info(pedindexIdx) != pedID)
// 				continue;
// 			pop.ind(*it).setSubPopID(newPedID);
// 			ps++;
// 			int spouse = pop.ind(*it).intInfo(spouseIdx);
// 			if (spouse < 0 || pop.ind(spouse).info(pedindexIdx) != pedID)
// 				continue;
// 			// if there is spouse, add it in
// 			pop.ind(spouse).setSubPopID(newPedID);
// 			ps++;
// 			for (size_t x = 0; x < m_maxOffspring; ++x) {
// 				int off = pop.ind(*it).intInfo(offspringIdx[x]);
// 				if (off != -1)
// 					children.push_back(off);
// 			}
// 		}
// 		// go to children
// 		pop.useAncestralGen(0);
// 		for (vectori::iterator it = children.begin(); it != children.end(); ++it) {
// 			if (pop.ind(*it).info(pedindexIdx) == pedID) {
// 				pop.ind(*it).setSubPopID(newPedID);
// 				ps++;
// 			}
// 		}
// 		DBG_FAILIF(ps != boost::get<1>(*ped), SystemError,
// 			"Pedigree sizes do not match, estimated " + toStr(boost::get<1>(*ped)) + " real: " + toStr(ps));
// 	}
// 	// just to make sure
// 	pop.useAncestralGen(0);
// 	// saving samples in a new population,
// 	population & newPop = pop.newPopByIndID(2);
// 	//
// 	resetParentalIndex(newPop, "father_idx", "mother_idx", "oldindex");
// 	return newPop;
// }
// 
// 
// bool nuclearFamilySample::prepareSample(population & pop)
// {
// 	DBG_FAILIF(pop.ancestralGens() < 1, ValueError,
// 		"At least one ancestral populations are needed to draw large pedigrees");
// 	//
// 	m_validPedigrees.clear();
// 	UINT nPed = 0;
// 	//
// 	vectorstr fields(3 + m_maxOffspring);
// 	vectori offspringIdx(m_maxOffspring);
// 	for (size_t i = 0; i < m_maxOffspring; ++i)
// 		fields[i] = "offspring" + toStr(i);
// 	fields[m_maxOffspring] = "pedindex";
// 	fields[m_maxOffspring + 1] = "spouse";
// 	fields[m_maxOffspring + 2] = "oldindex";
// 	// add info fields
// 	pop.addInfoFields(fields, -1);
// 	UINT pedindexIdx = pop.infoIdx("pedindex");
// 	UINT spouseIdx = pop.infoIdx("spouse");
// 	for (size_t i = 0; i < m_maxOffspring; ++i)
// 		offspringIdx[i] = pop.infoIdx(fields[i]);
// 	// save old index
// 	pop.locateRelatives(REL_Self, vectorstr(1, "oldindex"));
// 	pop.locateRelatives(REL_Offspring, vectorstr(fields.begin(), fields.begin() + m_maxOffspring));
// 	pop.locateRelatives(REL_Spouse, vectorstr(1, "spouse"));
// 	// offspring index
// 	m_validPedigrees.resize(pop.numSubPop());
// 	size_t pedIdx = 0;
// 	DBG_DO(DBG_SELECTOR, cout << "Finding all two-generation pedigrees" << endl);
// 	for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
// 		pop.useAncestralGen(1);
// 		vectori off;
// 		for (IndIterator it = pop.indBegin(sp); it.valid(); ++it) {
// 			// individual already belongs to another family
// 			if (it->info(pedindexIdx) != -1.)
// 				continue;
// 			int spouse = it->intInfo(spouseIdx);
// 			// has spouse, spouse does not belong to anther ped
// 			if (spouse != -1. && pop.ind(spouse).info(pedindexIdx) == -1.) {
// 				it->setInfo(pedIdx, pedindexIdx);
// 				pop.ind(spouse).setInfo(pedIdx, pedindexIdx);
// 				// many of the offspring may be -1.
// 				for (UINT oi = 0; oi < m_maxOffspring; ++oi)
// 					off.push_back(it->intInfo(offspringIdx[oi]));
// 				pedIdx++;
// 			}
// 		}
// 		// go to offspring generation and verify
// 		pop.useAncestralGen(0);
// 		m_validPedigrees[sp].clear();
// 		for (UINT i = 0; i < off.size() / m_maxOffspring; ++i) {
// 			UINT pedSize = 2;
// 			UINT pedAffected = 0;
// 			for (UINT oi = 0; oi < m_maxOffspring; ++oi) {
// 				// valid offspring
// 				if (off[i * m_maxOffspring + oi] != -1) {
// 					pop.ind(off[i * m_maxOffspring + oi]).setInfo(i, pedindexIdx);
// 					pedSize++;
// 					if (pop.ind(off[i * m_maxOffspring + oi]).affected())
// 						pedAffected++;
// 				}
// 			}
// 			if (pedSize >= m_minPedSize && pedAffected >= m_minAffected)
// 				m_validPedigrees[sp].push_back(boost::tie(i, pedSize));
// 		}
// 		nPed += m_validPedigrees[sp].size();
// 		DBG_DO(DBG_SELECTOR, cout << "Number of valid pedigrees in subpop " + toStr(sp) + " is "
// 			+ toStr(m_validPedigrees[sp].size()) << endl);
// 	}
// 	pop.useAncestralGen(0);
// 	pop.setIntVar("numPedigrees", nPed);
// 
// 	// do not do sampling if countOnly
// 	if (m_countOnly)
// 		return false;
// 	else
// 		return true;
// }
// 
// 
// population & nuclearFamilySample::drawsample(population & pop)
// {
// 	DBG_DO(DBG_SELECTOR, cout << "Generating nuclear family" << endl);
// 	pop.setIndSubPopID(vectori(1, -1));
// 	pedArray acceptedPeds;
// 
// 	if (m_size.size() <= 1) {                                         // draw from the whole population
// 		// collect all families
// 		pedArray allPeds;
// 
// 		for (UINT sp = 0; sp < pop.numSubPop(); ++sp)
// 			allPeds.insert(allPeds.end(), m_validPedigrees[sp].begin(), m_validPedigrees[sp].end());
// 
// 		if (!m_size.empty() && m_size[0] > allPeds.size())
// 			cout << "Warning: Not enough nuclear families to be sampled. Requested "
// 			     << m_size[0] << ", existing " << allPeds.size() << endl;
// 
// 		random_shuffle(allPeds.begin(), allPeds.end());
// 		UINT N = 0;
// 		// consider total size.
// 		if (m_size.empty() || m_minTotalSize > 0) {
// 			size_t totalSize = 0;
// 			for (pedArray::iterator it = allPeds.begin(); it != allPeds.end(); ++it) {
// 				totalSize += boost::get<1>(*it);
// 				N++;
// 				if (totalSize > m_minTotalSize)
// 					break;
// 			}
// 			if (totalSize < m_minTotalSize)
// 				cout << "Warning: can not reach min total size" << endl;
// 		} else
// 			N = m_size[0];
// 		acceptedPeds.insert(acceptedPeds.end(), allPeds.begin(), allPeds.begin() + N);
// 	} else {                                                                          // for each subpop
// 		for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
// 			pedArray & peds = m_validPedigrees[sp];
// 
// 			UINT N = peds.size();
// 			if (N > m_size[sp]) {
// 				N = m_size[sp];
// 				random_shuffle(peds.begin(), peds.end());
// 			}
// 			acceptedPeds.insert(acceptedPeds.end(), peds.begin(), peds.begin() + N);
// 		}                                                                                 // sp
// 	}
// 
// 	DBG_DO(DBG_SELECTOR, cout << "Sampling " << acceptedPeds.size() << " pedigrees " << endl);
// 
// 	//
// 	UINT spouseIdx = pop.infoIdx("spouse");
// 	UINT pedindexIdx = pop.infoIdx("pedindex");
// 	vectori offspringIdx(m_maxOffspring);
// 	for (size_t i = 0; i < m_maxOffspring; ++i)
// 		offspringIdx[i] = pop.infoIdx("offspring" + toStr(i));
// 	pop.useAncestralGen(1);
// 	vectori off;
// 	int pedIdx = 0;
// 
// 	for (size_t i = 0; i < pop.popSize(); ++i) {
// 		individual & ind = pop.ind(i);
// 		int infoPedIdx = ind.intInfo(pedindexIdx);
// 		if (infoPedIdx == -1)
// 			continue;
// 		int spouse = ind.intInfo(spouseIdx);
// 		// only look forward
// 		if (spouse < static_cast<int>(i))
// 			continue;
// 		// if this family is selected
// 		for (pedArray::iterator it = acceptedPeds.begin(); it != acceptedPeds.end(); ++it) {
// 			if (boost::get<0>(*it) == infoPedIdx) {
// 				ind.setSubPopID(pedIdx);
// 				pop.ind(spouse).setSubPopID(pedIdx);
// 				for (UINT oi = 0; oi < m_maxOffspring; ++oi)
// 					off.push_back(ind.intInfo(offspringIdx[oi]));
// 				pedIdx++;
// 				break;
// 			}
// 		}
// 	}
// 	pop.useAncestralGen(0);
// 	for (size_t i = 0; i < off.size(); ++i) {
// 		if (off[i] != -1)
// 			pop.ind(off[i]).setSubPopID(i / m_maxOffspring);
// 	}
// 	population & newPop = pop.newPopByIndID(1);
// 	//
// 	resetParentalIndex(newPop, "father_idx", "mother_idx", "oldindex");
// 	return newPop;
// }
// 
 
}
