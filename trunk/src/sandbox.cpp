
/**
 *  $File: sandbox.cpp $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
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

#include "sandbox.h"

namespace simuPOP {
namespace sandbox {

#ifdef LONGALLELE

bool RevertFixedSites::apply(Population & pop) const
{
	if (pop.popSize() == 0 || pop.totNumLoci() == 0)
		return true;

	bool chX = pop.chromType(0) == CHROMOSOME_X;

	RawIndIterator it = pop.rawIndBegin();
	RawIndIterator it_end = pop.rawIndEnd();
	std::set<Allele> commonAlleles(it->genoBegin(0), it->genoEnd(0));
	commonAlleles.erase(0);
	if (commonAlleles.size() == 0)
		return true;

	for (; it != it_end; ++it) {
		// common = commonAlleles & geno0
		std::set<Allele> common;
		std::set<Allele> alleles1(it->genoBegin(0), it->genoEnd(0));
		set_intersection(commonAlleles.begin(),
			commonAlleles.end(), alleles1.begin(), alleles1.end(),
			std::inserter(common, common.begin()));
		// commonAlleles = common & geno1
		if (chX && it->sex() == MALE) {
			// commonAlleles = common
			commonAlleles.swap(common);
		} else {
			commonAlleles.clear();
			std::set<Allele> alleles2(it->genoBegin(1), it->genoEnd(1));
			set_intersection(common.begin(),
				common.end(), alleles2.begin(), alleles2.end(),
				std::inserter(commonAlleles, commonAlleles.begin()));
		}
		if (commonAlleles.size() == 0)
			return true;
	}
	if (!noOutput()) {
		ostream & out = getOstream(pop.dict());
		out << pop.gen();
		std::set<Allele>::iterator beg = commonAlleles.begin();
		std::set<Allele>::iterator end = commonAlleles.end();
		for (; beg != end ; ++beg)
			out << '\t' << *beg;
		out << endl;
	}
	it = pop.rawIndBegin();
	vectora new_alleles(pop.totNumLoci());
	for (; it != it_end; ++it) {
		for (size_t p = 0; p < 2; ++p) {
			if (p == 1 && chX && it->sex() == MALE)
				continue;
			std::set<Allele> old_alleles(it->genoBegin(p), it->genoEnd(p));
			old_alleles.erase(0);
			std::fill(new_alleles.begin(), new_alleles.end(), Allele(0));
			set_difference(old_alleles.begin(), old_alleles.end(),
				commonAlleles.begin(), commonAlleles.end(), new_alleles.begin());
			std::copy(new_alleles.begin(), new_alleles.end(),
				it->genoBegin(p));
		}
	}
	return true;
}


double MutSpaceSelector::indFitness(Population & /* pop */, RawIndIterator ind) const
{
	if (ind->sex() == MALE && ind->chromType(0) == CHROMOSOME_X) {
		if (m_mode == MULTIPLICATIVE) {
			return randomSelMulFitnessExt(ind->genoBegin(0), ind->genoEnd(0), true);
		} else if (m_mode == ADDITIVE) {
			if (m_additive)
				return randomSelAddFitness(ind->genoBegin(0), ind->genoEnd(0), true);
			else
				return randomSelAddFitnessExt(ind->genoBegin(0), ind->genoEnd(0), true);
		} else if (m_mode == EXPONENTIAL) {
			if (m_additive)
				return randomSelExpFitness(ind->genoBegin(0), ind->genoEnd(0), true);
			else
				return randomSelExpFitnessExt(ind->genoBegin(0), ind->genoEnd(0), true);
		}
	} else {
		if (m_mode == MULTIPLICATIVE) {
			return randomSelMulFitnessExt(ind->genoBegin(), ind->genoEnd(), false);
		} else if (m_mode == ADDITIVE) {
			if (m_additive)
				return randomSelAddFitness(ind->genoBegin(), ind->genoEnd(), false);
			else
				return randomSelAddFitnessExt(ind->genoBegin(), ind->genoEnd(), false);
		} else if (m_mode == EXPONENTIAL) {
			if (m_additive)
				return randomSelExpFitness(ind->genoBegin(), ind->genoEnd(), false);
			else
				return randomSelExpFitnessExt(ind->genoBegin(), ind->genoEnd(), false);
		}
	}
	return 0;
}


bool MutSpaceSelector::apply(Population & pop) const
{
	m_newMutants.clear();
	if (!BaseSelector::apply(pop))
		return false;
	// output NEW mutant...
	if (!m_newMutants.empty() && !noOutput()) {
		ostream & out = getOstream(pop.dict());
		vectoru::const_iterator it = m_newMutants.begin();
		vectoru::const_iterator it_end = m_newMutants.end();
		for (; it != it_end; ++it) {
			SelCoef s = m_selFactory[*it];
			out << *it << '\t' << s.first << '\t' << s.second << '\n';
		}
		closeOstream();
	}
	return true;
}


MutSpaceSelector::SelCoef MutSpaceSelector::getFitnessValue(size_t mutant) const
{
	size_t sz = m_selDist.size();
	double s = 0;
	double h = 0.5;

	if (sz == 0) {
		// call a function
		const pyFunc & func = m_selDist.func();
		PyObject * res;
		if (func.numArgs() == 0)
			res = func("()");
		else {
			DBG_FAILIF(func.arg(0) != "loc", ValueError,
				"Only parameter loc is accepted for this user-defined function.");
			res = func("(i)", mutant);
		}
		if (PyNumber_Check(res)) {
			s = PyFloat_AsDouble(res);
		} else if (PySequence_Check(res)) {
			size_t sz = PySequence_Size(res);
			DBG_FAILIF(sz == 0, RuntimeError, "Function return an empty list.");
			PyObject * item = PySequence_GetItem(res, 0);
			s = PyFloat_AsDouble(item);
			Py_DECREF(item);
			if (sz > 1) {
				item = PySequence_GetItem(res, 1);
				h = PyFloat_AsDouble(item);
				Py_DECREF(item);
			}
		}
		Py_DECREF(res);
	} else {
		int mode = static_cast<int>(m_selDist[0]);
		if (mode == CONSTANT) {
			// constant
			s = m_selDist[1];
			if (m_selDist.size() > 2)
				h = m_selDist[2];
		} else {
			// a gamma distribution
			s = getRNG().randGamma(m_selDist[1], m_selDist[2]);
			if (m_selDist.size() > 3)
				h = m_selDist[3];
		}
	}
	m_selFactory[mutant] = SelCoef(s, h);
	m_newMutants.push_back(mutant);
	if (m_additive && h != 0.5)
		m_additive = false;
	return SelCoef(s, h);
}


double MutSpaceSelector::randomSelAddFitness(GenoIterator it, GenoIterator it_end, bool chrX) const
{
	double s = 0;

	for (; it != it_end; ++it) {
		if (*it == 0u)
			continue;
		SelMap::iterator sit = m_selFactory.find(static_cast<unsigned int>(*it));
		if (sit == m_selFactory.end())
			s += getFitnessValue(*it).first / 2.;
		else
			s += sit->second.first / 2;
	}
	if (chrX)
		// fitness of variant on chromosome X is as if it is homogeneous
		s += s;
	return 1 - s > 0 ? 1 - s : 0;
}


double MutSpaceSelector::randomSelExpFitness(GenoIterator it, GenoIterator it_end, bool chrX) const
{
	double s = 0;

	for (; it != it_end; ++it) {
		if (*it == 0u)
			continue;
		SelMap::iterator sit = m_selFactory.find(static_cast<unsigned int>(*it));
		if (sit == m_selFactory.end())
			s += getFitnessValue(*it).first / 2.;
		else
			s += sit->second.first / 2;
	}
	if (chrX)
		// fitness of variant on chromosome X is as if it is homogeneous
		s += s;
	return exp(-s);
}


double MutSpaceSelector::randomSelMulFitnessExt(GenoIterator it, GenoIterator it_end, bool chrX) const
{
	MutCounter cnt;

	for (; it != it_end; ++it) {
		if (*it == 0u)
			continue;
		MutCounter::iterator mit = cnt.find(*it);
		if (mit == cnt.end())
			cnt[*it] = 1;
		else
			++mit->second;
	}

	double s = 1;
	MutCounter::iterator mit = cnt.begin();
	MutCounter::iterator mit_end = cnt.end();
	for (; mit != mit_end; ++mit) {
		SelMap::iterator sit = m_selFactory.find(mit->first);
		if (sit == m_selFactory.end()) {
			SelCoef sf = getFitnessValue(mit->first);
			if (mit->second == 1 && !chrX)
				s *= 1 - sf.first * sf.second;
			else
				s *= 1 - sf.first;
		} else {
			if (mit->second == 1 && !chrX)
				s *= 1 - sit->second.first * sit->second.second;
			else
				s *= 1 - sit->second.first;
		}
	}
	return s;
}


double MutSpaceSelector::randomSelAddFitnessExt(GenoIterator it, GenoIterator it_end, bool chrX) const
{
	MutCounter cnt;

	for (; it != it_end; ++it) {
		if (*it == 0u)
			continue;
		MutCounter::iterator mit = cnt.find(*it);
		if (mit == cnt.end())
			cnt[*it] = 1;
		else
			++mit->second;
	}

	double s = 0;
	MutCounter::iterator mit = cnt.begin();
	MutCounter::iterator mit_end = cnt.end();
	for (; mit != mit_end; ++mit) {
		SelMap::iterator sit = m_selFactory.find(mit->first);
		if (sit == m_selFactory.end()) {
			SelCoef sf = getFitnessValue(mit->first);
			if (mit->second == 1 && !chrX)
				s += sf.first * sf.second;
			else
				s += sf.first;
		} else {
			if (mit->second == 1 && !chrX)
				s += sit->second.first * sit->second.second;
			else
				s += sit->second.first;
		}
	}
	return 1 - s > 0 ? 1 - s : 0;
}


double MutSpaceSelector::randomSelExpFitnessExt(GenoIterator it, GenoIterator it_end, bool chrX) const
{
	MutCounter cnt;

	for (; it != it_end; ++it) {
		if (*it == 0u)
			continue;
		MutCounter::iterator mit = cnt.find(*it);
		if (mit == cnt.end())
			cnt[*it] = 1;
		else
			++mit->second;
	}

	double s = 0;
	MutCounter::iterator mit = cnt.begin();
	MutCounter::iterator mit_end = cnt.end();
	for (; mit != mit_end; ++mit) {
		SelMap::iterator sit = m_selFactory.find(mit->first);
		if (sit == m_selFactory.end()) {
			SelCoef sf = getFitnessValue(mit->first);
			if (mit->second == 1 && !chrX)
				s += sf.first * sf.second;
			else
				s += sf.first;
		} else {
			if (mit->second == 1 && !chrX)
				s += sit->second.first * sit->second.second;
			else
				s += sit->second.first;
		}
	}
	return exp(-s);
}


size_t MutSpaceMutator::locateVacantLocus(Population & /* pop */, size_t beg, size_t end, std::set<size_t> & mutants) const
{
	size_t loc = getRNG().randInt(static_cast<ULONG>(end - beg)) + beg;

	std::set<size_t>::iterator it = std::find(mutants.begin(), mutants.end(), loc);

	if (it == mutants.end())
		return loc;
	// look forward and backward
	size_t loc1 = loc + 1;
	std::set<size_t>::iterator it1(it);
	++it1;
	for (; it1 != mutants.end() && loc1 != end; ++it1, ++loc1) {
		if (*it1 != loc1)
			return loc1;
	}
	size_t loc2 = loc - 1;
	std::set<size_t>::reverse_iterator it2(it);
	--it2;
	for (; it2 != mutants.rend() && loc2 != beg; --it2, --loc2) {
		if (*it2 != loc2)
			return loc2;
	}
	// still cannot find
	return 0;
}


bool MutSpaceMutator::apply(Population & pop) const
{
	const matrixi & ranges = m_ranges.elems();
	vectoru width(ranges.size());

	width[0] = ranges[0][1] - ranges[0][0];
	for (size_t i = 1; i < width.size(); ++i)
		width[i] = ranges[i][1] - ranges[i][0] + width[i - 1];

	size_t ploidyWidth = width.back();
	size_t indWidth = pop.ploidy() * ploidyWidth;

	ostream * out = NULL;
	if (!noOutput())
		out = &getOstream(pop.dict());

	// build a set of existing mutants
	std::set<size_t> mutants;
	bool saturated = false;

	subPopList subPops = applicableSubPops(pop);
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	bool chrX = pop.chromType(0) == CHROMOSOME_X;
	for (; sp != spEnd; ++sp) {
		DBG_FAILIF(sp->isVirtual(), ValueError, "This operator does not support virtual subpopulation.");
		for (size_t indIndex = 0; indIndex < pop.subPopSize(sp->subPop()); ++indIndex) {
			size_t loc = 0;
			while (true) {
				// using a geometric distribution to determine mutants
				loc += getRNG().randGeometric(m_rate);
				if (loc > indWidth)
					break;
				Individual & ind = pop.individual(indIndex);
				size_t p = (loc - 1) / ploidyWidth;
				// chromosome and position on chromosome?
				size_t mutLoc = (loc - 1) - p * ploidyWidth;
				// handle chromosome X
				if (p == 1 && chrX && ind.sex() == MALE) {
					if (out)
						(*out) << pop.gen() << '\t' << mutLoc << '\t' << indIndex << "\t4\n";
					continue;
				}
				size_t ch = 0;
				for (size_t reg = 0; reg < width.size(); ++reg) {
					if (mutLoc < width[reg]) {
						ch = reg;
						break;
					}
				}
				mutLoc += ranges[ch][0];
				if (ch > 0)
					mutLoc -= width[ch - 1];

				if (m_model == 2) {
					// under an infinite-site model
					if (saturated) {
						if (out)
							(*out)	<< pop.gen() << '\t' << mutLoc << '\t' << indIndex
							        << "\t3\n";
						continue;
					}
					bool ok = false;
					// if the first time
					if (mutants.empty()) {
						// first try our luck...
						ok = find(pop.genoBegin(false), pop.genoEnd(false), TO_ALLELE(mutLoc)) == pop.genoEnd(false);
						if (!ok) {
							std::set<size_t> existing(pop.genoBegin(false), pop.genoEnd(false));
							mutants.swap(existing);
							mutants.erase(0);
							saturated = mutants.size() == ploidyWidth;
							if (saturated)
								cerr << "Failed to introduce new mutants at generation " << pop.gen() << " because all loci have existing mutants." << endl;
						}
					}
					if (!ok && mutants.find(mutLoc) != mutants.end()) {
						size_t newLoc = locateVacantLocus(pop, ranges[ch][0], ranges[ch][1], mutants);
						// nothing is found
						if (out)
							(*out)	<< pop.gen() << '\t' << mutLoc << '\t' << indIndex
							        << (newLoc == 0 ? "\t3\n" : "\t2\n");
						if (newLoc != 0)
							mutLoc = newLoc;
						else {
							cerr << "Failed to introduce a new mutant at generation " << pop.gen() << " because all loci have existing mutants." << endl;
							// ignore this mutation, and subsequent mutations...
							saturated = true;
							continue;
						}
						// if there is no existing mutant, new mutant is allowed
					}
					mutants.insert(mutLoc);
				}
				GenoIterator geno = ind.genoBegin(p, ch);
				size_t nLoci = pop.numLoci(ch);
				if (*(geno + nLoci - 1) != 0u) {
					// if the number of mutants at this individual exceeds reserved numbers
					DBG_DO(DBG_MUTATOR, cerr << "Adding 10 loci to region " << ch << endl);
					vectorf added(10);
					for (size_t j = 0; j < 10; ++j)
						added[j] = static_cast<double>(nLoci + j + 1);
					vectoru addedChrom(10, ch);
					pop.addLoci(addedChrom, added);
					// individual might be shifted...
					ind = pop.individual(indIndex);
					geno = ind.genoBegin(p, ch);
					nLoci += 10;
				}
				// find the first non-zero location
				for (size_t j = 0; j < nLoci; ++j) {
					if (*(geno + j) == 0u) {
						// record mutation here
						DBG_FAILIF(mutLoc >= ModuleMaxAllele, RuntimeError,
							"Location can not be saved because it exceed max allowed allele.");
						*(geno + j) = TO_ALLELE(mutLoc);
						if (out)
							(*out) << pop.gen() << '\t' << mutLoc << '\t' << indIndex << "\t0\n";
						break;
					} else if (static_cast<size_t>(*(geno + j)) == mutLoc) {
						// back mutation
						//  from A b c d 0
						//  to   d b c d 0
						//  to   d b c 0 0
						for (size_t k = j + 1; k < nLoci; ++k)
							if (*(geno + k) == 0u) {
								*(geno + j) = *(geno + k - 1);
								*(geno + k - 1) = 0;
								if (out)
									(*out) << pop.gen() << '\t' << mutLoc << '\t' << indIndex << "\t1\n";
								break;
							}
						DBG_DO(DBG_MUTATOR, cerr << "Back mutation happens at generation " << pop.gen() << " on individual " << indIndex << endl);
						break;
					}
				}
			}   // while
		}       // each individual
	}           // each subpopulation
	if (out)
		closeOstream();
	return true;
}


void MutSpaceRecombinator::transmitGenotype0(Population & pop, Population & offPop, const Individual & parent,
                                             size_t offIndex, int ploidy) const
{
	size_t nCh = parent.numChrom();

	// chromosome X is not passed to male offspring

	// count duplicates...
	for (size_t ch = 0; ch < parent.numChrom(); ++ch) {
		MutCounter cnt;
		vectoru alleles;
		alleles.reserve(parent.numLoci(ch));
		if (nCh == 1) {
			// this is faster... for a most common case
			GenoIterator it = parent.genoBegin();
			GenoIterator it_end = parent.genoEnd();
			for (; it != it_end; ++it) {
				if (*it == 0u)
					continue;
				MutCounter::iterator mit = cnt.find(*it);
				if (mit == cnt.end())
					cnt[*it] = 1;
				else
					++mit->second;
			}
		} else {
			GenoIterator it = parent.genoBegin(0, ch);
			GenoIterator it_end = parent.genoEnd(0, ch);
			for (; it != it_end; ++it) {
				if (*it == 0u)
					break;
				MutCounter::iterator mit = cnt.find(*it);
				if (mit == cnt.end())
					cnt[*it] = 1;
				else
					++mit->second;
			}
			it = parent.genoBegin(1, ch);
			it_end = parent.genoEnd(1, ch);
			for (; it != it_end; ++it) {
				if (*it == 0u)
					break;
				MutCounter::iterator mit = cnt.find(*it);
				if (mit == cnt.end())
					cnt[*it] = 1;
				else
					++mit->second;
			}
		}
		// no valid allele
		if (cnt.empty()) {
			GenoIterator it = offPop.individual(offIndex).genoBegin(ploidy, ch);
			GenoIterator it_end = offPop.individual(offIndex).genoEnd(ploidy, ch);
			std::fill(it, it_end, 0);
			continue;
		}
		// keep 1 count with probability 0.5, keep 2 count with probability 1
		MutCounter::iterator mit = cnt.begin();
		MutCounter::iterator mit_end = cnt.end();
		for (; mit != mit_end; ++mit) {
			if (mit->second == 2 || getRNG().randBit())
				alleles.push_back(mit->first);
		}
		// not enough size
		if (alleles.size() + 1 > offPop.numLoci(ch)) {
			DBG_DO(DBG_TRANSMITTER, cerr << "Extending size of chromosome " << ch <<
				" to " << alleles.size() + 2 << endl);
			size_t sz = alleles.size() - offPop.numLoci(ch) + 2;
			vectorf added(sz);
			for (size_t j = 0; j < sz; ++j)
				added[j] = static_cast<double>(offPop.numLoci(ch) + j + 1);
			vectoru addedChrom(sz, ch);
			offPop.addLoci(addedChrom, added);
			pop.addLoci(addedChrom, added);
		}
		//
		GenoIterator it = offPop.individual(offIndex).genoBegin(ploidy, ch);
		GenoIterator it_end = offPop.individual(offIndex).genoEnd(ploidy, ch);
		for (size_t i = 0; i < alleles.size(); ++i, ++it)
			*it = TO_ALLELE(alleles[i]);
		// fill the rest with 0.
		std::fill(it, it_end, 0);
	}
}


void MutSpaceRecombinator::transmitGenotype1(Population & pop, Population & offPop, const Individual & parent,
                                             size_t offIndex, int ploidy) const
{
	const matrixi & ranges = m_ranges.elems();

	for (size_t ch = 0; ch < parent.numChrom(); ++ch) {
		size_t width = ranges[ch][1] - ranges[ch][0];
		size_t beg = 0;
		size_t end = getRNG().randGeometric(m_rate);
		int p = getRNG().randBit() ? 0 : 1;
		// no recombination
		if (end >= width) {
			copyChromosome(parent, p, offPop.individual(offIndex), ploidy, ch);
			continue;
		}
		// we are in trouble... get some properties of alleles to reduce comparison
		vectoru alleles;
		size_t minAllele[2];
		size_t maxAllele[2];
		size_t cnt[2];
		cnt[0] = 0;
		cnt[1] = 0;
		minAllele[0] = ranges[ch][1];
		minAllele[1] = ranges[ch][1];
		maxAllele[0] = ranges[ch][0];
		maxAllele[1] = ranges[ch][0];
		GenoIterator it = parent.genoBegin(0, ch);
		GenoIterator it_end = parent.genoEnd(0, ch);
		for (; it != it_end; ++it) {
			if (*it == 0u)
				break;
			++cnt[0];
			if (*it < minAllele[0])
				minAllele[0] = *it;
			if (*it > maxAllele[0])
				maxAllele[0] = *it;
		}
		it = parent.genoBegin(1, ch);
		it_end = parent.genoEnd(1, ch);
		for (; it != it_end; ++it) {
			if (*it == 0u)
				break;
			++cnt[1];
			if (*it < minAllele[1])
				minAllele[1] = *it;
			if (*it > maxAllele[1])
				maxAllele[1] = *it;
		}
		minAllele[0] -= ranges[ch][0];
		minAllele[1] -= ranges[ch][0];
		maxAllele[0] -= ranges[ch][0];
		maxAllele[1] -= ranges[ch][0];
		do {
			// copy piece
			// this algorithm is NOT efficient, but we count the rareness of recombination. :-)
			if (cnt[p] > 0 && end >= minAllele[p] && beg <= maxAllele[p]) {
				it = parent.genoBegin(p, ch);
				it_end = parent.genoEnd(p, ch);
				for (; it != it_end; ++it) {
					if (*it == 0u)
						break;
					if (*it >= beg + ranges[ch][0] && *it < end + ranges[ch][0]) {
						alleles.push_back(*it);
						--cnt[p];
					}
				}
			}
			// change ploidy
			p = (p + 1) % 2;
			// next step
			beg = end;
			end += getRNG().randGeometric(m_rate);
		} while (end < width && (cnt[0] > 0 || cnt[1] > 0));
		// last piece
		if (cnt[0] > 0 || cnt[1] > 0) {
			it = parent.genoBegin(p, ch);
			it_end = parent.genoEnd(p, ch);
			for (; it != it_end; ++it) {
				if (*it >= beg + static_cast<size_t>(ranges[ch][0]) &&
				    *it < static_cast<size_t>(ranges[ch][1]))
					alleles.push_back(*it);
			}
		}
		// set alleles
		// not enough size
		if (alleles.size() + 1 > offPop.numLoci(ch)) {
			DBG_DO(DBG_TRANSMITTER, cerr << "Extending size of chromosome " << ch <<
				" to " << alleles.size() + 2 << endl);
			size_t sz = alleles.size() - offPop.numLoci(ch) + 2;
			vectorf added(sz);
			for (size_t j = 0; j < sz; ++j)
				added[j] = static_cast<double>(offPop.numLoci(ch) + j + 1);
			vectoru addedChrom(sz, ch);
			offPop.addLoci(addedChrom, added);
			pop.addLoci(addedChrom, added);
		}
		//
		it = offPop.individual(offIndex).genoBegin(ploidy, ch);
		it_end = offPop.individual(offIndex).genoEnd(ploidy, ch);
		for (size_t i = 0; i < alleles.size(); ++i, ++it)
			*it = TO_ALLELE(alleles[i]);
		// fill the rest with 0.
		std::fill(it, it_end, 0);
	}
}


bool MutSpaceRecombinator::applyDuringMating(Population & pop, Population & offPop,
                                             RawIndIterator offspring,
                                             Individual * dad, Individual * mom) const
{
	// if offspring does not belong to subPops, do nothing, but does not fail.
	if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
		return true;

	initializeIfNeeded(*offspring);

	if (pop.chromType(0) == CHROMOSOME_X) {
		// for mom
		if (m_rate == 0)
			copyChromosome(*mom, getRNG().randBit(), *offspring, 0, 0);
		else if (m_rate == 0.5)
			transmitGenotype0(pop, offPop, *mom, offspring - offPop.rawIndBegin(), 0);
		else
			transmitGenotype1(pop, offPop, *mom, offspring - offPop.rawIndBegin(), 0);

		// for dad, pass X to daughter
		if (offspring->sex() == MALE)
			return true;
		else
			copyChromosome(*dad, 0, *offspring, 1, 0);
		return true;
	}

	// standard genotype transmitter
	if (m_rate == 0) {
		for (int ch = 0; static_cast<size_t>(ch) < pop.numChrom(); ++ch) {
			copyChromosome(*mom, getRNG().randBit(), *offspring, 0, ch);
			copyChromosome(*dad, getRNG().randBit(), *offspring, 1, ch);
		}
	} else if (m_rate == 0.5) {
		transmitGenotype0(pop, offPop, *mom, offspring - offPop.rawIndBegin(), 0);
		transmitGenotype0(pop, offPop, *dad, offspring - offPop.rawIndBegin(), 1);
	} else {
		transmitGenotype1(pop, offPop, *mom, offspring - offPop.rawIndBegin(), 0);
		transmitGenotype1(pop, offPop, *dad, offspring - offPop.rawIndBegin(), 1);
	}
	return true;
}


#endif
}
}
