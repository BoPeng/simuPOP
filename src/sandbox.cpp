
/**
 *  $File: sandbox.cpp $
 *  $LastChangedDate: 2010-06-04 13:29:09 -0700 (Fri, 04 Jun 2010) $
 *  $Rev: 3579 $
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


double InfSitesSelector::indFitness(Individual * ind, ULONG gen) const
{
	if (m_mode == MULTIPLICATIVE)
		return randomSelMulFitness(ind->genoBegin(), ind->genoEnd());
	else if (m_mode == ADDITIVE)
		return randomSelAddFitness(ind->genoBegin(), ind->genoEnd());
	else if (m_mode == EXPONENTIAL)
		return randomSelExpFitness(ind->genoBegin(), ind->genoEnd());
	return 0;
}


bool InfSitesSelector::apply(Population & pop) const
{
	m_newMutants.clear();
	if (!BaseSelector::apply(pop))
		return false;
	// output NEW mutant...
	if (!m_newMutants.empty() && !noOutput()) {
		ostream & out = getOstream(pop.dict());
		vectoru::const_iterator it = m_newMutants.begin();
		vectoru::const_iterator it_end = m_newMutants.end();
		for (; it != it_end; ++it)
			out << *it << '\t' << m_selFactory[*it] << '\n';
		closeOstream();
	}
	return true;
}


double InfSitesSelector::getFitnessValue(int mutant) const
{
	int sz = m_selDist.size();
	double s;

	if (sz == 0)
		// call a function
		s = m_selDist.func() (PyObj_As_Double, "()");
	else if (sz == 2)
		// constant
		s = m_selDist[1];
	else
		// a gamma distribution
		s = -getRNG().randGamma(m_selDist[1], m_selDist[2]);
	m_selFactory[mutant] = s;
	m_newMutants.push_back(mutant);
	return s;
}


double InfSitesSelector::randomSelMulFitness(GenoIterator it, GenoIterator it_end) const
{
	double s = 1;

	for (; it != it_end; ++it) {
		if (*it == 0)
			continue;
		SelMap::iterator sit = m_selFactory.find(*it);
		if (sit == m_selFactory.end())
			s *= 1 + getFitnessValue(static_cast<unsigned int>(*it));
		else
			s *= 1 + sit->second;
	}
	return s;
}


double InfSitesSelector::randomSelAddFitness(GenoIterator it, GenoIterator it_end) const
{
	double s = 0;

	for (; it != it_end; ++it) {
		if (*it == 0)
			continue;
		SelMap::iterator sit = m_selFactory.find(static_cast<unsigned int>(*it));
		if (sit == m_selFactory.end())
			s += getFitnessValue(*it);
		else
			s += sit->second;
	}
	return 1 + s > 0 ? 1 + s : 0;
}


double InfSitesSelector::randomSelExpFitness(GenoIterator it, GenoIterator it_end) const
{
	double s = 0;

	for (; it != it_end; ++it) {
		if (*it == 0)
			continue;
		SelMap::iterator sit = m_selFactory.find(static_cast<unsigned int>(*it));
		if (sit == m_selFactory.end())
			s += getFitnessValue(*it);
		else
			s += sit->second;
	}
	return exp(s);
}


bool InfSitesMutator::apply(Population & pop) const
{
	const matrixi & ranges = m_ranges.elems();
	vectoru width(ranges.size());

	width[0] = ranges[0][1] - ranges[0][0];
	for (size_t i = 1; i < width.size(); ++i)
		width[i] = ranges[i][1] - ranges[i][0] + width[i - 1];

	ULONG ploidyWidth = width.back();
	ULONG indWidth = pop.ploidy() * ploidyWidth;

	ostream * out = NULL;
	if (!noOutput())
		out = &getOstream(pop.dict());

	subPopList subPops = applicableSubPops(pop);
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for (; sp != spEnd; ++sp) {
		DBG_FAILIF(sp->isVirtual(), ValueError, "This operator does not support virtual subpopulation.");
		for (size_t indIndex = 0; indIndex < pop.subPopSize(sp->subPop()); ++indIndex) {
			ULONG loc = 0;
			while (true) {
				// using a geometric distribution to determine mutants
				loc += getRNG().randGeometric(m_rate);
				if (loc > indWidth)
					break;
				Individual & ind = pop.individual(indIndex);
				int p = (loc - 1) / ploidyWidth;
				// chromosome and position on chromosome?
				ULONG mutLoc = (loc - 1) - p * ploidyWidth;
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

				if (out)
					(*out) << pop.gen() << '\t' << mutLoc << '\t' << indIndex << '\n';

				GenoIterator geno = ind.genoBegin(p, ch);
				size_t nLoci = pop.numLoci(ch);
				if (*(geno + nLoci - 1) != 0) {
					// if the number of mutants at this individual exceeds reserved numbers
					DBG_DO(DBG_MUTATOR, cerr << "Adding 10 loci to region " << ch << endl);
					vectorf added(10);
					for (size_t j = 0; j < 10; ++j)
						added[j] = nLoci + j + 1;
					vectoru addedChrom(10, ch);
					pop.addLoci(addedChrom, added);
					// individual might be shifted...
					ind = pop.individual(indIndex);
					geno = ind.genoBegin(p, ch);
					nLoci += 10;
				}
				// find the first non-zero location
				for (size_t j = 0; j < nLoci; ++j) {

					if (*(geno + j) == 0) {
						// record mutation here
						*(geno + j) = mutLoc;
						break;
					} else if (static_cast<ULONG>(*(geno + j)) == mutLoc) {
						// back mutation
						//  from A b c d 0
						//  to   d b c d 0
						//  to   d b c 0 0
						for (size_t k = j + 1; k < nLoci; ++k)
							if (*(geno + k) == 0) {
								*(geno + j) = *(geno + k - 1);
								*(geno + k - 1) = 0;
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


}
