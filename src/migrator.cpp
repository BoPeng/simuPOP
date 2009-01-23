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

#include "virtualSubPop.h"
#include "migrator.h"

namespace simuPOP {

migrator::migrator(const matrix & rate, int mode, const uintList & toSubPops,
	int stage, int begin, int end, int step, const intList & at,
	const repList & rep, const subPopList & subPops, const vectorstr & infoFields)
	: baseOperator("", stage, begin, end, step, at, rep, subPops, infoFields),
	m_rate(rate), m_mode(mode), m_to(toSubPops.elems())
{
	DBG_FAILIF(!subPops.empty() && subPops.size() != rate.size(),
		ValueError, "Length of param fromSubPop must match rows of rate matrix.");

	DBG_FAILIF(!m_to.empty() && m_to.size() != rate[0].size(),
		ValueError, "Length of param toSubPop must match columns of rate matrix.");
}


void migrator::setRates(int mode, const subPopList & fromSubPops, const vectorlu & toSubPops)
{
	if (mode == ByIndInfo)
		return;

	UINT szFrom = m_rate.size();
	UINT szTo = m_rate[0].size();

	m_mode = mode;

	if (m_mode != ByProbability && m_mode != ByProportion && m_mode != ByCounts)
		throw ValueError("Migration mode can only be ByProbability, ByProportion or ByCounts");

	// check parameters
	for (UINT i = 0; i < szFrom; ++i) {
		DBG_FAILIF(m_rate[i].size() != szTo, ValueError,
			"Expecting a matrix of migration rate.");

		for (size_t j = 0; j < szTo; ++j)
			DBG_FAILIF(fcmp_lt(m_rate[i][j], 0.), ValueError,
				"Migration rate should be positive.");
			DBG_FAILIF(m_mode != ByCounts && fcmp_gt(m_rate[i][j], 1.), ValueError,
				"Migration rate should be in the range of [0,1]");
	}

	// set r[i][i]--- may need to extend rate (to add i->i)
	if (m_mode == ByProbability || m_mode == ByProportion) {
		for (UINT i = 0; i < szFrom; i++) {               // from
			// look for from=to cell.
			UINT spFrom = fromSubPops[i].subPop();
			double sum = accumulate(m_rate[i].begin(), m_rate[i].end(), 0.0);
			//
			vectorlu::const_iterator spTo = find(toSubPops.begin(), toSubPops.end(), spFrom);
			if (spTo == toSubPops.end() ) {                        // if no to, only check if sum <= 1
				if (fcmp_gt(sum, 1.0) )
					throw ValueError("Sum of migrate rate from one subPop should <= 1");
				// adding i->i item
				m_rate[i].push_back(1.0 - sum);
			} else {                                                          // if not, set r[i][i]
				double & self = m_rate[i][ spTo - toSubPops.begin() ];
				sum -= self;
				if (fcmp_gt(sum, 1.0) )
					throw ValueError("Sum of migrate rate from one subPop should <= 1");
				// reset to-my-self probability/proportion
				self = 1.0 - sum;
			}
		}
	}
}


bool migrator::apply(population & pop)
{
	// set info of individual
	UINT info = pop.infoIdx(infoField(0));

	subPopList fromSubPops = applicableSubPops();
	
	vectorinfo oldInfo;
	if (m_mode == ByIndInfo && !fromSubPops.empty()) {
		oldInfo.resize(pop.popSize());
		for (size_t i = 0; i < pop.popSize(); ++i)
			oldInfo[i] = pop.ind(i).info(info);
	}

	if (m_mode != ByIndInfo || !fromSubPops.empty()) {
		for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
			RawIndIterator it = pop.rawIndBegin(sp);
			RawIndIterator it_end = pop.rawIndEnd(sp);
			for (; it != it_end; ++it)
				it->setInfo(sp, info);
		}
	}

	DBG_FAILIF(pop.hasActivatedVirtualSubPop(), RuntimeError,
		"Migration can not be applied to activated virtual subpopulations");

	if (fromSubPops.empty() )
		for (UINT i = 0; i < pop.numSubPop(); ++i)
			fromSubPops.push_back(vspID(i));

	vectorlu toSubPops = m_to;
	if (m_to.empty() )
		for (UINT i = 0; i < pop.numSubPop(); ++i)
			toSubPops.push_back(i);

	setRates(m_mode, fromSubPops, toSubPops);

	for (UINT from = 0, fromEnd = fromSubPops.size(); from < fromEnd; ++from) {
		UINT spFrom = fromSubPops[from].subPop();
		// rateSize might be toSize + 1, the last one is from->from
		UINT toSize = toSubPops.size();
		UINT toIndex;

		// fromSubPops out of range....
		DBG_FAILIF(spFrom >= pop.numSubPop(), IndexError,
			"Subpopulation index " + toStr(spFrom) + " out of range");

		if (fromSubPops[from].isVirtual())
			pop.activateVirtualSubPop(fromSubPops[from]);

		// if subpopulation spFrom is activated, subPopSize can
		// get the correct size.
		ULONG spSize = pop.subPopSize(spFrom);

		if (m_mode == ByIndInfo) {
			// restore information fields set by user so that other individuals
			// can stay at their original subpopulation.
			if (!oldInfo.empty()) {
				for (IndIterator ind = pop.indBegin(spFrom); ind.valid();  ++ind)
					ind->setInfo(oldInfo[&*ind - &*pop.rawIndBegin()], info);
			}
		} else if (m_mode == ByProbability) {
			weightedSampler ws(rng(), m_rate[from]);

			// for each individual, migrate according to migration probability
			for (IndIterator ind = pop.indBegin(spFrom); ind.valid();  ++ind) {
				//toIndex = rng().randIntByFreq( rateSize, &m_rate[from][0] ) ;
				toIndex = ws.get();

				DBG_ASSERT(toIndex < m_rate[from].size(), ValueError,
					"Return index out of range.");

				// rateSize == toSize (no i->i addition)
				//     toIndex < toSize
				// rateSize = toSize + 1, ignore i->1 (last one)
				//     toIndex < toSize
				if (toIndex < toSize && toSubPops[toIndex] != spFrom)
					ind->setInfo(toSubPops[toIndex], info);
			}
		} else {
			// 2nd, or 3rd method
			// first find out how many people will move to other subPop
			// then randomly assign individuals to move
			vectorlu toNum(toSize);
			if (m_mode == ByProportion) {
				// in case that to sub is not in from sub, the last added
				// element is not used. sum of toNum is not spSize.
				for (UINT i = 0; i < toSize; ++i)
					toNum[i] = static_cast<ULONG>(spSize * m_rate[from][i]);
			} else {                                                                      // by count
				for (UINT i = 0; i < toSize; ++i)
					toNum[i] = static_cast<ULONG>(m_rate[from][i]);
			}
			// create a vector and assign indexes, then random shuffle
			// and assign info
			vectorlu toIndices(spSize);
			UINT k = 0;
			for (UINT i = 0; i < toSize && k < spSize; ++i)
				for (UINT j = 0; j < toNum[i] && k < spSize; ++j)
					toIndices[k++] = toSubPops[i];

			// the rest of individuals will stay in their original subpopulation.
			while (k < spSize)
				toIndices[k++] = spFrom;

			random_shuffle(toIndices.begin(), toIndices.end());
			IndIterator ind = pop.indBegin(spFrom);
			// set info
			for (UINT i = 0; ind.valid(); ++i, ++ind)
				// The previous migration_to value, if set by a previous vsp, will be overridden.
				ind->setInfo(static_cast<InfoType>(toIndices[i]), info);
		}
		if (fromSubPops[from].isVirtual())
			pop.deactivateVirtualSubPop(spFrom);
	}   // for all subPop.

	// do migration.
	UINT oldNumSubPop = pop.numSubPop();
	pop.setSubPopByIndInfo(infoField(0));
	// try to keep number of subpopulations.
	if (pop.numSubPop() < oldNumSubPop) {
		vectorf split(oldNumSubPop - pop.numSubPop() + 1, 0);
		split[0] = 1.;
		pop.splitSubPop(pop.numSubPop() - 1, split);
	}
	DBG_ASSERT(pop.numSubPop() >= oldNumSubPop, RuntimeError,
		"Migrator should not decrease number of subpopulations.");

	return true;
}


bool splitSubPop::apply(population & pop)
{
	UINT info = pop.infoIdx(infoField(0));

	// randomize indiviudlas
	if (m_randomize) {
		// random shuffle individuals
		for (ULONG it = 0; it < pop.subPopSize(m_which); ++it)
			pop.ind(it, m_which).setInfo(static_cast<SubPopID>(rng().randInt(MaxSubPopID)), info);
		std::sort(pop.indBegin(m_which), pop.indEnd(m_which), indCompare(info));
		// not actully required since spliSubPop will do this.
		// this is to remind myself this step is important.
		pop.setIndOrdered(false);
	}
	if (!m_subPopSizes.empty()) {
		vectorf sizes;
		for (size_t i = 0; i < m_subPopSizes.size(); ++i)
			sizes.push_back(m_subPopSizes[i]);
		pop.splitSubPop(m_which, sizes);
	} else
		pop.splitSubPop(m_which, m_proportions);
	return true;
}


bool resizeSubPops::apply(population & pop)
{
	vectorlu newSizes = pop.subPopSizes();

	subPopList subPops = applicableSubPops();

	DBG_FAILIF(subPops.empty() && subPops.size() != pop.numSubPop(), ValueError,
		"Please specify new subpopulation size for all subpopulations.");

	if (subPops.empty()) {
		for (size_t i = 0; i < pop.numSubPop(); ++i)
			newSizes[i] = m_newSizes[i];
	} else {
		for (size_t i = 0; i < subPops.size(); ++i) {
			DBG_FAILIF(subPops[i].subPop() >= pop.numSubPop(), IndexError,
				"Subpopulation index " + toStr(subPops[i].subPop()) + " out of range of 0 ~ "
				+ toStr(pop.numSubPop() - 1));
			newSizes[subPops[i].subPop()] = m_newSizes[i];
		}
	}
	DBG_DO(DBG_MIGRATOR, cout << "Resize subpopulations to size " << newSizes << endl);
	pop.resize(newSizes, m_propagate);
	return true;
}


}
