/**
 *  $File: Migrator.cpp $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
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

#include "virtualSubPop.h"
#include "migrator.h"

namespace simuPOP {

Migrator::Migrator(const matrix & rate, int mode, const uintList & toSubPops,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops, const stringList & infoFields)
	: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
	m_rate(rate), m_mode(mode), m_to(toSubPops)
{
	DBG_FAILIF(!subPops.empty() && subPops.size() != rate.size(),
		ValueError, "Length of param fromSubPop must match rows of rate matrix.");

	DBG_FAILIF(!m_to.elems().empty() && m_to.elems().size() != rate[0].size(),
		ValueError, "Length of param toSubPop must match columns of rate matrix.");
}


string Migrator::describe(bool format) const
{
	return "<simuPOP.Migrator>";
}


bool Migrator::apply(Population & pop) const
{
	// set info of individual
	UINT info = pop.infoIdx(infoField(0));

	subPopList fromSubPops = applicableSubPops(pop);
	DBG_FAILIF(m_mode != BY_IND_INFO && fromSubPops.size() != m_rate.size(),
		ValueError, "Number of 'from' subpopulations should match number of rows of migration rate matrix.");

	vectorinfo oldInfo;

	if (m_mode == BY_IND_INFO && !fromSubPops.empty()) {
		oldInfo.resize(pop.popSize());
		for (size_t i = 0; i < pop.popSize(); ++i)
			oldInfo[i] = pop.individual(i).info(info);
	}

	if (m_mode != BY_IND_INFO || !fromSubPops.empty()) {
		for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
			RawIndIterator it = pop.rawIndBegin(sp);
			RawIndIterator it_end = pop.rawIndEnd(sp);
			for (; it != it_end; ++it)
				it->setInfo(sp, info);
		}
	}

	DBG_FAILIF(pop.hasActivatedVirtualSubPop(), RuntimeError,
		"Migration can not be applied to activated virtual subpopulations");

	// This one does not support ALL_AVAIL
	vectoru toSubPops = m_to.elems();
	if (m_to.allAvail())
		for (UINT i = 0; i < pop.numSubPop(); ++i)
			toSubPops.push_back(i);

	// real migration matrix might change from population to population because
	// of different number of subpopulations, and toSubPops can be ALL_AVAIL, and
	// then does not have to match subPops.
	matrix migrationRate = m_rate;
	if (m_mode != BY_IND_INFO) {
		UINT szFrom = migrationRate.size();
		UINT szTo = migrationRate[0].size();

		// check parameters
		for (UINT i = 0; i < szFrom; ++i) {
			DBG_FAILIF(migrationRate[i].size() != szTo, ValueError,
				"Expecting a matrix of migration rate.");

			for (size_t j = 0; j < szTo; ++j) {
				DBG_FAILIF(fcmp_lt(migrationRate[i][j], 0.), ValueError,
					"Migration rate should be positive.");
				DBG_FAILIF(m_mode != BY_COUNTS && fcmp_gt(migrationRate[i][j], 1.), ValueError,
					"Migration rate should be in the range of [0,1]");
			}
		}

		// set r[i][i]--- may need to extend rate (to add i->i)
		if (m_mode == BY_PROBABILITY || m_mode == BY_PROPORTION) {
			for (UINT i = 0; i < szFrom; i++) {               // from
				// look for from=to cell.
				UINT spFrom = fromSubPops[i].subPop();
				double sum = accumulate(migrationRate[i].begin(), migrationRate[i].end(), 0.0);
				//
				vectoru::const_iterator spTo = find(toSubPops.begin(), toSubPops.end(), spFrom);
				if (spTo == toSubPops.end()) {                        // if no to, only check if sum <= 1
					if (fcmp_gt(sum, 1.0))
						throw ValueError("Sum of migrate rate from one subPop should <= 1");
					// adding i->i item
					migrationRate[i].push_back(1.0 - sum);
				} else {                                                          // if not, set r[i][i]
					double & self = migrationRate[i][ spTo - toSubPops.begin() ];
					sum -= self;
					if (fcmp_gt(sum, 1.0))
						throw ValueError("Sum of migrate rate from one subPop should <= 1");
					// reset to-my-self probability/proportion
					self = 1.0 - sum;
				}
			}
		}
	}


	for (UINT from = 0, fromEnd = fromSubPops.size(); from < fromEnd; ++from) {
		UINT spFrom = fromSubPops[from].subPop();
		// rateSize might be toSize + 1, the last one is from->from
		UINT toSize = toSubPops.size();
		UINT toIndex;

		// fromSubPops out of range....
		DBG_FAILIF(spFrom >= pop.numSubPop(), IndexError,
			"Subpopulation index " + toStr(spFrom) + " out of range");

		ULONG spSize = pop.subPopSize(fromSubPops[from]);

		if (fromSubPops[from].isVirtual())
			pop.activateVirtualSubPop(fromSubPops[from]);

		if (m_mode == BY_IND_INFO) {
			// restore information fields set by user so that other individuals
			// can stay at their original subpopulation.
			if (!oldInfo.empty()) {
				for (IndIterator ind = pop.indIterator(spFrom); ind.valid(); ++ind)
					ind->setInfo(oldInfo[&*ind - &*pop.rawIndBegin()], info);
			}
		} else if (m_mode == BY_PROBABILITY) {
			Weightedsampler ws(getRNG(), migrationRate[from]);

			// for each individual, migrate according to migration probability
			for (IndIterator ind = pop.indIterator(spFrom); ind.valid(); ++ind) {
				//toIndex = getRNG().randIntByFreq( rateSize, &migrationRate[from][0] ) ;
				toIndex = ws.get();

				DBG_ASSERT(toIndex < migrationRate[from].size(), ValueError,
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
			vectoru toNum(toSize);
			if (m_mode == BY_PROPORTION) {
				// in case that to sub is not in from sub, the last added
				// element is not used. sum of toNum is not spSize.
				for (UINT i = 0; i < toSize; ++i)
					toNum[i] = static_cast<ULONG>(spSize * migrationRate[from][i]);
			} else {                                                                      // by count
				for (UINT i = 0; i < toSize; ++i)
					toNum[i] = static_cast<ULONG>(migrationRate[from][i]);
			}
			// create a vector and assign indexes, then random shuffle
			// and assign info
			vectoru toIndices(spSize);
			UINT k = 0;
			for (UINT i = 0; i < toSize && k < spSize; ++i)
				for (UINT j = 0; j < toNum[i] && k < spSize; ++j)
					toIndices[k++] = toSubPops[i];

			// the rest of individuals will stay in their original subpopulation.
			while (k < spSize)
				toIndices[k++] = spFrom;

			getRNG().randomShuffle(toIndices.begin(), toIndices.end());
			IndIterator ind = pop.indIterator(spFrom);
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
		vectoru split(oldNumSubPop - pop.numSubPop() + 1, 0);
		split[0] = pop.subPopSize(pop.numSubPop() - 1);
		pop.splitSubPop(pop.numSubPop() - 1, split);
	}
	DBG_ASSERT(pop.numSubPop() >= oldNumSubPop, RuntimeError,
		"Migrator should not decrease number of subpopulations.");

	return true;
}


struct compareVSP
{
	int operator()(const vspID & v1, const vspID & v2)
	{
		return v1.subPop() > v2.subPop();
	}


};

bool SplitSubPops::apply(Population & pop) const
{
	subPopList subPops = applicableSubPops(pop);

	// we have to split from top to bottom subpopulations
	// because split will change subpopulation index
	std::sort(subPops.begin(), subPops.end(), compareVSP());

	for (size_t i = 0; i < subPops.size(); ++i) {
		SubPopID sp = subPops[i].subPop();
		if (pop.subPopSize(sp) == 0)
			continue;

		// randomize indiviudlas
		if (m_randomize) {
			getRNG().randomShuffle(pop.rawIndBegin(sp), pop.rawIndEnd(sp));
			pop.setIndOrdered(false);
		}

		if (!m_subPopSizes.empty())
			pop.splitSubPop(sp, m_subPopSizes, m_names);
		else if (!m_proportions.empty()) {
			vectoru sizes;
			for (size_t i = 0; i < m_proportions.size(); ++i)
				sizes.push_back(static_cast<ULONG>(pop.subPopSize(sp) * m_proportions[i]));
			// floating point problem.
			sizes[m_proportions.size() - 1] += pop.subPopSize(sp) - accumulate(sizes.begin(), sizes.end(), 0LU);
			pop.splitSubPop(sp, sizes, m_names);
		} else {
			// using an information field
			UINT idx = pop.infoIdx(infoField(0));
			std::sort(pop.rawIndBegin(sp), pop.rawIndEnd(sp), indCompare(idx));
			ConstRawIndIterator it = pop.rawIndBegin(sp);
			ConstRawIndIterator it_end = pop.rawIndEnd(sp);
			vectoru spSize(1, 1);
			double value = it->info(idx);
			for (++it; it != it_end; ++it) {
				if (fcmp_ne(it->info(idx), value)) {
					spSize.push_back(1);
					value = it->info(idx);
				} else
					++spSize.back();
			}
			pop.splitSubPop(sp, spSize, m_names);
		}
	}
	return true;
}


bool MergeSubPops::apply(Population & pop) const
{
	const subPopList sp = applicableSubPops(pop);

	vectoru subPops(sp.size());

	for (size_t i = 0; i < sp.size(); ++i)
		subPops[i] = sp[i].subPop();
	pop.mergeSubPops(subPops, m_name);
	return true;
}


bool ResizeSubPops::apply(Population & pop) const
{
	vectoru newSizes = pop.subPopSizes();

	const subPopList subPops = applicableSubPops(pop);

	DBG_FAILIF(subPops.empty() && m_sizes.size() != pop.numSubPop()
		&& m_proportions.size() != pop.numSubPop(), ValueError,
		"Please specify new subpopulation size for all subpopulations.");

	for (size_t i = 0; i < subPops.size(); ++i) {
		DBG_FAILIF(static_cast<UINT>(subPops[i].subPop()) >= pop.numSubPop(), IndexError,
			"Subpopulation index " + toStr(subPops[i].subPop()) + " out of range of 0 ~ "
			+ toStr(pop.numSubPop() - 1));
		if (m_sizes.empty())
			newSizes[subPops[i].subPop()] = static_cast<ULONG>(newSizes[subPops[i].subPop()] * m_proportions[i]);
		else
			newSizes[subPops[i].subPop()] = m_sizes[i];
	}
	DBG_DO(DBG_MIGRATOR, cerr << "Resize subpopulations to size " << newSizes << endl);
	pop.resize(newSizes, m_propagate);
	return true;
}


}
